'''A collection of different strategies for iterating over streams
of :class:`~.Scan`-like objects.
'''
import warnings
from collections import deque, defaultdict

from . import ScanBunch


def _null_scan_validator(scan):
    return True


def _null_scan_cacher(scan):
    return None


ITERATION_MODE_GROUPED = "grouped"
ITERATION_MODE_SINGLE = "single"


class _ScanIteratorImplBase(object):
    """Internal base class for scan iteration strategies

    Attributes
    ----------
    iterator : :class:`Iterator`
        An iterator that produces a raw scan which will be packaged by :attr:`scan_packer`
    scan_cacher : :class:`Callable`
        A callable that will cache the scan, or a no-op
    scan_packer : :class:`Callable`
        A callable that will package a raw scan object into a :class:`Scan` object
    scan_validator : :class:`Callable`
        A callable that will check if a raw scan object is valid or not for filtering
        out unwanted entries, or a no-op
    """

    def __init__(self, iterator, scan_packer, scan_validator=None, scan_cacher=None, **kwargs):
        if scan_validator is None:
            scan_validator = _null_scan_validator
        if scan_cacher is None:
            scan_cacher = _null_scan_cacher
        self.iterator = iterator
        self.scan_packer = scan_packer
        self.scan_validator = scan_validator
        self.scan_cacher = scan_cacher
        self._producer = None

    def __iter__(self):
        return self

    def next(self):
        if self._producer is None:
            self._producer = self._make_producer()
        return next(self._producer)

    def __next__(self):
        return self.next()

    def _make_producer(self):
        raise NotImplementedError()

    @property
    def iteration_mode(self):
        """The iteration mode of the strategy, "grouped" or "single"
        """
        raise NotImplementedError()

    @classmethod
    def from_scan_source(cls, iterator, scan_source, **kwargs):
        """Create an iterator strategy for `iterator` from `scan_source`

        Parameters
        ----------
        iterator : :class:`Iterable`
            The iterator over raw scan data to consume
        scan_source : :class:`~.ScanIterator`
            The data extraction wrapper to provide with an iteration strategy
        """
        return cls(iterator, scan_source._make_scan, scan_source._validate, scan_source._cache_scan, **kwargs)


class _SingleScanIteratorImpl(_ScanIteratorImplBase):
    """Iterate over individual scans.

    The default strategy when MS1 scans are missing.
    """

    @property
    def iteration_mode(self):
        return ITERATION_MODE_SINGLE

    def _make_producer(self):
        _make_scan = self.scan_packer
        _validate = self.scan_validator
        _cache_scan = self.scan_cacher
        for scan in self.iterator:
            packed = _make_scan(scan)
            if not _validate(packed):
                continue
            _cache_scan(packed)
            yield packed


class _FakeGroupedScanIteratorImpl(_SingleScanIteratorImpl):
    '''Mimics the interface of :class:`_GroupedScanIteratorImpl` for
    scan sequences which only support single scans, or which do not
    guarantee sequential access to precursor/product collections.
    '''
    @property
    def iteration_mode(self):
        return ITERATION_MODE_GROUPED

    def _make_producer(self):
        generator = super(_FakeGroupedScanIteratorImpl, self)._make_producer()
        for scan in generator:
            if scan.ms_level == 1:
                yield ScanBunch(scan, [])
            else:
                yield ScanBunch(None, [scan])


class _GroupedScanIteratorImpl(_ScanIteratorImplBase):
    """Iterate over related scan bunches.

    The default strategy when MS1 scans are known to be
    present, even if MSn scans are not.
    """

    @property
    def iteration_mode(self):
        return ITERATION_MODE_GROUPED

    def _make_producer(self):
        _make_scan = self.scan_packer
        _validate = self.scan_validator
        _cache_scan = self.scan_cacher

        precursor_scan = None
        product_scans = []

        current_level = 1

        for scan in self.iterator:
            packed = _make_scan(scan)
            if not _validate(packed):
                continue
            _cache_scan(packed)
            if packed.ms_level > 1:
                # inceasing ms level
                if current_level < packed.ms_level:
                    current_level = packed.ms_level
                # decreasing ms level
                elif current_level > packed.ms_level:
                    current_level = packed.ms_level
                product_scans.append(packed)
            elif packed.ms_level == 1:
                if current_level > 1:
                    if precursor_scan is not None:
                        precursor_scan.product_scans = list(product_scans)
                    yield ScanBunch(precursor_scan, product_scans)
                else:
                    if precursor_scan is not None:
                        precursor_scan.product_scans = list(product_scans)
                        yield ScanBunch(precursor_scan, product_scans)
                precursor_scan = packed
                product_scans = []
            else:
                raise ValueError("Could not interpret MS Level %r" % (packed.ms_level,))
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)


class _InterleavedGroupedScanIteratorImpl(_GroupedScanIteratorImpl):
    """Iterate over related scan bunches.

    The default strategy when MS1 scans are known to be
    present, even if MSn scans are not.
    """

    def __init__(self, iterator, scan_packer, scan_validator=None, scan_cacher=None, buffering=5):
        super(_InterleavedGroupedScanIteratorImpl, self).__init__(iterator, scan_packer, scan_validator, scan_cacher)
        if buffering < 2:
            raise ValueError("Interleaved buffering must be greater than 1")
        self.buffering = buffering
        self.ms1_buffer = deque()
        self.product_mapping = defaultdict(list)
        self.passed_first_ms1 = False

    def deque_group(self):
        precursor = self.ms1_buffer.popleft()
        products = self.product_mapping.pop(precursor.id, [])
        precursor.product_scans = products
        if not self.passed_first_ms1 and self.ms1_buffer:
            current_ms1_time = precursor.scan_time
            next_ms1_time = self.ms1_buffer[0].scan_time
            for prec_id, prods in list(self.product_mapping.items()):
                for prod in prods:
                    if current_ms1_time <= prod.scan_time <= next_ms1_time:
                        products.append(prod)
            self.passed_first_ms1 = True
        return ScanBunch(precursor, products)

    def add_product(self, scan):
        pinfo = scan.precursor_information
        if pinfo is None:
            precursor_id = None
        else:
            precursor_id = pinfo.precursor_scan_id
        if precursor_id is None:
            precursor_id = self.ms1_buffer[-1].id
        self.product_mapping[precursor_id].append(scan)

    def add_precursor(self, scan):
        self.ms1_buffer.append(scan)
        return len(self.ms1_buffer) >= self.buffering

    def _make_producer(self):
        _make_scan = self.scan_packer
        _validate = self.scan_validator
        _cache_scan = self.scan_cacher

        current_level = 1

        for scan in self.iterator:
            packed = _make_scan(scan)
            if not _validate(packed):
                continue
            _cache_scan(packed)
            if packed.ms_level > 1:
                # inceasing ms level
                if current_level < packed.ms_level:
                    current_level = packed.ms_level
                # decreasing ms level
                elif current_level > packed.ms_level:
                    current_level = packed.ms_level
                self.add_product(packed)
            elif packed.ms_level == 1:
                do_emit = self.add_precursor(packed)
                if do_emit:
                    yield self.deque_group()
            else:
                raise ValueError("Could not interpret MS Level %r" %
                                 (packed.ms_level,))
        while self.ms1_buffer:
            yield self.deque_group()
        if self.product_mapping:
            warnings.warn("Lingering Product Sets For %r!" % (list(self.product_mapping), ))


class MSEIterator(_GroupedScanIteratorImpl):
    def __init__(self, iterator, scan_packer, low_energy_config, lock_mass_config, scan_validator=None, scan_cacher=None):
        super(MSEIterator, self).__init__(iterator, scan_packer, scan_validator, scan_cacher)
        self.low_energy_config = low_energy_config
        self.lock_mass_config = lock_mass_config

    def _make_producer(self):
        _make_scan = self.scan_packer
        _validate = self.scan_validator
        _cache_scan = self.scan_cacher

        precursor_scan = None
        product_scans = []

        current_level = 'low'

        for scan in self.iterator:
            packed = _make_scan(scan)
            if not _validate(packed):
                continue
            _cache_scan(packed)
            config = packed.acquisition_information[0].scan_configuration
            if config == self.lock_mass_config:
                continue

            if config != self.low_energy_config:
                if current_level == 'low':
                    current_level = 'high'
                # decreasing ms level
                product_scans.append(packed)
            elif config == self.low_energy_config:
                if current_level != 'low':
                    yield ScanBunch(precursor_scan, product_scans)
                else:
                    if precursor_scan is not None:
                        yield ScanBunch(precursor_scan, product_scans)
                precursor_scan = packed
                product_scans = []
            else:
                raise ValueError("Could not interpret MS Level %r" % (packed.ms_level,))
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)
