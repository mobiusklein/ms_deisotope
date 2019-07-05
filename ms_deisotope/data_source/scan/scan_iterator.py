'''A collection of different strategies for iterating over streams
of :class:`~.Scan`-like objects.
'''
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

    def __init__(self, iterator, scan_packer, scan_validator=None, scan_cacher=None):
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
    def from_scan_source(cls, iterator, scan_source):
        """Create an iterator strategy for `iterator` from `scan_source`

        Parameters
        ----------
        iterator : :class:`Iterable`
            The iterator over raw scan data to consume
        scan_source : :class:`~.ScanIterator`
            The data extraction wrapper to provide with an iteration strategy
        """
        return cls(iterator, scan_source._make_scan, scan_source._validate, scan_source._cache_scan)



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
