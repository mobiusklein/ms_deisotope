'''A collection of different strategies for iterating over streams
of :class:`~.Scan`-like objects.
'''
import bisect
import warnings
from collections import deque, defaultdict
from six import PY2

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
        '''Py2 compatible iterator
        '''
        if self._producer is None:
            self._producer = self._make_producer()
        return next(self._producer)

    def __next__(self):
        return self.next() # pylint: disable=not-callable

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
                raise ValueError("Could not interpret MS Level %r" %
                                 (packed.ms_level,))
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)


class GenerationTracker(object):
    _generations_type = list if PY2 else deque

    def __init__(self):
        self.generation_to_id = defaultdict(set)
        self.id_to_generation = dict()
        self.generations = self._generations_type()

    def _add_generation(self, generation):
        bisect.insort_left(self.generations, generation)

    def clear(self):
        self.generation_to_id.clear()
        self.id_to_generation.clear()
        self.generations = self._generations_type()

    def add(self, identifier, generation):
        if generation not in self.generation_to_id:
            self._add_generation(generation)
        self.generation_to_id[generation].add(identifier)
        self.id_to_generation[identifier] = generation

    def remove(self, identifier):
        if identifier not in self.id_to_generation:
            return False
        generation = self.id_to_generation[identifier]
        self.generation_to_id[generation].remove(identifier)
        del self.id_to_generation[identifier]
        if len(self.generation_to_id[generation]) == 0:
            self.generations.remove(generation)
        return True

    def older_than(self, generation):
        result = []
        purged = []
        for gen in self.generations:
            if gen < generation:
                members = self.generation_to_id[gen]
                result.extend(members)
                purged.append(gen)
            else:
                break
        for r in result:
            del self.id_to_generation[r]
        for gen in purged:
            del self.generation_to_id[gen]
            self.generations.remove(gen)
        return result



class _InterleavedGroupedScanIteratorImpl(_GroupedScanIteratorImpl):
    """Iterate over related scan bunches.

    The default strategy when MS1 scans are known to be
    present, even if MSn scans are not.
    """

    def __init__(self, iterator, scan_packer, scan_validator=None, scan_cacher=None, buffering=5):
        super(_InterleavedGroupedScanIteratorImpl, self).__init__(
            iterator, scan_packer, scan_validator, scan_cacher)
        if buffering < 2:
            raise ValueError("Interleaved buffering must be greater than 1")
        self.buffering = buffering
        self.ms1_buffer = deque()
        self.product_mapping = defaultdict(list)
        self.generation_tracker = GenerationTracker()
        self.orphans = []
        self.passed_first_ms1 = False
        self.highest_ms_level = 0
        self.generation = 0

    def pop_precursor(self, precursor_id):
        self.generation_tracker.remove(precursor_id)
        return self.product_mapping.pop(precursor_id, [])

    def deque_group(self, flush_products=False):
        '''Remove the next scan from the MS1 queue, grouped with
        any associated MSn scans.

        Parameters
        ----------
        flush_products : bool
            Whether to flush all the remaining product scans with this
            group.

        Returns
        -------
        ScanBunch
        '''
        precursor = self.ms1_buffer.popleft()
        _empty = []
        products = self.pop_precursor(precursor.id)
        if None in self.product_mapping:
            products += self.product_mapping.pop(None, _empty)

        # Flush out older precursors' products that haven't turned up yet, they
        # probably aren't coming soon.
        for prec_id in self.generation_tracker.older_than(self.generation - self.buffering):
            products.extend(self.product_mapping.pop(prec_id, _empty))

        # Look for MSn for n > 2
        if self.highest_ms_level > 2:
            extra_blocks = [products]
            for _ in range(self.highest_ms_level - 2):
                new_block = []
                for prod in extra_blocks[-1]:
                    new_block.extend(self.pop_precursor(prod.id))
                if not new_block:
                    break
                extra_blocks.append(new_block)
            if len(extra_blocks) > 1:
                products = extra_blocks[0]
                for block in extra_blocks[1:]:
                    products.extend(block)

        if flush_products:
            if self.product_mapping:
                lingering = {s.id for ss in self.product_mapping.values() for s in ss}
                missing = set(self.product_mapping)
                unclaimed_precursors = missing - lingering
                warnings.warn("Lingering Product Sets For %r!" %
                              (sorted(unclaimed_precursors), ))
            for _, value in self.product_mapping.items():
                products += value
            self.product_mapping.clear()
            self.generation_tracker.clear()
        precursor.product_scans = products

        # Collect any MSn spectra which pre-date the first precursor if they are encountered before
        # the first precursor is found.
        if not self.passed_first_ms1 and self.ms1_buffer:
            current_ms1_time = precursor.scan_time
            for prec_id, prods in list(self.product_mapping.items()):
                masked = set()
                for i, prod in enumerate(prods):
                    if prod.scan_time <= current_ms1_time:
                        products.append(prod)
                        masked.add(i)
                # We've only removed some of the products under this precursor, so just
                # remove those products from mapping.
                if len(masked) < len(prods):
                    prods = [v for i, v in enumerate(prods) if i not in masked]
                    self.product_mapping[prec_id] = prods
                else:
                    # Otherwise we must have completely consumed the products of this
                    # precursor, so we need to remove it from the tracking.
                    self.pop_precursor(prec_id)
            self.passed_first_ms1 = True


        self.generation += 1
        return ScanBunch(precursor, products)

    def add_product(self, scan):
        '''Add MSn scan to :attr:`product_mapping` for the associated
        precursor scan ID.

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The scan to track.
        '''
        pinfo = scan.precursor_information
        if pinfo is None:
            precursor_id = None
        else:
            precursor_id = pinfo.precursor_scan_id
        if precursor_id is None:
            try:
                precursor_id = self.ms1_buffer[-1].id
            except IndexError:
                precursor_id = None
        if precursor_id not in self.product_mapping:
            self.generation_tracker.add(precursor_id, self.generation)
        self.product_mapping[precursor_id].append(scan)

    def add_precursor(self, scan):
        '''Add MS1 scan to :attr:`ms1_buffer`

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The scan to track.

        Returns
        -------
        buffer_full : bool
            Whether or not :attr:`ms1_buffer` is full.
        '''
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
                    if current_level > self.highest_ms_level:
                        self.highest_ms_level = current_level

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

        while len(self.ms1_buffer) > 1:
            yield self.deque_group()

        if self.ms1_buffer:
            yield self.deque_group(flush_products=True)


class MSEIterator(_GroupedScanIteratorImpl):
    '''A scan iterator implementation for grouping MS^E spectra according
    to the specified functions.

    Attributes
    ----------
    low_energy_config : int
        The function corresponding to lower energy. These correspond
        to the MS1 equivalent.
    lock_mass_config : int
        The function corresponding to the lockmass. Lockmass scans
        will be skipped.
    '''
    def __init__(self, iterator, scan_packer, low_energy_config, lock_mass_config, scan_validator=None,
                 scan_cacher=None):
        super(MSEIterator, self).__init__(
            iterator, scan_packer, scan_validator, scan_cacher)
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
                raise ValueError("Could not interpret MS Level %r" %
                                 (packed.ms_level,))
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)
