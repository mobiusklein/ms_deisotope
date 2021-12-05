from collections import defaultdict, deque
from ms_deisotope.data_source.scan.mobility_frame import IonMobilityFrame
from ms_deisotope.data_source.scan import ScanBunch, ScanBase

from .scan.scan_iterator import ITERATION_MODE_GROUPED, ITERATION_MODE_SINGLE
from .scan.loader import ScanIterator
from .scan.base import ScanBunch
from .metadata.scan_traits import FAIMS_compensation_voltage


class ScanIteratorProxyBase(object):
    """A wrapper around a :class:`~.RandomAccessScanSource` that alters its behavior
    in some way, either limiting the range it iterates over, or filtering the items
    it would yield.

    This is an abstract base class. It should be extended to be used.

    This class re-uses the same file handle that as :attr:`scan_source`, and should
    not be used interleaved with the original source.

    Attributes
    ----------
    scan_source: :class:`~.RandomAccessScanSource`
        The object to use to load :class:`~.Scan` and :class:`~.ScanBunch` instances.
    iteration_mode: str
        A string denoting :const:`~.ITERATION_MODE_GROUPED` or :const:`~.ITERATION_MODE_SINGLE`
        that controls whether :class:`~.ScanBunch` or :class:`~.Scan` are produced
        by iteration.
    """

    def __init__(self, scan_source, *args, **kwargs):
        self.scan_source = scan_source

    @property
    def iteration_mode(self):
        '''A string denoting :const:`~.ITERATION_MODE_GROUPED` or :const:`~.ITERATION_MODE_SINGLE`
        that controls whether :class:`~.ScanBunch` or :class:`~.Scan` are produced
        by iteration.

        Returns
        -------
        str
        '''
        return self.scan_source.iteration_mode

    def has_ms1_scans(self):
        '''Checks if this :class:`ScanDataSource` contains MS1 spectra.

        Returns
        -------
        :class:`bool` or :const:`None`
            Returns a boolean value if the presence of MS1 scans is known for certain, or :const:`None`
            if it cannot be determined in the case of missing metadata.
        '''
        return self.scan_source.has_ms1_scans()

    def has_msn_scans(self):
        '''Checks if this :class:`ScanDataSource` contains MSn spectra.

        Returns
        -------
        :class:`bool` or :const:`None`
            Returns a boolean value if the presence of MSn scans is known for certain, or :const:`None`
            if it cannot be determined in the case of missing metadata.
        '''
        return self.scan_source.has_msn_scans()

    def next(self):
        '''Advance the iterator, fetching the next :class:`~.ScanBunch` or :class:`~.ScanBase`
        depending upon iteration strategy.

        Returns
        -------
        :class:`~.ScanBunch` or :class:`~.ScanBase`
        '''
        raise NotImplementedError()

    def __next__(self):
        '''Advance the iterator, fetching the next :class:`~.ScanBunch` or :class:`~.ScanBase`
        depending upon iteration strategy.

        Returns
        -------
        :class:`~.ScanBunch` or :class:`~.ScanBase`
        '''
        return self.next()

    def __iter__(self):
        return self


ScanIterator.register(ScanIteratorProxyBase)


class QueryIteratorBase(ScanIteratorProxyBase):
    """A base class for types which iterate over a subset of a :class:`~.RandomAccessScanSource`.

    Attributes
    ----------
    grouped: bool
        Whether or not to produce grouped or single scans.

    """
    def __init__(self, scan_source, grouped=True, *args, **kwargs):
        super(QueryIteratorBase, self).__init__(scan_source, *args, **kwargs)
        self.grouped = grouped

    @property
    def iteration_mode(self):
        if self.grouped:
            return ITERATION_MODE_GROUPED
        else:
            return ITERATION_MODE_SINGLE

    def make_iterator(self, iterator=None, grouped=False):
        """Configure the :class:`ScanIterator`'s behavior, selecting it's iteration strategy over
        either its default iterator or the provided ``iterator`` argument.

        Parameters
        ----------
        iterator : Iterator, optional
            Unused. Included for compatibility with :class:`ScanIterator` API
        grouped : bool, optional
            Whether the iterator should be grouped and produce :class:`.ScanBunch` objects
            or single :class:`.Scan`. If :const:`None` is passed, :meth:`has_ms1_scans` will be
            be used instead. Defaults to :const:`None`.
        """
        self.grouped = grouped
        self._make_iterator(grouped=grouped)

    def _make_iterator(self, grouped=False):
        self.scan_source.reset()
        self.scan_source.make_iterator()

    def reset(self):
        '''Reset the iterator, if possible, and clear any caches.

        Resets the underlying :class:`~.RandomAccessScanSource`
        '''
        self.scan_source.reset()
        self.make_iterator(grouped=self.grouped)


class TimeIntervalIterator(QueryIteratorBase):
    """Query over a retention time interval.

    Attributes
    ----------
    start: float
        The time to start the query iterator from
    end: float
        The time to end the query iterator at

    """
    def __init__(self, scan_source, start=None, end=None, grouped=True, *args, **kwargs):
        super(TimeIntervalIterator, self).__init__(
            scan_source, grouped, * args, **kwargs)
        self.start = start
        self.end = end
        self._update_bounds()
        self.make_iterator(grouped=grouped)

    def _update_bounds(self):
        if self.start is None:
            try:
                first_entry = self.scan_source[0]
                self.start = first_entry.scan_time
            except (AttributeError, TypeError):
                self.start = 0
        if self.end is None:
            try:
                last_entry = self.scan_source[-1]
                self.end = last_entry.scan_time
            except (AttributeError, TypeError):
                self.end = float('inf')

    def _make_iterator(self, grouped=False):
        self.scan_source.start_from_scan(rt=self.start, grouped=self.grouped)

    def next(self):
        result = next(self.scan_source)
        if self.grouped:
            if result.precursor.scan_time > self.end:
                raise StopIteration()
            return result
        else:
            if result.scan_time > self.end:
                raise StopIteration()
            return result


class IndexIntervalIterator(QueryIteratorBase):
    """Query over a scan index interval.

    Attributes
    ----------
    start: int
        The index to start the query iterator from
    end: int
        The index to end the query iterator at

    """

    def __init__(self, scan_source, start=None, end=None, grouped=True, *args, **kwargs):
        super(IndexIntervalIterator, self).__init__(
            scan_source, grouped, * args, **kwargs)
        self.start = start
        self.end = end
        self._update_bounds()
        self.make_iterator(grouped=grouped)

    def _update_bounds(self):
        if self.start is None:
            try:
                first_entry = self.scan_source[0]
                self.start = first_entry.index
            except (AttributeError, TypeError):
                self.start = 0
        if self.end is None:
            try:
                last_entry = self.scan_source[-1]
                self.end = last_entry.index
            except (AttributeError, TypeError):
                self.end = float('inf')

    def _make_iterator(self, grouped=False):
        self.scan_source.start_from_scan(index=self.start, grouped=self.grouped)

    def next(self):
        result = next(self.scan_source)
        if self.grouped:
            if result.precursor.index > self.end:
                raise StopIteration()
            return result
        else:
            if result.index > self.end:
                raise StopIteration()
            return result


def scan_range(scan_source, time=None, index=None, grouped=True, *args, **kwargs):
    """Create an iterator proxy over `scan_source` spanning a specified range in time
    or index.

    If neither `time` nor `index` is provided, `scan_source` is returned unchanged.
    If both are provided, an error is thrown.

    Parameters
    ----------
    scan_source : :class:`~.RandomAccessScanSource`
        The scan source to iterate over
    time : tuple[float, float], optional
        The start and stop times to use
    index : tuple[int, int], optional
        The start and top scan indices to use
    grouped : bool, optional
        Whether or not to create a :class:`~.ScanBunch` iterator or a :class:`~.Scan` iterator
    *args:
    **kwargs:
        Forwarded to the iterator proxy.

    Returns
    -------
    :class:`ScanIterator`

    Raises
    ------
    ValueError:
        If both `time` and `index` are provided.
    """
    if time and index:
        raise ValueError("Only one of time and index intervals may be specified")
    if time:
        return TimeIntervalIterator(
            scan_source, time[0], time[1], grouped=grouped, *args, **kwargs)
    elif index:
        return IndexIntervalIterator(
            scan_source, index[0], index[1], grouped=grouped, *args, **kwargs)
    else:
        return scan_source


class ScanIteratorFilterBase(ScanIteratorProxyBase):
    """A base class for types which filter or transform scans produced
    by a :class:`~.ScanIterator` or :class:`~.RandomAccessScanSource`.

    """
    def __init__(self, scan_source, *args, **kwargs):
        super(ScanIteratorFilterBase, self).__init__(scan_source, *args, **kwargs)
        self._generator = None

    def filter_scan_bunch(self, scan_group, **kwargs):
        """Filter a scan bunch

        Parameters
        ----------
        scan : :class:`~.ScanBunch`
            The scan bunch to filter

        Yields
        ------
        :class:`~.ScanBunch`
        """
        raise NotImplementedError()

    def filter_scan(self, scan, **kwargs):
        """Filter a single scan

        Parameters
        ----------
        scan : :class:`~.Scan`
            The single scan to filter

        Yields
        ------
        :class:`~.Scan`
        """
        raise NotImplementedError()

    def filter(self, result, **kwargs):
        """Filter scans, producing only those that passed the criteria.

        This returns a generator which yields the same type as `result`

        This method will dispatch to :meth:`filter_scan` or :meth:`filter_scan_bunch`,
        whichever is appropriate for `result`

        Parameters
        ----------
        result : :class:`~.Scan` or :class:`~.ScanBunch`
            The scan or bunch to filter.

        Yields
        ------
        :class:`~.Scan` or :class:`~.ScanBunch`
        """
        if self.scan_source.iteration_mode == ITERATION_MODE_GROUPED:
            return self.filter_scan_bunch(result, **kwargs)
        else:
            return self.filter_scan(result, **kwargs)

    def _generate(self):
        for batch in self.scan_source:
            for result in self.filter(batch):
                yield result

    def next(self):
        if self._generator is None:
            self._generator = self._generate()
        return next(self._generator)


class _PredicateFilterIterator(ScanIteratorFilterBase):
    def _check(self, scan):
        raise NotImplementedError()

    def filter_scan_bunch(self, scan_group, **kwargs):
        if scan_group.products:
            yield scan_group.__class__(
                scan_group.precursor if self._check(
                    scan_group.precursor) else None,
                [scan for scan in scan_group.products if self._check(scan)])

    def filter_scan(self, scan, **kwargs):
        if self._check(scan):
            yield scan

    @property
    def key(self):
        raise NotImplementedError()


class PolarityFilter(_PredicateFilterIterator):
    """A Scan Iterator that filters out scans with a polarity
    which do not match :attr:`polarity`

    Attributes
    ----------
    polarity: int
        The polarity to match

    """

    def __init__(self, scan_source, polarity, *args, **kwargs):
        super(PolarityFilter, self).__init__(scan_source, *args, **kwargs)
        self.polarity = polarity

    def _check(self, scan):
        return scan.polarity == self.polarity

    @property
    def key(self):
        return self.polarity


class MSLevelFilter(_PredicateFilterIterator):
    """A Scan Iterator that filters out scans with MS levels
    which do not match :attr:`ms_level`

    Attributes
    ----------
    ms_level: int
        The MS level to require matching

    """
    def __init__(self, scan_source, ms_level, *args, **kwargs):
        super(MSLevelFilter, self).__init__(scan_source, *args, **kwargs)
        self.ms_level = ms_level

    def _check(self, scan):
        return scan.ms_level == self.ms_level

    @property
    def key(self):
        return self.ms_level


class MassAnalyzerFilter(_PredicateFilterIterator):
    """A Scan Iterator that filters out scans which were not performed using
    :attr:`mass_analyzer`

    Attributes
    ----------
    mass_analyzer: :class:`~.Component`
        The mass analyzer to match
    """
    def __init__(self, scan_source, mass_analyzer, *args, **kwargs):
        super(MassAnalyzerFilter, self).__init__(scan_source, *args, **kwargs)
        from .metadata.instrument_components import analyzer_types, Component
        if not isinstance(mass_analyzer, Component):
            mass_analyzer = analyzer_types[mass_analyzer]
        self.mass_analyzer = mass_analyzer

    def _check(self, scan):
        ic = scan.instrument_configuration
        if ic is None:
            return False
        for analyzer in ic.analyzers:
            if analyzer.is_a(self.mass_analyzer):
                return True
        return False

    @property
    def key(self):
        return self.mass_analyzer


class CallableFilter(_PredicateFilterIterator):
    """A Scan Iterator that filters out scans which do not pass
    a user-provided callable

    Attributes
    ----------
    filter_fn: :class:`Callable`
        The callable object to use to test if a scan should be used
    """

    def __init__(self, scan_source, filter_fn, *args, **kwargs):
        super(CallableFilter, self).__init__(scan_source, *args, **kwargs)
        self.filter_fn = filter_fn

    def _check(self, scan):
        return self.filter_fn(scan)

    @property
    def key(self):
        return self.filter_fn


filter_scans = CallableFilter


class MS1MergingTransformer(ScanIteratorProxyBase):

    def __init__(self, scan_source, ms1_bin_width=10, mz_resolution=None, *args, **kwargs):
        super(MS1MergingTransformer, self).__init__(scan_source, *args, **kwargs)
        self._generator = None
        self.ms1_bin_width = ms1_bin_width
        self.mz_resolution = mz_resolution
        self.ms1_buffer = []
        self.msn_buffer = []

    def _redirect_precursor_ids(self, products, precursor):
        level_to_redirect = precursor.ms_level + 1
        # Don't re-direct MS3+ scans
        for product in products:
            if product.ms_level == level_to_redirect:
                product.precursor_information.precursor_scan_id = precursor.id
        return products

    def _generate(self):
        half_width = self.ms1_bin_width // 2
        for batch in self.scan_source:
            if len(self.ms1_buffer) == self.ms1_bin_width:
                ref_scan = self.ms1_buffer[0]
                merged_scan = ref_scan.average_with(self.ms1_buffer[1:], self.mz_resolution)
                # Scale up to be the sum not the mean
                merged_scan.arrays *= self.ms1_bin_width
                bunch = ScanBunch(merged_scan, self._redirect_precursor_ids(self.msn_buffer, merged_scan))
                yield bunch
                self.msn_buffer = []
                self.ms1_buffer = self.ms1_buffer[half_width:]

            self.ms1_buffer.append(batch.precursor)
            self.msn_buffer.extend(batch.products)

        if self.ms1_buffer:
            ref_scan = self.ms1_buffer[0]
            merged_scan = ref_scan.average_with(self.ms1_buffer[1:], self.mz_resolution)
            # Scale up to be the sum not the mean
            merged_scan.arrays *= self.ms1_bin_width
            bunch = ScanBunch(merged_scan, self._redirect_precursor_ids(self.msn_buffer, merged_scan))
            yield bunch


    def next(self):
        if self._generator is None:
            self._generator = self._generate()
        return next(self._generator)


class TimeOrderMergingIterator(object):
    def __init__(self, sources):
        self.sources = list(sources)
        self._init_heads()

    def _init_heads(self):
        n = len(self.sources)
        self.heads = [None] * n
        for i in range(n):
            self._advance(i)

    def __next__(self):
        result = self.get_next_head()
        if result is None:
            raise StopIteration()
        return result

    next = __next__

    def __iter__(self):
        return self

    def count_exhausted(self):
        n = 0
        for h in self.heads:
            n += h is None
        return n

    def _advance(self, i):
        try:
            self.heads[i] = next(self.sources[i])
        except StopIteration:
            self.heads[i] = None

    def get_next_head(self):
        best = None
        best_i = None
        best_time = float('inf')
        for i, head in enumerate(self.heads):
            if head is None:
                continue
            time = self.get_time(head)
            assert time is not None
            if time < best_time:
                best = head
                best_i = i
                best_time = time
        if best is None:
            return best
        i = best_i
        self._advance(i)
        return best

    def get_time(self, obj):
        if isinstance(obj, ScanBunch):
            if obj.precursor is not None:
                return self.get_time(obj.precursor)
            elif obj.products:
                return self.get_time(obj.products[0])
        elif isinstance(obj, ScanBase):
            return obj.scan_time
        elif isinstance(obj, IonMobilityFrame):
            return obj.time
        else:
            # Just guessing at this point
            return obj.scan_time


class FAIMSFilter(_PredicateFilterIterator):
    def __init__(self, scan_source, compensation_voltage, *args, **kwargs):
        super(FAIMSFilter, self).__init__(scan_source, *args, **kwargs)
        self.compensation_voltage = compensation_voltage

    def _check(self, scan):
        if scan.has_ion_mobility():
            if scan.ion_mobility_type == FAIMS_compensation_voltage:
                if scan.drift_time == self.compensation_voltage:
                    return scan

    @property
    def key(self):
        return self.compensation_voltage


class QueueIterator(ScanIteratorProxyBase):
    def __init__(self, scan_source, data=None):
        self.scan_source = scan_source
        self.data = deque(data or [])

    def __next__(self):
        return self.data.popleft()

    def has_value(self):
        return bool(self.data)

    def _feed(self, value):
        self.data.append(value)

    def __len__(self):
        return len(self.data)


class DemultiplexingIteratorBase(ScanIteratorProxyBase):
    def __init__(self, scan_source, buffer_size=10, *args, **kwargs):
        super(DemultiplexingIteratorBase, self).__init__(
            scan_source, *args, **kwargs)
        self.channels = dict()
        self.buffer_size = buffer_size
        self.has_more = self.feed()

    @property
    def iteration_mode(self):
        return ITERATION_MODE_SINGLE

    def _new_channel(self, scan):
        raise NotImplementedError()

    def __getitem__(self, key):
        return self.channels[key]

    def __iter__(self):
        return self

    def _get_channel_name(self, scan_filter):
        return scan_filter.key

    def feed(self):
        loading = sum([len(v[1]) for v in self.channels.values()])
        while loading < self.buffer_size:
            try:
                batch = next(self.scan_source)
            except StopIteration:
                return False
            if self.scan_source.iteration_mode == ITERATION_MODE_SINGLE:
                for channel_name, (scan_filter, queue) in self.channels.items():
                    if scan_filter._check(batch):
                        queue._feed(batch)
                        break
                else:
                    self._new_channel(batch)
            else:
                for scan in [batch.precursor] + batch.products:
                    for channel_name, (scan_filter, queue) in self.channels.items():
                        if scan_filter._check(scan):
                            queue._feed(scan)
                            break
                    else:
                        self._new_channel(scan)
            loading = sum([len(v[1]) for v in self.channels.values()])
        return True

    def __next__(self):
        if not self.has_more:
            raise StopIteration()
        meta = {}
        for channel_name, (scan_filter, queue) in self.channels.items():
            if queue.has_value():
                meta[channel_name] = next(scan_filter)
            else:
                meta[channel_name] = None
        self.has_more = self.feed()
        return meta


class FAIMSDemultiplexingIterator(DemultiplexingIteratorBase):
    def __init__(self, scan_source, buffer_size=10, *args, **kwargs):
        super(FAIMSDemultiplexingIterator, self).__init__(scan_source, buffer_size, *args, **kwargs)

    def _new_channel(self, scan):
        queue = QueueIterator(self)
        scan_filter = FAIMSFilter(queue, scan.drift_time)
        queue._feed(scan)
        self.channels[self._get_channel_name(scan_filter)] = (scan_filter, queue)

