from .scan.scan_iterator import ITERATION_MODE_GROUPED, ITERATION_MODE_SINGLE


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

    def _make_iterator(self, grouped=False):
        raise NotImplementedError()

    def reset(self):
        self.scan_source.reset()
        self._make_iterator(self.grouped)


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
        self._make_iterator(self, grouped)

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


class MSLevelFilter(ScanIteratorFilterBase):
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

    def filter_scan_bunch(self, scan_group, **kwargs):
        if scan_group.products:
            yield scan_group.__class__(
                scan_group.precursor if self._check(scan_group.precursor) else None,
                [scan for scan in scan_group.products if self._check(scan)])

    def filter_scan(self, scan, **kwargs):
        if self._check(scan):
            yield scan
