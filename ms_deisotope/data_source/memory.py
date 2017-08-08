from collections import OrderedDict
from .common import (
    ScanBunch, ScanIterator, ScanDataSource, RandomAccessScanSource,
    Scan)

try:
    range = xrange
except NameError:
    pass


class MemoryScanInterface(ScanDataSource):
    def _scan_arrays(self, scan):
        """Returns raw data arrays for m/z and intensity

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        mz: np.array
            An array of m/z values for this scan
        intensity: np.array
            An array of intensity values for this scan
        """
        raise NotImplementedError()

    def _precursor_information(self, scan):
        """Returns information about the precursor ion,
        if any, that this scan was derived form.

        Returns `None` if this scan has no precursor ion

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        PrecursorInformation
        """
        raise NotImplementedError()

    def _scan_title(self, scan):
        """Returns a verbose name for this scan, if one
        were stored in the file. Usually includes both the
        scan's id string, as well as information about the
        original file and format.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        str
        """
        raise NotImplementedError()

    def _scan_id(self, scan):
        """Returns the scan's id string, a unique
        identifier for this scan in the context of
        the data file it is recordered in

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        str
        """
        raise NotImplementedError()

    def _scan_index(self, scan):
        """Returns the base 0 offset from the start
        of the data file in number of scans to reach
        this scan.

        If the original format does not natively include
        an index value, this value may be computed from
        the byte offset index.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        int
        """
        raise NotImplementedError()

    def _ms_level(self, scan):
        """Returns the degree of exponential fragmentation
        used to produce this scan.

        1 refers to a survey scan of unfragmented ions, 2
        refers to a tandem scan derived from an ms level 1
        ion, and so on.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        int
        """
        raise NotImplementedError()

    def _scan_time(self, scan):
        """Returns the time in minutes from the start of data
        acquisition to when this scan was acquired.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        float
        """
        raise NotImplementedError()

    def _is_profile(self, scan):
        """Returns whether the scan contains profile data (`True`)
        or centroided data (`False`).
        
        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`
        
        Returns
        -------
        bool
        """
        raise NotImplementedError()

    def _polarity(self, scan):
        """Returns whether this scan was acquired in positive mode (+1)
        or negative mode (-1).
        
        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`
        
        Returns
        -------
        int
        """
        raise NotImplementedError()

    def _activation(self, scan):
        """Returns information about the activation method used to
        produce this scan, if any.

        Returns `None` for MS1 scans
        
        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`
        
        Returns
        -------
        ActivationInformation
        """
        raise NotImplementedError()


class MemoryScanLoader(MemoryScanInterface, ScanIterator, RandomAccessScanSource):

    def __init__(self, scans, **kwargs):
        self._scans = tuple(sorted(scans, key=lambda x: x.scan_time))
        self._producer = None
        self._scan_id_map = dict()
        self._scan_index_map = OrderedDict()
        self._index = OrderedDict()
        self._build_indices()
        self.make_iterator()

    @classmethod
    def build(cls, scans):
        scans = tuple(scans)
        try:
            scan = scans[0]
            if isinstance(scan, ScanBunch):
                def generate():
                    for bunch in scans:
                        yield bunch.precursor
                        for product in bunch.products:
                            yield product
                return cls.build(generate())
            elif isinstance(scan, Scan):
                scans = tuple(s.pack() for s in scans)
            else:
                scans = tuple(s.clone() for s in scans)
        except IndexError:
            raise ValueError("Must pass a non-empty iterable")
        for i, scan in enumerate(scans):
            scan.index = i
        return cls(scans)

    def reset(self):
        self.make_iterator(None)

    def _build_indices(self):
        self._scan_id_map = dict()
        self._scan_index_map = OrderedDict()
        self._index = OrderedDict()

        for scan in self._scans:
            self._scan_id_map[scan.id] = scan
            self._scan_index_map[scan.index] = scan
            self._index[scan.id] = scan.index

    @property
    def index(self):
        return self._index

    def _make_scan_index_producer(self, start_index=None, start_time=None):
        if start_index is not None:
            return range(start_index, len(self._scans))
        elif start_time is not None:
            start_index = self.get_scan_by_time(start_time).index
            while start_index != 0:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            return range(start_index, len(self._scans))
        else:
            return range(0, len(self._scans))

    def _make_default_iterator(self):
        return self._make_scan_index_producer()

    def next(self):
        return next(self._producer)

    def get_scan_by_id(self, scan_id):
        return self._scan_id_map[scan_id]

    def get_scan_by_time(self, time):
        scan_ids = tuple(self.index)
        lo = 0
        hi = len(scan_ids)
        while hi != lo:
            mid = (hi + lo) // 2
            sid = scan_ids[mid]
            scan = self.get_scan_by_id(sid)
            scan_time = scan.scan_time
            if scan_time == time:
                return scan
            elif (hi - lo) == 1:
                return scan
            elif scan_time > time:
                hi = mid
            else:
                lo = mid

    def get_scan_by_index(self, index):
        return self._scan_index_map[index]

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True):
        if scan_id is not None:
            index = self.get_scan_by_id(scan_id).index
        elif rt is not None:
            index = self.get_scan_by_time(rt).index
        scan = self.get_scan_by_index(index)
        if require_ms1:
            while scan.ms_level > 1:
                index -= 1
                scan = self.get_scan_by_index(index)

        self.make_iterator(self._make_scan_index_producer(index))
        return self
