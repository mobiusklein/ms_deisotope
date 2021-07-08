from collections import OrderedDict
import bisect
from .common import (
    ScanBunch, ScanDataSource, RandomAccessScanSource,
    Scan, ProcessedScan, WrappedScan)

try:
    range = xrange
except NameError:
    pass


class MemoryScanInterface(ScanDataSource):
    '''A dummy :class:`~.ScanDataSource` that stores as much data as
    possible on the :class:`~.Scan` instance, with the exception of
    :attr:`~.Scan.index`.
    '''

    @classmethod
    def make_scan(cls, arrays=None, ms_level=None, id=None, index=None, scan_time=None,
                  is_profile=False, polarity=1, precursor_information=None, activation=None,
                  isolation_window=None, annotations=None, acquisition_information=None,
                  instrument_configuration=None, peak_set=None, deconvoluted_peak_set=None,
                  **kwargs):
        '''
        Parameters
        ----------
        arrays: :class:`RawDataArrays`
            A pair of :class:`numpy.ndarray` objects corresponding to the raw m/z and intensity data points
        ms_level: int
            The degree of fragmentation performed. 1 corresponds to a MS1 or "Survey" scan, 2 corresponds
            to MS/MS, and so on. If :attr:`ms_level` > 1, the scan is considered a "tandem scan" or "MS^n" scan
        id: str
            The unique identifier for this scan as given by the source
        index: int
            The integer number indicating how many scans were acquired prior to this scan.
        scan_time: float
            The time the scan was acquired during data acquisition. The unit of time will always be minutes.
        is_profile: bool
            Whether this scan's raw data points corresponds to a profile scan or whether the raw data was
            pre-centroided.
        polarity: int
            If the scan was acquired in positive mode, the value ``+1``.  If the scan was acquired in negative
            mode, the value ``-1``. May be used to indicating how to calibrate charge state determination methods.
        precursor_information: :class:`PrecursorInformation` or None
            Descriptive metadata for the ion which was chosen for fragmentation, and a reference to
            the precursor scan
        activation: :class:`.ActivationInformation` or None
            If this scan is an MS^n scan, this attribute will contain information about the process
            used to produce it from its parent ion.
        instrument_configuration: :class:`~.InstrumentInformation`
            The instrument configuration used to acquire this scan.
        acquisition_information: :class:`.ScanAcquisitionInformation` or None
            Describes the type of event that produced this scan, as well as the scanning method
            used.
        isolation_window: :class:`.IsolationWindow` or None
            Describes the range of m/z that were isolated from a parent scan to create this scan
        annotations: dict
            A set of key-value pairs describing the scan not part of the standard interface
        peak_set : :class:`ms_peak_picker.PeakSet` or None
            Picked peaks and (possibly) associated raw data points as produced by :meth:`pick_peaks`.
            Will be `None` if peak picking has not been done.
        deconvoluted_peak_set : :class:`ms_deisotope.DeconvolutedPeakSet` or None
            Deconvoluted peaks resulting from charge state deconvolution and deisotoping. Will
            be `None` if deconvolution has not been done.
        '''
        if annotations is None:
            annotations = {}
        if arrays is None:
            arrays = [[], []]
        annotations.update(kwargs)
        scan = WrappedScan({}, cls(), arrays, None, annotations=annotations)
        scan.arrays = arrays
        scan.ms_level = ms_level
        scan.id = str(id)
        scan.title = str(id)
        scan.index = index or -1
        scan.scan_time = scan_time or -1
        scan.is_profile = bool(is_profile)
        scan.polarity = int(polarity)
        scan.activation = activation
        scan.isolation_window = isolation_window
        scan.acquisition_information = acquisition_information
        scan.instrument_configuration = instrument_configuration
        scan.peak_set = peak_set
        scan.deconvoluted_peak_set = deconvoluted_peak_set
        scan.precursor_information = precursor_information
        return scan


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
        return None

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
        return None

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
        return None

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
        try:
            return self._scan_index_map[scan].id
        except AttributeError:
            return None

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
        return None

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
        try:
            return self._scan_index_map[scan].ms_level
        except AttributeError:
            return None

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
        return None

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
        return None

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
        return None

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
        return None


make_scan = MemoryScanInterface.make_scan


class ScanCollection(MemoryScanInterface, RandomAccessScanSource):
    '''A :class:`~.RandomAccessScanSource` implementation which contains scan objects
    materialized from other sources or which have already been fully specified in memory.
    '''

    def __init__(self, scans, binds=True, **kwargs):
        self._scans = sorted(scans, key=lambda x: x.scan_time)
        self._binds = binds
        self._producer = None
        self._scan_id_map = dict()
        self._scan_index_map = OrderedDict()
        self._index = OrderedDict()
        self._build_indices()
        self.make_iterator()

    @property
    def _scan_cache(self):
        return self._scan_id_map

    @classmethod
    def build(cls, scans, binds=True, **kwargs):
        '''Construct an :class:`ScanCollection`
        '''
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
                scans = tuple(s.clone()._load() for s in scans)
            elif isinstance(scan, ProcessedScan):
                scans = tuple(s.clone() for s in scans)
            else:
                raise TypeError("Cannot build an in-memory scan source from {}".format(
                    type(scan)))
        except IndexError:
            raise ValueError("Must pass a non-empty iterable")
        # This alters the index of the scan while it may still be bound
        # to another scan source for attribute lookup. At the time of
        # writing this is needed to properly check bounds on ScanIterator
        # methods, but may lead to problems if the index is needed to resolve
        # other attributes.
        for i, scan in enumerate(scans):
            scan.index = i
        return cls(scans, binds=binds, **kwargs)

    def reset(self):
        self.make_iterator(None)

    def _build_indices(self):
        self._scan_id_map = dict()
        self._scan_index_map = OrderedDict()
        self._index = OrderedDict()

        for scan in self._scans:
            if self._binds:
                scan.bind(self)
            self._scan_id_map[scan.id] = scan
            self._scan_index_map[scan.index] = scan
            self._index[scan.id] = scan.index

    @property
    def index(self):
        return self._index

    def __len__(self):
        return len(self.index)

    def _validate(self, scan):
        return True

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

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        if scan_id is not None:
            index = self.get_scan_by_id(scan_id).index
        elif rt is not None:
            index = self.get_scan_by_time(rt).index
        scan = self.get_scan_by_index(index)
        if require_ms1:
            while scan.ms_level > 1:
                index -= 1
                scan = self.get_scan_by_index(index)

        self.make_iterator(self._make_scan_index_producer(index), grouped=grouped)
        return self


make_scan_collection = ScanCollection.build
