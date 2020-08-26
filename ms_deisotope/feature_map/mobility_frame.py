from abc import ABCMeta, abstractmethod


class IonMobilitySource(ABCMeta):
    @abstractmethod
    def _frame_id(self, data):
        raise NotImplementedError()

    @abstractmethod
    def _frame_index(self, data):
        raise NotImplementedError()

    def _frame_time(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.scan_time

    def _frame_ms_level(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.ms_level

    @abstractmethod
    def _frame_start_scan_index(self, data):
        raise NotImplementedError()

    @abstractmethod
    def _frame_end_scan_index(self, data):
        raise NotImplementedError()

    def _frame_precursor_information(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.precursor_information

    def _frame_activation(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.activation

    def _frame_isolation_window(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.isolation_window

    def _frame_polarity(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.polarity


class IonMobilitySourceRandomAccessFrameSource(IonMobilitySource):
    @abstractmethod
    def get_frame_by_index(self, index):
        raise NotImplementedError()

    @abstractmethod
    def get_frame_by_time(self, time):
        raise NotImplementedError()

    @abstractmethod
    def _validate_frame(self, data):
        raise NotImplementedError()

    @abstractmethod
    def _make_frame(self, data):
        raise NotImplementedError()

    @abstractmethod
    def _cache_frame(self, frame):
        raise NotImplementedError()

    @abstractmethod
    def _default_frame_iterator(self, start_index=None):
        raise NotImplementedError()

    @abstractmethod
    def make_frame_iterator(self, iterator=None, grouped=False):
        raise NotImplementedError()

    @abstractmethod
    def start_from_frame(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        raise NotImplementedError()


class IonMobilityFrame(object):
    def __init__(self, data, source):
        self._data = data
        self.source = source
        self.features = None
        self.deconvoluted_features = None

    def __repr__(self):
        template = ("{self.__class__.__name__}({self.id}, index={self.index},"
                    " time={self.time}, ms_level={self.ms_level})")
        return template.format(self=self)

    @property
    def id(self):
        return self.source._frame_id(self._data)

    @property
    def index(self):
        return self.source._frame_index(self._data)

    @property
    def time(self):
        return self.source._frame_time(self._data)

    @property
    def ms_level(self):
        return self.source._frame_ms_level(self._data)

    @property
    def polarity(self):
        return self.source._frame_polarity(self._data)

    @property
    def start_scan_index(self):
        return self.source._frame_start_scan_index(self._data)

    @property
    def end_scan_index(self):
        return self.source._frame_end_scan_index(self._data)

    @property
    def precursor_information(self):
        return self.source._frame_precursor_information(self._data)

    @property
    def activation(self):
        return self.source._frame_activation(self._data)

    @property
    def isolation_window(self):
        return self.source._frame_isolation_window(self._data)

    def scans(self):
        scans = []
        for i in range(self.start_scan_index, self.end_scan_index):
            scan = self.source.get_scan_by_index(i)
            scans.append(scan)
        return scans

    def extract_features(self, error_tolerance=1.5e-5, max_gap_size=0.25, min_size=2, **kwargs):
        from ms_deisotope.feature_map import feature_map
        scans = self.scans()
        lff = feature_map.LCMSFeatureForest(error_tolerance=error_tolerance)
        for scan in scans:
            scan.pick_peaks(**kwargs)
            for peak in scan:
                lff.handle_peak(peak, scan.drift_time)
        lff.split_sparse(max_gap_size, min_size).smooth_overlaps(
            error_tolerance)
        self.features = lff
        return self

    def deconvolute_features(self, averagine=None, scorer=None, truncate_after=0.95,
                             minimum_intensity=5, minimum_size=2, maximum_gap_size=0.25, **kwargs):
        from ms_deisotope.feature_map import feature_processor
        if averagine is None:
            from ms_deisotope.averagine import peptide
            averagine = peptide
        if scorer is None:
            from ms_deisotope.scoring import MSDeconVFitter
            scorer = MSDeconVFitter(1)
        if self.features is None:
            raise ValueError(
                "IM-MS Features must be extracted before they can be charge state deconvoluted")
        decon = feature_processor.LCMSFeatureProcessor(
            self.features, averagine, scorer, minimum_size=minimum_size,
            maximum_time_gap=maximum_gap_size)
        self.deconvoluted_features = decon.deconvolute(
            minimum_intensity=minimum_intensity, truncate_after=truncate_after, **kwargs)
        return self
