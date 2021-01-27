from abc import abstractmethod
from weakref import WeakValueDictionary
from itertools import chain

from ms_deisotope.peak_set import IonMobilityDeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.peak_dependency_network import IntervalTreeNode, Interval

class IonMobilitySource(object):
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

    def _frame_drift_times(self, data):
        scans = []
        for i in range(self._frame_start_scan_index(data), self._frame_end_scan_index(data)):
            scan = self.get_scan_by_index(i)
            scans.append(scan.drift_time)
        return scans

    def _frame_acquisition_information(self, data):
        scan = self.get_scan_by_index(self._frame_start_scan_index(data))
        acq = scan.acquisition_information
        if acq is None:
            return acq
        acq = acq.copy()
        for event in acq:
            im = event._ion_mobility
            mt = im.ion_mobility_type()
            if mt is not None:
                im.remove_ion_mobility_type(mt)
        return acq


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

    def initialize_frame_cache(self):
        '''Initialize a cache which keeps track of which :class:`~.IonMobilityFrame`
        objects are still in memory using a :class:`weakref.WeakValueDictionary`.

        When a frame is requested, if the frame object is found in the cache, the
        existing object is returned rather than re-read from disk.
        '''
        self._frame_cache = WeakValueDictionary()

    @property
    def frame_cache(self):
        '''A :class:`weakref.WeakValueDictionary` mapping used to retrieve
        frames from memory if available before re-reading them from disk.
        '''
        return self._frame_cache

    @frame_cache.setter
    def frame_cache(self, value):
        self._frame_cache = value


class IonMobilityFrame(object):
    '''A :class:`IonMobilityFrame` represents an single time point acquisition of
    multiple mass spectra across multiple ion mobility drift time points. Because
    of its drift time by m/z structure, it does not have 1-D peak sets, but pseudo-
    :class:`~.LCMSFeatureMap`-like data structures which conserve the over-time
    property.
    '''
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

    @property
    def acquisition_information(self):
        return self.source._frame_acquisition_information(self._data)

    def scans(self):
        scans = []
        for i in range(self.start_scan_index, self.end_scan_index):
            scan = self.source.get_scan_by_index(i)
            scans.append(scan)
        return scans

    @property
    def drift_times(self):
        return self.source._frame_drift_times(self._data)

    def get_scan_by_drift_time(self, drift_time):
        dt_axis = self.drift_times
        lo = 0
        hi = n = len(dt_axis)
        best_match = None
        best_error = float('inf')

        # Handles infinities
        if drift_time > dt_axis[-1]:
            return self.source.get_scan_by_index(self.end_scan_index - 1)

        while hi != lo:
            mid = (hi + lo) // 2
            dt = dt_axis[mid]
            err = abs(dt - drift_time)
            if err < best_error:
                best_error = err
                best_match = mid
            if (dt == drift_time) or (hi - 1 == lo):
                i = mid - 1
                while i >= 0:
                    err = abs(dt_axis[i] - drift_time)
                    if err < best_error:
                        best_error = err
                        best_match = i
                        i -= 1
                    else:
                        break
                i = mid + 1
                while i < n:
                    err = abs(dt_axis[i] - drift_time)
                    if err < best_error:
                        best_error = err
                        best_match = i
                        i += 1
                    else:
                        break
                return self.source.get_scan_by_index(self.start_scan_index + best_match)
            elif dt > drift_time:
                hi = mid
            else:
                lo = mid

    def extract_features(self, error_tolerance=1.5e-5, max_gap_size=0.25, min_size=2, average=0, dx=0.001, **kwargs):
        from ms_deisotope.feature_map import feature_map
        scans = self.scans()
        lff = feature_map.LCMSFeatureForest(error_tolerance=error_tolerance)
        n = len(scans)
        for i, scan in enumerate(scans):
            if average:
                acc = []
                if i > 0:
                    for j in range(max((i - average, 0)), i):
                        acc.append(scans[j])
                    for j in range(i + 1, min(i + average + 1, n)):
                        acc.append(scans[j])
                scan = scan.average_with(acc, dx=dx)
            scan.pick_peaks(**kwargs)
            for peak in scan:
                lff.handle_peak(peak, scan.drift_time)
        lff.split_sparse(max_gap_size, min_size).smooth_overlaps(
            error_tolerance)
        self.features = lff
        return self

    def deconvolute_features(self, averagine=None, scorer=None, truncate_after=0.95,
                             minimum_intensity=5, min_size=2, max_gap_size=0.25, **kwargs):
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
            self.features, averagine, scorer, minimum_size=min_size,
            maximum_time_gap=max_gap_size)
        self.deconvoluted_features = decon.deconvolute(
            minimum_intensity=minimum_intensity, truncate_after=truncate_after, **kwargs)
        return self

    def to_deconvoluted_peak_set(self, time_bound=None):
        if self.deconvoluted_features is None:
            raise ValueError("Must first deconvolute IMS features before converting to a time-binned peak set!")
        if time_bound is not None:
            cft = IntervalTreeNode.build(
                [Interval(f.start_time, f.end_time, [f]) for f in self.deconvoluted_features])
            features = chain.from_iterable(cft.overlaps(*time_bound))
        else:
            features = self.deconvoluted_features
        return features_to_peak_set(features)


def weighted_centroid(feature):
    total = 0
    normalizer = 0
    for node in feature:
        weight = node.total_intensity()
        total += node.time * weight
        normalizer += weight
    return total / normalizer


def merge_envelopes(envelopes):
    base = envelopes[0].clone()
    for env in envelopes[1:]:
        for i, p in enumerate(env):
            base[i].intensity += p.intensity

    return base


def feature_to_peak(feature):
    peak_cluster = feature.peaks
    peak_cluster = [pi for p in peak_cluster for pi in p]
    total_intensity = sum(p.intensity for p in peak_cluster)
    mz = sum(p.mz * p.intensity for p in peak_cluster) / total_intensity
    neutral_mass = sum(
        p.neutral_mass * p.intensity for p in peak_cluster) / total_intensity
    most_abundant_mass = sum(
        p.most_abundant_mass * p.intensity for p in peak_cluster) / total_intensity
    a_to_a2_ratio = sum(
        p.a_to_a2_ratio * p.intensity for p in peak_cluster) / total_intensity
    average_mass = sum(
        p.average_mass * p.intensity for p in peak_cluster) / total_intensity
    signal_to_noise = sum(p.signal_to_noise *
                          p.intensity for p in peak_cluster) / total_intensity
    fwhm = sum(p.full_width_at_half_max *
               p.intensity for p in peak_cluster) / total_intensity
    area = sum(p.area * p.intensity for p in peak_cluster) / total_intensity
    score = sum(p.score * p.intensity for p in peak_cluster) / total_intensity
    charge = peak_cluster[0].charge
    envelope = merge_envelopes([p.envelope for p in peak_cluster])
    return IonMobilityDeconvolutedPeak(
        neutral_mass=neutral_mass, intensity=total_intensity, charge=charge, signal_to_noise=signal_to_noise,
        full_width_at_half_max=fwhm, index=-1, a_to_a2_ratio=a_to_a2_ratio, most_abundant_mass=most_abundant_mass,
        average_mass=average_mass, score=score, envelope=envelope, mz=mz, fit=None, chosen_for_msms=False,
        area=area, drift_time=weighted_centroid(feature))


def features_to_peak_set(features):
    peak_set = DeconvolutedPeakSet(tuple(map(feature_to_peak, features)))
    peak_set.reindex()
    return peak_set
