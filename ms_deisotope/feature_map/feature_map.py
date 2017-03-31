from collections import defaultdict
from .lcms_feature import LCMSFeature
from ms_deisotope.data_source.common import ProcessedScan
from ms_deisotope import DeconvolutedPeakSet
from ms_peak_picker import PeakSet


class LCMSFeatureMap(object):

    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.mz)

    def __len__(self):
        return len(self.features)

    def __iter__(self):
        return iter(self.features)

    def __getitem__(self, i):
        if isinstance(i, (int, slice)):
            return self.features[i]
        else:
            return [self.features[j] for j in i]

    def search(self, mz, error_tolerance=2e-5):
        i = binary_search(self.features, mz, error_tolerance)
        if i is None:
            return None
        match = self[i]
        return match

    def find_all(self, mz, error_tolerance=2e-5):
        bounds = search_sweep(self.features, mz, error_tolerance)
        if bounds is not None:
            lo, hi = bounds
            return self[lo:hi]
        else:
            return []

    def index_range(self, lo, hi, error_tolerance=2e-5):
        return (
            binary_search_with_flag(
                self.features, lo, error_tolerance)[0][0],
            binary_search_with_flag(
                self.features, hi, error_tolerance)[0][0])

    def between(self, lo, hi, error_tolerance):
        lo_ix = binary_search_with_flag(
            self.features, lo, error_tolerance)[0][0]
        hi_ix = binary_search_with_flag(
            self.features, hi, error_tolerance)[0][0]
        if self[lo_ix].mz < lo:
            lo_ix += 1
        if self[hi_ix] > hi:
            hi_ix -= 1
        if hi_ix < 0:
            hi_ix = 0
        if lo_ix > len(self):
            lo_ix = len(self) - 1
        return self[lo_ix:hi_ix]

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))


try:
    from ms_deisotope._c.feature_map.feature_map import LCMSFeatureMap
except ImportError:
    pass


class LCMSFeatureForest(LCMSFeatureMap):
    """An an algorithm for aggregating features from peaks of close mass
    weighted by intensity.

    This algorithm assumes that mass accuracy is correlated with intensity, so
    the most intense peaks should most accurately reflect their true neutral mass.
    The expected input is a list of (scan id, peak) pairs. This list is sorted by
    descending peak intensity. For each pair, using binary search, locate the nearest
    existing feature in :attr:`features`. If the nearest feature is within
    :attr:`error_tolerance` ppm of the peak's neutral mass, add this peak to that
    feature, otherwise create a new feature containing this peak and insert
    it into :attr:`features` while preserving the overall sortedness. This algorithm
    is carried out by :meth:`aggregate_unmatched_peaks`

    This process may produce features with large gaps in them, which
    may or may not be acceptable. To break gapped features into separate
    entities, the :class:`LCMSFeatureFilter` type has a method :meth:`split_sparse`.

    Attributes
    ----------
    features : list of LCMSFeature
        A list of growing LCMSFeature objects, ordered by neutral mass
    count : int
        The number of peaks accumulated
    error_tolerance : float
        The mass error tolerance between peaks and possible features (in ppm)
    """
    def __init__(self, features=None, error_tolerance=1e-5):
        if features is None:
            features = []
        self.features = sorted(features, key=lambda x: x.mz)
        self.error_tolerance = error_tolerance
        self.count = 0

    def find_insertion_point(self, peak):
        index, matched = binary_search_with_flag(
            self.features, peak.mz, self.error_tolerance)
        return index, matched

    def find_minimizing_index(self, peak, indices):
        best_index = None
        best_error = float('inf')
        for index_case in indices:
            feature = self[index_case]
            err = abs(feature.mz - peak.mz) / peak.mz
            if err < best_error:
                best_index = index_case
                best_error = err
        return best_index

    def handle_peak(self, peak, scan_time):
        if len(self) == 0:
            index = [0]
            matched = False
        else:
            index, matched = self.find_insertion_point(peak)
        if matched:
            feature = self.features[self.find_minimizing_index(peak, index)]
            feature.insert(peak, scan_time)
        else:
            feature = LCMSFeature()
            feature.created_at = "forest"
            feature.insert(peak, scan_time)
            self.insert_feature(feature, index)
        self.count += 1

    def insert_feature(self, feature, index):
        if index[0] != 0:
            self.features.insert(index[0] + 1, feature)
        else:
            if len(self) == 0:
                new_index = index[0]
            else:
                x = self.features[index[0]]
                if x.mz < feature.mz:
                    new_index = index[0] + 1
                else:
                    new_index = index[0]
            self.features.insert(new_index, feature)

    def aggregate_peaks(self, scans, minimum_mz=160, minimum_intensity=500.):
        for scan in scans:
            for peak in scan.peak_set:
                if peak.mz < minimum_mz or peak.intensity < minimum_intensity:
                    continue
                self.handle_peak(peak, scan.scan_time)
        self.features = smooth_overlaps(self.features, self.error_tolerance)


def smooth_overlaps(feature_list, error_tolerance=1e-5):
    feature_list = sorted(feature_list, key=lambda x: x.mz)
    out = []
    last = feature_list[0]
    i = 1
    while i < len(feature_list):
        current = feature_list[i]
        mass_error = abs((last.mz - current.mz) / current.mz)
        if mass_error <= error_tolerance:
            if last.overlaps_in_time(current):
                last = last.merge(current)
                last.created_at = "smooth_overlaps"
            else:
                out.append(last)
                last = current
        else:
            out.append(last)
            last = current
        i += 1
    out.append(last)
    return out


def binary_search_with_flag(array, mz, error_tolerance=1e-5):
        lo = 0
        n = hi = len(array)
        while hi != lo:
            mid = (hi + lo) / 2
            x = array[mid]
            err = (x.mz - mz) / mz
            if abs(err) <= error_tolerance:
                i = mid - 1
                # Begin Sweep forward
                while i > 0:
                    x = array[i]
                    err = (x.mz - mz) / mz
                    if abs(err) <= error_tolerance:
                        i -= 1
                        continue
                    else:
                        break
                low_end = i
                i = mid + 1

                # Begin Sweep backward
                while i < n:
                    x = array[i]
                    err = (x.mz - mz) / mz
                    if abs(err) <= error_tolerance:
                        i += 1
                        continue
                    else:
                        break
                high_end = i
                return list(range(low_end, high_end)), True
            elif (hi - lo) == 1:
                return [mid], False
            elif err > 0:
                hi = mid
            elif err < 0:
                lo = mid
        return 0, False


def binary_search(array, mz, error_tolerance=1e-5):
    """Binary search an ordered array of objects with :attr:`mz`
    using a PPM error tolerance of `error_toler

    Parameters
    ----------
    array : list
        An list of objects, sorted over :attr:`mz` in increasing order
    mz : float
        The mz to search for
    error_tolerance : float, optional
        The PPM error tolerance to use when deciding whether a match has been found

    Returns
    -------
    int:
        The index in `array` of the best match
    bool:
        Whether or not a match was actually found, used to
        signal behavior to the caller.
    """
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) / 2
        x = array[mid]
        err = (x.mz - mz) / mz
        if abs(err) <= error_tolerance:
            best_index = mid
            best_error = err
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err < best_error:
                    best_error = err
                    best_index = i
                i -= 1

            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err < best_error:
                    best_error = err
                    best_index = i
                i += 1
            return best_index
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


def search_sweep(array, mz, error_tolerance=1e-5):
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) / 2
        x = array[mid]
        err = (x.mz - mz) / mz
        if abs(err) <= error_tolerance:
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err > error_tolerance:
                    break
                i -= 1
            start = i + 1
            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err > error_tolerance:
                    break
                i += 1
            end = i
            return (start, end)
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


def binary_search_with_flag_neutral(array, neutral_mass, error_tolerance=1e-5):
        lo = 0
        n = hi = len(array)
        while hi != lo:
            mid = (hi + lo) / 2
            x = array[mid]
            err = (x.neutral_mass - neutral_mass) / neutral_mass
            if abs(err) <= error_tolerance:
                i = mid - 1
                # Begin Sweep forward
                while i > 0:
                    x = array[i]
                    err = (x.neutral_mass - neutral_mass) / neutral_mass
                    if abs(err) <= error_tolerance:
                        i -= 1
                        continue
                    else:
                        break
                low_end = i
                i = mid + 1

                # Begin Sweep backward
                while i < n:
                    x = array[i]
                    err = (x.neutral_mass - neutral_mass) / neutral_mass
                    if abs(err) <= error_tolerance:
                        i += 1
                        continue
                    else:
                        break
                high_end = i
                return list(range(low_end, high_end)), True
            elif (hi - lo) == 1:
                return [mid], False
            elif err > 0:
                hi = mid
            elif err < 0:
                lo = mid
        return 0, False


def binary_search_neutral(array, neutral_mass, error_tolerance=1e-5):
    """Binary search an ordered array of objects with :attr:`neutral_mass`
    using a PPM error tolerance of `error_toler

    Parameters
    ----------
    array : list
        An list of objects, sorted over :attr:`neutral_mass` in increasing order
    neutral_mass : float
        The neutral_mass to search for
    error_tolerance : float, optional
        The PPM error tolerance to use when deciding whether a match has been found

    Returns
    -------
    int:
        The index in `array` of the best match
    bool:
        Whether or not a match was actually found, used to
        signal behavior to the caller.
    """
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) / 2
        x = array[mid]
        err = (x.neutral_mass - neutral_mass) / neutral_mass
        if abs(err) <= error_tolerance:
            best_index = mid
            best_error = err
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err < best_error:
                    best_error = err
                    best_index = i
                i -= 1

            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err < best_error:
                    best_error = err
                    best_index = i
                i += 1
            return best_index
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


def search_sweep_neutral(array, neutral_mass, error_tolerance=1e-5):
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) / 2
        x = array[mid]
        err = (x.neutral_mass - neutral_mass) / neutral_mass
        if abs(err) <= error_tolerance:
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err > error_tolerance:
                    break
                i -= 1
            start = i + 1
            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err > error_tolerance:
                    break
                i += 1
            end = i
            return (start, end)
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


class DeconvolutedLCMSFeatureMap(object):

    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.neutral_mass)
        self._by_mz = sorted(features, key=lambda x: x.mz)

    def __len__(self):
        return len(self.features)

    def __iter__(self):
        return iter(self.features)

    def __getitem__(self, i):
        if isinstance(i, (int, slice)):
            return self.features[i]
        else:
            return [self.features[j] for j in i]

    def search(self, mz, error_tolerance=2e-5, use_mz=False):
        if use_mz:
            i = binary_search(self._by_mz, mz, error_tolerance)
            if i is None:
                return None
            return self._by_mz[i]
        i = binary_search_neutral(self.features, mz, error_tolerance)
        if i is None:
            return None
        match = self[i]
        return match

    def find_all(self, mz, error_tolerance=2e-5, use_mz=False):
        if use_mz:
            bounds = search_sweep(self._by_mz, mz, error_tolerance)
            collection = self._by_mz
        else:
            bounds = search_sweep_neutral(self.features, mz, error_tolerance)
            collection = self
        if bounds is not None:
            lo, hi = bounds
            return collection[lo:hi]
        else:
            return []

    def index_range(self, lo, hi, error_tolerance=2e-5):
        return (
            binary_search_with_flag_neutral(
                self.features, lo, error_tolerance)[0][0],
            binary_search_with_flag_neutral(
                self.features, hi, error_tolerance)[0][0])

    def between(self, lo, hi, error_tolerance):
        lo_ix = binary_search_with_flag_neutral(
            self.features, lo, error_tolerance)[0][0]
        hi_ix = binary_search_with_flag_neutral(
            self.features, hi, error_tolerance)[0][0]
        if self[lo_ix].mz < lo:
            lo_ix += 1
        if self[hi_ix] > hi:
            hi_ix -= 1
        if hi_ix < 0:
            hi_ix = 0
        if lo_ix > len(self):
            lo_ix = len(self) - 1
        return self[lo_ix:hi_ix]

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))


def convert_map_to_scan_peak_list(feature_map, peak_loader, time_precision=4, deconvoluted=True):
    metadata_map = {}
    scan_accumulator = defaultdict(list)
    for scan_id, metadata in peak_loader.extended_index.ms1_ids.items():
        metadata_map[round(metadata["scan_time"], time_precision)] = metadata
    for feature in feature_map:
        for node in feature:
            scan_accumulator[round(node.time, time_precision)].extend(node.members)

    packed = []
    for key, peaks in sorted(scan_accumulator.items(), key=lambda x: x[0]):
        template = peak_loader.get_scan_by_time(key)
        if deconvoluted:
            peak_set = PeakSet([])
            deconvoluted_peak_set = DeconvolutedPeakSet(peaks)
        else:
            peak_set = PeakSet(peaks)
            deconvoluted_peak_set = DeconvolutedPeakSet([])
        peak_set.reindex()
        deconvoluted_peak_set.reindex()
        scan = ProcessedScan(
            template.id, template.title, None, template.ms_level, template.scan_time,
            template.index, peak_set, deconvoluted_peak_set, template.polarity,
            None)
        packed.append(scan)
    return packed


class _FeatureIndex(object):
    def __iter__(self):
        return iter(self.features)

    def __getitem__(self, i):
        return self.features[i]

    def __len__(self):
        return len(self.features)

    def __nonzero__(self):
        return bool(self.features)

    def __bool__(self):
        return bool(self.features)

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))


class MZIndex(_FeatureIndex):
    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.mz)

    def find_all(self, mz, error_tolerance=2e-5):
        bounds = search_sweep(self.features, mz, error_tolerance)
        if bounds is not None:
            lo, hi = bounds
            return self[lo:hi]
        else:
            return []


class NeutralMassIndex(_FeatureIndex):
    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.neutral_mass)

    def find_all(self, mass, error_tolerance=2e-5):
        bounds = search_sweep_neutral(self.features, mass, error_tolerance)
        if bounds is not None:
            lo, hi = bounds
            return self[lo:hi]
        else:
            return []
