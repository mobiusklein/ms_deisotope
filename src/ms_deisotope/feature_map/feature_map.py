import logging
import array

from collections import defaultdict
from typing import List, Optional, Tuple

import numpy as np

from ms_peak_picker import PeakSet, PeakLike as PeakBase

from ms_deisotope.data_source.common import ProcessedScan
from ms_deisotope import DeconvolutedPeakSet
from ms_deisotope.peak_dependency_network.intervals import Interval, IntervalTreeNode
from ms_deisotope.peak_set import DeconvolutedPeak, IonMobilityProfileDeconvolutedPeakSolution

from .lcms_feature import LCMSFeature
from .feature_fit import DeconvolutedLCMSFeature, IonMobilityDeconvolutedLCMSFeature, IonMobilityProfileDeconvolutedLCMSFeature
from ._mz_feature_search import (
    search_sweep,
    binary_search,
    binary_search_exact,
    binary_search_with_flag,

    search_sweep_neutral,
    binary_search_neutral,
    binary_search_exact_neutral,
    binary_search_with_flag_neutral,

    _FeatureIndex,
    MZIndex,
    NeutralMassIndex)

from .feature_graph import GapAwareFeatureSmoother, GapAwareDeconvolutedFeatureSmoother, GapAwareIonMobilityDeconvolutedFeatureSmoother

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class _QueryMixin(object):
    __slots__ = ()

    def between_times(self, start: float, end: float) -> list:
        itree: IntervalTreeNode = build_rt_interval_tree(self)
        return [f for i in itree.overlaps(start, end) for f in i]


class _FeatureCollection(object):
    __slots__ = ()

    def __len__(self):
        return len(self.features)

    def __iter__(self):
        return iter(self.features)

    def __getitem__(self, i):
        if isinstance(i, (int, slice)):
            return self.features[i]
        else:
            return [self.features[j] for j in i]

    def __eq__(self, other):
        if other is None:
            return False
        return self.features == other.features

    def __ne__(self, other):
        return not self == other


class LCMSFeatureMap(_FeatureCollection, _QueryMixin):
    features: List[LCMSFeature]

    def __init__(self, features):
        self.features = sorted(features, key=lambda x: (x.mz, x.start_time))

    def search(self, mz: float, error_tolerance: Optional[float]=2e-5) -> Optional[LCMSFeature]:
        """
        Search for a single feature within `error_tolerance` of `mz`.

        Parameters
        ----------
        mz : float
            The m/z to search for
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM

        Returns
        -------
        :class:`~.LCMSFeature`
            The a feature matching the query m/z, or :const:`None`.

        See Also
        --------
        find_all
        between
        """
        i = binary_search(self.features, mz, error_tolerance)
        if i is None:
            return None
        match = self[i]
        return match

    def find_all(self, mz: float, error_tolerance: Optional[float]=2e-5) -> List[LCMSFeature]:
        """
        Search for all features within `error_tolerance` of `mz`.

        Parameters
        ----------
        mz : float
            The m/z to search for
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM

        Returns
        -------
        list of :class:`~.LCMSFeature`
            A list of all features matching the query m/z, or the empty list.

        See Also
        --------
        between
        """
        bounds = search_sweep(self.features, mz, error_tolerance)
        if bounds is not None:
            lo, hi = bounds
            return self[lo:hi]
        else:
            return []

    def spanning_time(self, time_point: float) -> List[LCMSFeature]:
        return [feature for feature in self if feature.spans_in_time(time_point)]

    def index_range(self, lo: float, hi: float, error_tolerance: Optional[float]=2e-5) -> Tuple[int, int]:
        return (
            binary_search_with_flag(
                self.features, lo, error_tolerance)[0][0],
            binary_search_with_flag(
                self.features, hi, error_tolerance)[0][0])

    def between(self, lo: float, hi: float, error_tolerance: Optional[float]=2e-5) -> List[LCMSFeature]:
        """
        Search for all features between `lo` and `hi`, allowing `error_tolerance`
        around the edges.

        Parameters
        ----------
        lo : float
            The m/z to search for the lower bound
        hi : float
            The m/z to search for the upper bound
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM

        Returns
        -------
        list of :class:`~.LCMSFeature`
            The features between the query bounds with the permitted error, or an empty list.

        See Also
        --------
        find_all
        """
        n = len(self)
        if n == 0:
            return []
        lo_ix = binary_search_with_flag(
            self.features, lo, error_tolerance)[0][0]
        if self[lo_ix].mz < lo:
            lo_ix += 1
        if lo_ix > len(self):
            lo_ix = len(self) - 1
        i = lo_ix
        while i < n:
            f = self[i]
            if f.mz > hi:
                break
            i += 1
        hi_ix = i
        return self[lo_ix:hi_ix]

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))

    def as_arrays(self):
        mz_array = array.array('d')
        intensity_array = array.array('d')
        ion_mobility_array = array.array('d')
        feature_id_array = array.array('L')
        for i, feature in enumerate(self):
            for node in feature:
                time = node.time
                for peak in node.members:
                    ion_mobility_array.append(time)
                    mz_array.append(peak.mz)
                    intensity_array.append(peak.intensity)
                    feature_id_array.append(i)
        mz_array = np.array(mz_array, copy=False)
        intensity_array = np.array(intensity_array, copy=False)
        ion_mobility_array = np.array(ion_mobility_array, copy=False)
        feature_id_array = np.array(feature_id_array, copy=False)
        mask = np.lexsort(np.stack((ion_mobility_array, mz_array)))
        return (mz_array[mask], intensity_array[mask], ion_mobility_array[mask], feature_id_array[mask])


try:
    from ms_deisotope._c.feature_map.feature_map import LCMSFeatureMap as _CLCMSFeatureMap

    class LCMSFeatureMap(_CLCMSFeatureMap, _QueryMixin):
        __slots__ = ()

except ImportError:
    pass


class LCMSFeatureForest(LCMSFeatureMap, _QueryMixin):
    """
    An an algorithm for aggregating features from peaks of similar m/z.

    Attributes
    ----------
    features : list of LCMSFeature
        A list of growing LCMSFeature objects, ordered by m/z
    count : int
        The number of peaks accumulated
    error_tolerance : float
        The mass error tolerance between peaks and possible features (in ppm)
    """

    error_tolerance: float
    count: int

    def __init__(self, features=None, error_tolerance=1e-5):
        if features is None:
            features = []
        self.features = sorted(features, key=lambda x: x.mz)
        self.error_tolerance = error_tolerance
        self.count = 0

    def find_insertion_point(self, peak: PeakBase) -> Tuple[int, bool]:
        """
        Find the position in `self` where the peak's matching feature
        should be, and whether there was a feature that matched.

        Parameters
        ----------
        peak : PeakBase
            The peak to search for

        Returns
        -------
        index : int
            The position in `self` where the feature should be.
        matched : bool
            Whether there was a feature matched at that position
        """
        index, matched = binary_search_with_flag(
            self.features, peak.mz, self.error_tolerance)
        return index, matched

    def find_minimizing_index(self, peak: PeakBase, indices: List[int]) -> int:
        """
        Amongst the set of `indices`, find the one that
        minimizes the distance to `peak`

        Parameters
        ----------
        peak : PeakBase
            The peak to match
        indices : Iterable of int
            The feature inidices in `self` to choose between

        Returns
        -------
        int:
            The best index out of `indices`, with the smallest distance
            to `peak`
        """
        best_index = None
        best_error = float('inf')
        for index_case in indices:
            feature = self[index_case]
            err = abs(feature.mz - peak.mz) / peak.mz
            if err < best_error:
                best_index = index_case
                best_error = err
        return best_index

    def handle_peak(self, peak: PeakBase, scan_time: float):
        """
        Add `peak` at `scan_time` to the feature forest.

        If `peak` matches an existing feature, insert it into that feature
        with `scan_time`. If no feature is matched, create a new feature wit
        `peak` and add it to the feature forest.

        Parameters
        ----------
        peak : PeakBase
            The peak to add.
        scan_time : float
            The time the peak was observed.

        """
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

    def insert_feature(self, feature: LCMSFeature, index: int):
        """
        Insert `feature` at `index` in :attr:`self.features`

        This method does minimal order-correctness checking.

        Parameters
        ----------
        feature : LCMSFeature
            The feature to be added.
        index : tuple of (int, bool)
            Where to insert the feature, and whether that feature had
            a match or not. By construction, index[1] should be :const`False`
            if this method is called.
        """
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

    def split_sparse(self, delta_rt=1.0, min_size=2):
        """
        Split features at gaps and eliminate features short features.

        This modifies `self` in-place.

        Parameters
        ----------
        delta_rt : float, optional
            The maximum gap size in the time dimension to tolerate. Defaults to 1.0 minute
        min_size : int, optional
            The minimum number of time points a feature must have in order to be retained.
            Defaults to 2.

        Returns
        -------
        self : LCMSFeatureForest
            This object, modified in-place.

        See Also
        --------
        LCMSFeature.split_sparse
        """
        features = []
        for f in self.features:
            for fi in f.split_sparse(delta_rt):
                if len(fi) >= min_size:
                    features.append(fi)
        self.features = features
        return self

    def smooth_overlaps(self, error_tolerance=None, time_bridge: float=0.25):
        """
        Use a single pass of :func:`smooth_overlaps` to merge features which
        are nearby in the m/z dimension and overlap in the time dimension.

        This modifies `self` in-place.

        Parameters
        ----------
        error_tolerance : float, optional
            The mass error tolerance to use when finding similar features to merge.
            Defaults to :attr:`error_tolerance`
        time_bridge : float, optional
            The maximum gap in retention time to tolerate

        Returns
        -------
        self : LCMSFeatureForest
            This object, modified in-place.

        See Also
        --------
        :func:`smooth_overlaps`
        :class:`LCMSFeatureOverlapSmoother`
        """
        if error_tolerance is None:
            error_tolerance = self.error_tolerance
        # self.features = smooth_overlaps(self.features, error_tolerance)

        self.features = GapAwareFeatureSmoother.smooth(
            self.features,
            time_bridge=time_bridge,
            mass_error_tolerance=error_tolerance).features

        return self

    def aggregate_peaks(self, scans, minimum_mz=160, minimum_intensity=500., maximum_mz=float('inf')):
        """
        Aggregate peaks from `scans` into features in this collection.

        After all peaks are added, :meth:`smooth_overlaps` is used to merge
        nearby features.

        Parameters
        ----------
        scans : :class:`Iterable`[:class:`~.ScanBase`]
            The scan collection to aggregate from
        minimum_mz : float, optional
            The minimum peak m/z to be included. Defaults to 160 m/z
        minimum_intensity : float, optional
            The minimum peak intensity to be included. Defaults to 500.
        maximum_mz : float, optional
            The maximum peak m/z to be included. Defaults to infinity.

        Returns
        -------
        self : LCMSFeatureForest
            This object, modified in-place.
        """
        for scan in scans:
            peak_set = scan.peak_set
            if peak_set is None:
                scan.pick_peaks()
                peak_set = scan.peak_set
            for peak in scan.peak_set:
                if peak.mz < minimum_mz or peak.mz > maximum_mz or peak.intensity < minimum_intensity:
                    continue
                self.handle_peak(peak, scan.scan_time)
        return self.smooth_overlaps()

    @classmethod
    def from_reader(cls, reader, error_tolerance=1e-5, minimum_mz=160, minimum_intensity=500., maximum_mz=float('inf'),
                    start_time=None, end_time=None):

        def generate():
            if start_time is not None:
                reader.start_from_scan(rt=start_time, grouped=False)
            else:
                reader.reset()
                reader.make_iterator(grouped=False)
            i = 0
            n = len(reader)
            for i, scan in enumerate(reader):
                if end_time is not None:
                    if scan.scan_time > end_time:
                        break
                logger.info("... Processed Scan %d/%d (%0.2f%%)" % (i, n, i * 100.0 / n))
                if scan.ms_level == 1:
                    yield scan

        self = cls(error_tolerance=error_tolerance)
        self.aggregate_peaks(
            generate(),
            minimum_mz=minimum_mz,
            minimum_intensity=minimum_intensity,
            maximum_mz=maximum_mz)
        return self


class DeconvolutedLCMSFeatureMap(_FeatureCollection, _QueryMixin):

    def __init__(self, features):
        self.features = sorted(features, key=lambda x: (round(x.neutral_mass, 6), x.start_time))
        self._by_mz = sorted(features, key=lambda x: (x.mz, x.start_time))

    def spanning_time(self, time_point):
        return [feature for feature in self if feature.spans_in_time(time_point)]

    def search(self, mass, error_tolerance=2e-5, use_mz=False):
        """
        Search for a single feature within `error_tolerance` of `mass`.

        By default, this searches in the neutral mass domain, to search in
        the m/z domain, pass `use_mz=True`.

        Parameters
        ----------
        mass : float
            The mass or m/z to search for
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM
        use_mz : bool
            Whether `mass` is interpreted as neutral or mass to charge ratio

        Returns
        -------
        :class:`~.LCMSFeature`
            The a feature matching the query mass, or :const:`None`.

        See Also
        --------
        find_all
        between
        """
        if use_mz:
            i = binary_search(self._by_mz, mass, error_tolerance)
            if i is None:
                return None
            return self._by_mz[i]
        i = binary_search_neutral(self.features, mass, error_tolerance)
        if i is None:
            return None
        match = self[i]
        return match

    def find_all(self, mass, error_tolerance=2e-5, use_mz=False):
        """
        Search for all features within `error_tolerance` of `mass`.

        By default, this searches in the neutral mass domain, to search in
        the m/z domain, pass `use_mz=True`.

        Parameters
        ----------
        mass : float
            The mass or m/z to search for
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM
        use_mz : bool
            Whether `mass` is interpreted as neutral or mass to charge ratio

        Returns
        -------
        list of :class:`~.LCMSFeature`
            The features matching the query mass with the permitted error, or an empty list.

        See Also
        --------
        between
        search
        """
        if use_mz:
            bounds = search_sweep(self._by_mz, mass, error_tolerance)
            collection = self._by_mz
        else:
            bounds = search_sweep_neutral(self.features, mass, error_tolerance)
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

    def between(self, lo, hi, error_tolerance=2e-5, use_mz=False):
        """
        Search for all features between `lo` and `hi`, allowing `error_tolerance`
        around the edges.

        By default, this searches in the neutral mass domain, to search in
        the m/z domain, pass `use_mz=True`.

        Parameters
        ----------
        lo : float
            The mass or m/z to search for the lower bound
        hi : float
            The mass or m/z to search for the upper bound
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM
        use_mz : bool
            Whether `lo` and `hi` are interpreted as neutral or mass to charge ratio

        Returns
        -------
        list of :class:`~.LCMSFeature`
            The features between the query bounds with the permitted error, or an empty list.

        See Also
        --------
        find_all
        """
        n = len(self)
        if n == 0:
            return []
        if use_mz:
            lo_ix = binary_search_with_flag(
                self._by_mz, lo, error_tolerance)[0][0]
            if self._by_mz[lo_ix].mz < lo:
                lo_ix += 1
            if lo_ix > len(self):
                lo_ix = len(self) - 1
            i = lo_ix
            while i < n:
                f = self._by_mz[i]
                if f.mz > hi:
                    break
                i += 1
            hi_ix = i
            return self._by_mz[lo_ix:hi_ix]
        else:
            lo_ix = binary_search_with_flag_neutral(
                self.features, lo, error_tolerance)[0][0]
            if self[lo_ix].neutral_mass < lo:
                lo_ix += 1
            if lo_ix > len(self):
                lo_ix = len(self) - 1
            i = lo_ix
            while i < n:
                f = self[i]
                if f.neutral_mass > hi:
                    break
                i += 1
            hi_ix = i
            return self[lo_ix:hi_ix]

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))

    def as_arrays(self):
        mz_array = array.array('d')
        intensity_array = array.array('d')
        charge_array = array.array('i')
        score_array = array.array('d')
        ion_mobility_array = array.array('d')
        feature_id_array = array.array('L')
        envelopes = []

        point_count = 0
        for i, feature in enumerate(self):
            for node in feature:
                time = node.time
                for peak in node.members:
                    ion_mobility_array.append(time)
                    mz_array.append(peak.mz)
                    intensity_array.append(peak.intensity)
                    score_array.append(peak.score)
                    charge_array.append(peak.charge)
                    envelopes.append(peak.envelope)
                    point_count += (len(peak.envelope) + 1) * 2
                    feature_id_array.append(i)

        mz_array = np.array(mz_array, copy=False)
        intensity_array = np.array(intensity_array, copy=False)
        charge_array = np.array(charge_array, copy=False)
        score_array = np.array(score_array, copy=False)
        ion_mobility_array = np.array(ion_mobility_array, copy=False)
        feature_id_array = np.array(feature_id_array, copy=False)
        mask = np.lexsort(np.stack((ion_mobility_array, mz_array)))
        envelope_array = np.zeros(point_count, dtype=np.float32)

        k = 0
        for j in mask:
            for point in envelopes[j]:
                envelope_array[k] = point.mz
                envelope_array[k+1] = point.intensity
                k += 2
            k += 2
        return (mz_array[mask], intensity_array[mask], charge_array[mask],
                score_array[mask], ion_mobility_array[mask], envelope_array,
                feature_id_array[mask])


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


try:
    from ms_deisotope._c.feature_map.feature_map import binary_search_with_flag
except ImportError:
    pass


class FeatureRetentionTimeInterval(Interval):
    def __init__(self, chromatogram):
        super(FeatureRetentionTimeInterval, self).__init__(
            chromatogram.start_time, chromatogram.end_time, [chromatogram])
        self.neutral_mass = chromatogram.neutral_mass
        self.start_time = self.start
        self.end_time = self.end


def build_rt_interval_tree(chromatogram_list, interval_tree_type=IntervalTreeNode):
    intervals = list(map(FeatureRetentionTimeInterval, chromatogram_list))
    interval_tree = interval_tree_type.build(intervals)
    return interval_tree


class LCMSFeatureMerger(_FeatureCollection):
    def __init__(self, features=None, error_tolerance=1e-5):
        if features is None:
            features = []
        self.features = sorted(
            features, key=self._mass_coordinate)
        self.error_tolerance = error_tolerance
        self.count = 0

    def _mass_coordinate(self, x):
        return x.mz

    # Mass Coordinate Explicit
    def find_candidates(self, new_feature):
        index, matched = binary_search_with_flag(
            self.features, new_feature.mz, self.error_tolerance)
        return index, matched

    def merge_overlaps(self, new_feature, feature_range):
        has_merged = False
        query_mass = self._mass_coordinate(new_feature)

        for chroma in feature_range:
            cond = chroma.overlaps_in_time(new_feature)
            cond = (cond and abs(
                    (self._mass_coordinate(chroma) - query_mass) / query_mass) < self.error_tolerance)
            if cond:
                chroma.merge(new_feature)
                has_merged = True
                break
        return has_merged

    # Mass Coordinate Explicit
    def find_insertion_point(self, new_feature):
        return binary_search_exact(
            self.features, new_feature.mz)

    def handle_new_feature(self, new_feature):
        if len(self) == 0:
            index = [0]
            matched = False
        else:
            index, matched = self.find_candidates(new_feature)
        if matched:

            chroma = self[index]
            has_merged = self.merge_overlaps(new_feature, chroma)
            if not has_merged:
                insertion_point = self.find_insertion_point(new_feature)
                self.insert_feature(new_feature, insertion_point)
        else:
            self.insert_feature(new_feature, index[0])
        self.count += 1

    def insert_feature(self, feature, index):
        if index != 0:
            self.features.insert(index + 1, feature)
        else:
            if len(self) == 0:
                new_index = index
            else:
                x = self.features[index]
                if self._mass_coordinate(x) < self._mass_coordinate(feature):
                    new_index = index + 1
                else:
                    new_index = index
            self.features.insert(new_index, feature)

    def aggregate_features(self, features):
        unmatched = sorted(
            features, key=lambda x: x.total_signal, reverse=True)
        for chroma in unmatched:
            self.handle_new_feature(chroma)


def flatten_tree(tree):
    """
    Perform a FIFO infix right-to-left traversal of the tree
    to flatten an :class:`~.IntervalTreeNode` hierarchy into a
    list.

    This still iterates over interval tree nodes, unlike the
    IntervalTreeNode class's own `__iter__` method, which flattens
    out the tree to iterate over contained intervals.

    Parameters
    ----------
    tree: :class:`~.IntervalTreeNode`
        The interval tree to flatten

    Returns
    -------
    list
    """
    output_queue = []
    input_queue = [tree]
    while input_queue:
        next_node = input_queue.pop()
        output_queue.append(next_node)

        next_right = next_node.right
        if next_right is not None:
            input_queue.append(next_right)

        next_left = next_node.left
        if next_left is not None:
            input_queue.append(next_left)
    return output_queue[::-1]


def layered_traversal(nodes):
    return sorted(nodes, key=lambda x: (x.level, x.center), reverse=True)


class DeconvolutedLCMSFeatureForest(DeconvolutedLCMSFeatureMap):
    """
    An algorithm for aggregating features from peaks of close mass.

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
        self.features = sorted(features, key=lambda x: x.neutral_mass)
        self._by_mz = sorted(features, key=lambda x: x.mz)
        self.error_tolerance = error_tolerance
        self.count = 0


    def find_insertion_point(self, peak: DeconvolutedPeak) -> Tuple[int, bool]:
        """
        Find the position in `self` where the peak's matching feature
        should be, and whether there was a feature that matched.

        Parameters
        ----------
        peak : DeconvolutedPeak
            The peak to search for

        Returns
        -------
        index : int
            The position in `self` where the feature should be.
        matched : bool
            Whether there was a feature matched at that position
        """
        index, matched = binary_search_with_flag_neutral(
            self.features, peak.neutral_mass, self.error_tolerance)
        return index, matched

    def find_minimizing_index(self, peak: DeconvolutedPeak, indices: List[int]) -> int:
        """
        Amongst the set of `indices`, find the one that
        minimizes the distance to `peak`

        Parameters
        ----------
        peak : PeakBase
            The peak to match
        indices : Iterable of int
            The feature inidices in `self` to choose between

        Returns
        -------
        int:
            The best index out of `indices`, with the smallest distance
            to `peak`
        """
        best_index = None
        best_error = float('inf')
        for index_case in indices:
            feature = self[index_case]
            if feature.charge != peak.charge:
                continue
            err = abs(feature.neutral_mass - peak.neutral_mass) / \
                peak.neutral_mass
            if err < best_error:
                best_index = index_case
                best_error = err
        return best_index

    def handle_peak(self, peak: DeconvolutedPeak, scan_time: float):
        """
        Add `peak` at `scan_time` to the feature forest.

        If `peak` matches an existing feature, insert it into that feature
        with `scan_time`. If no feature is matched, create a new feature wit
        `peak` and add it to the feature forest.

        Parameters
        ----------
        peak : DeconvolutedPeak
            The peak to add.
        scan_time : float
            The time the peak was observed.

        """
        if len(self) == 0:
            index = [0]
            matched = False
        else:
            index, matched = self.find_insertion_point(peak)
        if matched:
            minimized_feature_index = self.find_minimizing_index(peak, index)
            if minimized_feature_index is not None:
                feature = self.features[minimized_feature_index]
                feature.insert(peak, scan_time)
            else:
                feature = DeconvolutedLCMSFeature(charge=peak.charge)
                feature.created_at = "forest"
                feature.insert(peak, scan_time)
                self.insert_feature(feature, index)
        else:
            feature = DeconvolutedLCMSFeature(charge=peak.charge)
            feature.created_at = "forest"
            feature.insert(peak, scan_time)
            self.insert_feature(feature, index)
        self.count += 1

    def insert_feature(self, feature: DeconvolutedLCMSFeature, index: int):
        """
        Insert `feature` at `index` in :attr:`self.features`

        This method does minimal order-correctness checking.

        Parameters
        ----------
        feature : DeconvolutedLCMSFeature
            The feature to be added.
        index : tuple of (int, bool)
            Where to insert the feature, and whether that feature had
            a match or not. By construction, index[1] should be :const`False`
            if this method is called.
        """
        if index[0] != 0:
            self.features.insert(index[0] + 1, feature)
        else:
            if len(self) == 0:
                new_index = index[0]
            else:
                x = self.features[index[0]]
                if x.neutral_mass < feature.neutral_mass:
                    new_index = index[0] + 1
                else:
                    new_index = index[0]
            self.features.insert(new_index, feature)

    def split_sparse(self, delta_rt=1.0, min_size=2):
        """
        Split features at gaps and eliminate features short features.

        This modifies `self` in-place.

        Parameters
        ----------
        delta_rt : float, optional
            The maximum gap size in the time dimension to tolerate. Defaults to 1.0 minute
        min_size : int, optional
            The minimum number of time points a feature must have in order to be retained.
            Defaults to 2.

        Returns
        -------
        self : DeconvolutedLCMSFeatureForest
            This object, modified in-place.

        See Also
        --------
        DeconvolutedLCMSFeature.split_sparse
        """
        features = [fi for f in self.features for fi in f.split_sparse(
            delta_rt) if len(fi) >= min_size]
        self.features = features
        return self

    def smooth_overlaps(self, error_tolerance=None):
        """
        Use a single pass of :func:`smooth_overlaps_neutral` to merge features which
        are nearby in the neutral mass dimension and overlap in the time dimension.

        This modifies `self` in-place.

        Parameters
        ----------
        error_tolerance : float, optional
            The mass error tolerance to use when finding similar features to merge.
            Defaults to :attr:`error_tolerance`

        Returns
        -------
        self : DeconvolutedLCMSFeatureForest
            This object, modified in-place.

        See Also
        --------
        :func:`smooth_overlaps_neutral`
        :class:`NeutralMassLCMSFeatureOverlapSmoother`
        """
        if error_tolerance is None:
            error_tolerance = self.error_tolerance

        self.features = GapAwareDeconvolutedFeatureSmoother.smooth(self.features, mass_error_tolerance=error_tolerance)
        self._by_mz = sorted(self.features, key=lambda x: x.mz)
        return self

    def aggregate_peaks(self, scans, minimum_mass=160, minimum_intensity=10., maximum_mass=float('inf')):
        """
        Aggregate peaks from `scans` into features in this collection.

        After all peaks are added, :meth:`smooth_overlaps` is used to merge
        nearby features.

        Parameters
        ----------
        scans : :class:`Iterable`[:class:`~.ScanBase`]
            The scan collection to aggregate from
        minimum_mass : float, optional
            The minimum peak neutral mass to be included. Defaults to 160 Da
        minimum_intensity : float, optional
            The minimum peak intensity to be included. Defaults to 10.
        maximum_mass : float, optional
            The maximum peak neutral mass to be included. Defaults to infinity.

        Returns
        -------
        self : DeconvolutedLCMSFeatureForest
            This object, modified in-place.
        """
        for scan in scans:
            peak_set = scan.deconvoluted_peak_set
            if peak_set is None:
                scan.pick_peaks()
                peak_set = scan.deconvoluted_peak_set
            for peak in peak_set:
                if peak.neutral_mass < minimum_mass or peak.neutral_mass > maximum_mass or peak.intensity < minimum_intensity:
                    continue
                self.handle_peak(peak, scan.scan_time)
        self.smooth_overlaps()
        return self

    @classmethod
    def from_reader(cls, reader, error_tolerance=1e-5, minimum_mass=160, minimum_intensity=10., maximum_mass=float('inf'),
                    start_time=None, end_time=None, ms_level=1):

        def generate():
            if start_time is not None:
                reader.start_from_scan(rt=start_time, grouped=False)
            else:
                reader.reset()
                reader.make_iterator(grouped=False)
            i = 0
            n = len(reader)
            for i, scan in enumerate(reader):
                if end_time is not None:
                    if scan.scan_time > end_time:
                        break
                    logger.info("... Processed Scan %d/%d (%0.2f%%)" %
                                (i, n, i * 100.0 / n))
                if scan.ms_level == ms_level:
                    yield scan

        self = cls(error_tolerance=error_tolerance)
        self.aggregate_peaks(
            generate(),
            minimum_mass=minimum_mass,
            minimum_intensity=minimum_intensity,
            maximum_mass=maximum_mass)
        return self


class IonMobilityDeconvolutedLCMSFeatureForest(DeconvolutedLCMSFeatureForest):
    def __init__(self, features=None, error_tolerance=1e-5, drift_error_tolerance=0.1):
        self.drift_error_tolerance = drift_error_tolerance
        super(IonMobilityDeconvolutedLCMSFeatureForest, self).__init__(features, error_tolerance)

    def smooth_overlaps(self, error_tolerance=None):
        if error_tolerance is None:
            error_tolerance = self.error_tolerance
        self.features = GapAwareIonMobilityDeconvolutedFeatureSmoother.smooth(
            self.features, error_tolerance, ion_mobility_error_tolerance=self.drift_error_tolerance).features
        return self

    def find_minimizing_index(self, peak, indices):
        best_index = None
        best_error = float('inf')
        for index_case in indices:
            feature = self[index_case]
            if feature.charge != peak.charge:
                continue
            if abs(feature.drift_time - peak.drift_time) >= self.drift_error_tolerance:
                continue
            err = abs(feature.neutral_mass - peak.neutral_mass) / peak.neutral_mass
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
            minimized_feature_index = self.find_minimizing_index(peak, index)
            if minimized_feature_index is not None:
                feature = self.features[minimized_feature_index]
                feature.insert(peak, scan_time)
            else:
                feature = IonMobilityDeconvolutedLCMSFeature(charge=peak.charge)
                feature.created_at = "forest"
                feature.insert(peak, scan_time)
                self.insert_feature(feature, index)
        else:
            feature = IonMobilityDeconvolutedLCMSFeature(charge=peak.charge)
            feature.created_at = "forest"
            feature.insert(peak, scan_time)
            self.insert_feature(feature, index)
        self.count += 1

    @classmethod
    def from_reader(cls, reader, error_tolerance=1e-5, drift_error_tolerance=0.1, minimum_mass=160,
                    minimum_intensity=10., maximum_mass=float('inf'), start_time=None, end_time=None,
                    ms_level=1):

        def generate():
            if start_time is not None:
                reader.start_from_scan(rt=start_time, grouped=False)
            else:
                reader.reset()
                reader.make_iterator(grouped=False)
            i = 0
            n = len(reader)
            for i, scan in enumerate(reader):
                if end_time is not None:
                    if scan.scan_time > end_time:
                        break
                logger.info("... Processed Scan %d/%d (%0.2f%%)" % (i, n, i * 100.0 / n))
                if scan.ms_level == ms_level:
                    yield scan

        self = cls(error_tolerance=error_tolerance,
                   drift_error_tolerance=drift_error_tolerance)
        self.aggregate_peaks(
            generate(),
            minimum_mass=minimum_mass,
            minimum_intensity=minimum_intensity,
            maximum_mass=maximum_mass)
        return self


class IonMobilityProfileDeconvolutedLCMSFeatureForest(DeconvolutedLCMSFeatureForest):
    def find_minimizing_index(self, peak: IonMobilityProfileDeconvolutedPeakSolution, indices):
        """
        Amongst the set of `indices`, find the one that
        minimizes the distance to `peak`

        Parameters
        ----------
        peak : PeakBase
            The peak to match
        indices : Iterable of int
            The feature inidices in `self` to choose between

        Returns
        -------
        int:
            The best index out of `indices`, with the smallest distance
            to `peak`
        """
        best_index = None
        best_error = float('inf')
        for index_case in indices:
            feature = self[index_case]
            if feature.charge != peak.charge:
                continue
            if not feature.ion_mobility_interval.overlaps(peak.ion_mobility_interval):
                continue
            err = abs(feature.neutral_mass - peak.neutral_mass) / \
                peak.neutral_mass
            if err < best_error and err < self.error_tolerance:
                best_index = index_case
                best_error = err
        return best_index

    def handle_peak(self, peak, scan_time):
        """
        Add `peak` at `scan_time` to the feature forest.

        If `peak` matches an existing feature, insert it into that feature
        with `scan_time`. If no feature is matched, create a new feature wit
        `peak` and add it to the feature forest.

        Parameters
        ----------
        peak : DeconvolutedPeak
            The peak to add.
        scan_time : float
            The time the peak was observed.

        """
        if len(self) == 0:
            index = [0]
            matched = False
        else:
            index, matched = self.find_insertion_point(peak)
        if matched:
            minimized_feature_index = self.find_minimizing_index(peak, index)
            if minimized_feature_index is not None:
                feature = self.features[minimized_feature_index]
                feature.insert(peak, scan_time)
            else:
                feature = IonMobilityProfileDeconvolutedLCMSFeature(charge=peak.charge)
                feature.created_at = "forest"
                feature.insert(peak, scan_time)
                self.insert_feature(feature, index)
        else:
            feature = IonMobilityProfileDeconvolutedLCMSFeature(charge=peak.charge)
            feature.created_at = "forest"
            feature.insert(peak, scan_time)
            self.insert_feature(feature, index)
        self.count += 1

    def smooth_overlaps(self, error_tolerance=None):
        from .feature_graph import GapAwareIonMobilityProfileDeconvolutedFeatureSmoother
        if error_tolerance is None:
            error_tolerance = self.error_tolerance
        self.features = GapAwareIonMobilityProfileDeconvolutedFeatureSmoother.smooth(
            self.features, mass_error_tolerance=error_tolerance).features
        return self


def remove_peaks_below_threshold(feature_map: _FeatureCollection, minimum_intensity: float, min_size: int = 2, maximum_time_gap: float = 1.0) -> _FeatureCollection:
    out = []
    feature: LCMSFeature
    for feature in feature_map:
        keep = []
        for node in feature:
            if node.total_intensity() > minimum_intensity:
                keep.append(node)

        if keep:
            filtered_feature: LCMSFeature = feature._copy_chunk(
                keep)
            filtered_feature.feature_id = feature.feature_id
            out.extend(filter(
                lambda f: len(f) >= min_size,
                filtered_feature.split_sparse(
                    delta_rt=maximum_time_gap)))
    feature_map = feature_map.__class__(out)
    return feature_map
