import logging

from collections import defaultdict
from typing import DefaultDict, Dict, Iterator, List, Sized, Union
from ms_deisotope.data_source.metadata.activation import ActivationInformation
from ms_deisotope.data_source.scan.mobility_frame import IonMobilityFrame, IonMobilitySourceRandomAccessFrameSource
from ms_deisotope.data_source.scan.scan import ProcessedScan, WrappedScan
from ms_deisotope.data_source.metadata.scan_traits import ion_mobility_drift_time

import numpy as np

from ms_deisotope.task import LogUtilsMixin
from ms_deisotope.data_source.scan.base import ScanBase, ScanBunch

from ms_deisotope.peak_set import IonMobilityProfileDeconvolutedPeakSolution, DeconvolutedPeakSet, IonMobilityDeconvolutedPeak
from ms_deisotope.peak_dependency_network.intervals import SpanningMixin

from ms_deisotope.data_source import PrecursorInformation, make_scan

from ms_deisotope.feature_map.profile_transform import sliding_mean

from ms_deisotope.feature_map.feature_fit import (
    IonMobilityProfileDeconvolutedLCMSFeature,
    DeconvolutedLCMSFeature)

from ms_deisotope.feature_map.feature_graph import (
    FeatureGraphEdge,
    IonMobilityProfileDeconvolutedFeatureGraph,
    IonMobilityProfileFeatureGraphNode)


from ms_deisotope._c.utils import correlation
from ms_deisotope._c.feature_map.profile_transform import interpolate_point

logger = logging.getLogger(__name__)


class PseudoXIC(object):
    __slots__ = ("x", "y", "source", )

    x: np.ndarray
    y: np.ndarray
    source: IonMobilityProfileDeconvolutedLCMSFeature

    def __init__(self, x: np.ndarray, y: np.ndarray, source: IonMobilityProfileDeconvolutedLCMSFeature):
        self.x = x
        self.y = y
        self.source = source

    def smooth(self) -> "PseudoXIC":
        return self.__class__(self.x, sliding_mean(self.y), self.source)

    def interpolate(self, t: float) -> float:
        return interpolate_point(self.x, self.y, t)

    @classmethod
    def from_feature(cls, source: IonMobilityProfileDeconvolutedLCMSFeature, smooth: bool = True) -> 'PseudoXIC':
        x, y = source.as_arrays(np.float32)
        if smooth:
            y = sliding_mean(y)
        return cls(x, y, source)

    def as_arrays(self):
        return self.x, self.y

    @property
    def neutral_mass(self) -> float:
        return self.source.neutral_mass

    @property
    def charge(self) -> int:
        return self.source.charge

    @property
    def mz(self) -> float:
        return self.source.mz

    @property
    def start_time(self) -> float:
        return self.x[0]

    @property
    def end_time(self) -> float:
        return self.x[-1]


def interpolate_feature(f: DeconvolutedLCMSFeature, Xp: np.ndarray) -> np.ndarray:
    x, y = f.as_arrays(np.float32)
    return np.interp(Xp, x, y, 0, 0)


def features_to_aligned_arrays(f1: DeconvolutedLCMSFeature, f2: DeconvolutedLCMSFeature, points_per_minute: int=150):
    start = f1.start_time
    end = f1.end_time
    n_points = int((end - start) * points_per_minute)
    X = np.linspace(start, end, n_points)
    Y1 = interpolate_feature(f1, X)
    Y2 = interpolate_feature(f2, X)
    return X, Y1, Y2


def reduce_feature(feature: IonMobilityProfileDeconvolutedLCMSFeature,
                   match: IonMobilityProfileDeconvolutedLCMSFeature,
                   include_mz: bool = True, include_envelope: bool = True) -> IonMobilityProfileDeconvolutedLCMSFeature:
    submatch = match.between_ion_mobilities(
        feature.ion_mobility_interval.start,
        feature.ion_mobility_interval.end,
        feature.start_time,
        feature.end_time,
        include_mz=include_mz,
        include_envelope=include_envelope)
    return submatch


class PrecursorProductCorrelationGraph(LogUtilsMixin):
    edge_cls = FeatureGraphEdge

    precursor_graph: IonMobilityProfileDeconvolutedFeatureGraph
    product_graph: IonMobilityProfileDeconvolutedFeatureGraph
    max_edge_count_per_node: int
    edges: set
    _is_built: bool

    def __init__(self, precursor_graph, product_graph, max_edge_count_per_node=1500):
        self.precursor_graph = precursor_graph
        self.product_graph = product_graph
        self.max_edge_count_per_node = max_edge_count_per_node
        self.edges = set()
        self._is_built = False

    def score_edge(self, feature: IonMobilityProfileDeconvolutedLCMSFeature,
                   match: IonMobilityProfileDeconvolutedLCMSFeature,
                   minimum_overlap_size: int=2) -> float:
        submatch: IonMobilityProfileDeconvolutedLCMSFeature = match.between_ion_mobilities(
            feature.ion_mobility_interval.start,
            feature.ion_mobility_interval.end,
            feature.start_time, feature.end_time,
            include_envelope=False,
            include_mz=False)
        if len(submatch) < minimum_overlap_size:
            return -1.0
        _x, y1, y2 = features_to_aligned_arrays(feature, submatch)
        corr = correlation(
            sliding_mean(y1),
            sliding_mean(y2),
        )
        if np.isnan(corr):
            corr = -1.0
        return corr

    def find_edges(self, precursor_node: IonMobilityProfileDeconvolutedFeatureGraph.node_cls, minimum_score=0.1, **kwargs):
        nodes = self.product_graph.rt_tree.overlaps_2d(
            np.array([precursor_node.start_time, 0.0]),
            np.array([precursor_node.end_time, precursor_node.neutral_mass])
        )
        charge = precursor_node.charge
        edges = set()
        for match in nodes:
            if match.charge > charge:
                continue
            if not match.ion_mobility_interval.overlaps(precursor_node.ion_mobility_interval):
                continue
            score = self.score_edge(precursor_node.feature, match.feature)
            if score < minimum_score:
                continue
            edge = self.edge_cls(
                precursor_node, match, score)
            edges.add(edge)
        edges = sorted(edges, key=lambda x: x.transition, reverse=True)
        keepers = edges[:self.max_edge_count_per_node]
        rest = edges[self.max_edge_count_per_node:]
        for edge in rest:
            edge.node_a.edges.remove(edge)
            edge.node_b.edges.remove(edge)
        self.edges.update(keepers)


    def build(self, min_charge=2, min_size=3):
        n = len(self.precursor_graph)
        last_logged = 0
        batch_size = 500
        for i, precursor in enumerate(self.precursor_graph):
            if precursor.charge < min_charge:
                continue
            if len(precursor.feature) < min_size:
                continue
            if i % 500 == 0 and i or i - last_logged >= batch_size:
                last_logged = i
                self.log(f"... Processed {i}/{n} precursors ({i / n * 100.0:0.2f}%)")
            self.find_edges(precursor)
        self._is_built = True

    def iterspectra(self, weight_scale_factor: float=1.0, **kwargs):
        if not self._is_built:
            self.build(**kwargs)
        iterator = PrecursorProductCorrelatingIterator(self.precursor_graph, self.product_graph, weight_scale_factor=weight_scale_factor)
        for batch in iterator:
            yield batch


PrecursorProductCorrelationGraph.log_with_logger(logger)


def edge_to_pseudopeak(edge: FeatureGraphEdge) -> IonMobilityProfileDeconvolutedPeakSolution:
    precursor = edge.node_a.feature
    product = edge.node_b.feature
    sub_product = reduce_feature(precursor, product)
    peak = IonMobilityProfileDeconvolutedPeakSolution.from_feature(sub_product)
    peak.intensity *= edge.transition ** 2
    return peak


def get_weighted_output_xic(edge: FeatureGraphEdge, rt: float, return_weight=False, weight_scale_factor: float=1.0):
    weights = []
    index = None
    o: FeatureGraphEdge
    for i, o in enumerate(edge.node_b.edges):
        if o.node_a.index == edge.node_a.index:
            index = i
        f = reduce_feature(edge.node_a.feature, o.node_a.feature)
        if not f:
            weights.append(0)
            continue
        f = PseudoXIC.from_feature(f)
        v = f.interpolate(rt)
        weights.append(v * o.transition)
    weights = np.array(weights)
    weights /= weights.max()
    if weight_scale_factor != 1.0:
        weights *= weight_scale_factor
    product = reduce_feature(edge.node_a.feature, edge.node_b.feature)
    f = PseudoXIC.from_feature(product)
    if return_weight:
        return (weights[index], edge.transition ** 2, weights[index] * f.interpolate(rt) * edge.transition ** 2, )
    return (weights[index] * f.interpolate(rt) * edge.transition ** 2)


def edge_at_time_to_pseudopeak(edge: FeatureGraphEdge, rt: float):
    precursor = edge.node_a.feature
    product = edge.node_b.feature
    sub_product = reduce_feature(precursor, product)
    peak = IonMobilityProfileDeconvolutedPeakSolution.from_feature(sub_product)
    peak.intensity = get_weighted_output_xic(edge, rt)
    return peak


class EdgeToXICCache(object):
    store: Dict[FeatureGraphEdge, PseudoXIC]
    ref_count: DefaultDict[FeatureGraphEdge, int]

    def __init__(self):
        self.store = dict()
        self.ref_count = defaultdict(int)

    def __getitem__(self, key):
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[key] = value

    def __contains__(self, key):
        return key in self.store

    def __len__(self):
        return len(self.store)

    def checkout(self, key: FeatureGraphEdge):
        value = self.store[key]
        self.ref_count[key] += 1
        return value

    def checkin(self, key: FeatureGraphEdge):
        # Ref-count may be negative if someone
        # belatedly checks in an edge that has
        # already been evicted directly.
        count = self.ref_count[key]
        count -= 1
        if count <= 0:
            self.store.pop(key, None)
            self.ref_count.pop(key)
            return True
        else:
            self.ref_count[key] = count
        return False

    def __delitem__(self, key):
        self.store.pop(key, None)
        self.ref_count.pop(key, None)

    def clear(self):
        self.store.clear()
        self.ref_count.clear()

    def session(self) -> 'EdgeToXICCacheSession':
        return EdgeToXICCacheSession(self)


class EdgeToXICCacheSession(object):
    cache: EdgeToXICCache
    store: Dict[FeatureGraphEdge, PseudoXIC]

    def __init__(self, cache: EdgeToXICCache):
        self.cache = cache
        self.store = {}

    def __contains__(self, key):
        return key in self.cache

    def __getitem__(self, key: FeatureGraphEdge) -> PseudoXIC:
        try:
            return self.store[key]
        except KeyError:
            value = self.cache.checkout(key)
            self.store[key] = value
            return value

    def __setitem__(self, key, value):
        self.cache[key] = value
        self.store[key] = self.cache.checkout(key)

    def __delitem__(self, key):
        del self.store[key]
        self.cache.checkin(key)

    def remove(self, key: FeatureGraphEdge):
        if key in self.store:
            del self[key]

    def clear(self):
        for key in self.store:
            self.cache.checkin(key)
        self.store.clear()

    def __del__(self):
        self.clear()

    def __len__(self):
        return len(self.store)


class XICPseudoSpectrumGenerator(SpanningMixin):
    __slots__ = ('node', 'edge_to_pseudo_xic', 'edge_to_pseudo_xic_children', 'minimum_intensity', 'precursor_information')

    node: IonMobilityProfileFeatureGraphNode
    minimum_intensity: float
    precursor_information: PrecursorInformation
    edge_to_pseudo_xic: Union[EdgeToXICCacheSession, Dict[FeatureGraphEdge, PseudoXIC]]
    edge_to_pseudo_xic_children: Union[EdgeToXICCacheSession, Dict[FeatureGraphEdge, PseudoXIC]]

    def __init__(self, node: IonMobilityProfileFeatureGraphNode, minimum_intensity: float = 1, edge_to_pseudo_xic=None,
                 edge_to_pseudo_xic_children=None):
        if edge_to_pseudo_xic is None:
            edge_to_pseudo_xic = {}
        if edge_to_pseudo_xic_children is None:
            edge_to_pseudo_xic_children = {}

        self.node = node
        self.edge_to_pseudo_xic = edge_to_pseudo_xic
        self.edge_to_pseudo_xic_children = edge_to_pseudo_xic_children
        self.minimum_intensity = minimum_intensity
        self.start = self.node.start
        self.end = self.node.end
        self.precursor_information = None
        self._make_precursor_information()

    def _make_precursor_information(self):
        self.precursor_information = PrecursorInformation(
            self.node.feature.mz,
            self.node.feature.intensity,
            self.node.charge)

        dt = (self.node.ion_mobility_interval.end +
              self.node.ion_mobility_interval.start) / 2
        self.precursor_information._ion_mobility.add_ion_mobility(
            'ion mobility drift time', dt)

    def get_xic_for_edge(self, edge: FeatureGraphEdge) -> PseudoXIC:
        try:
            return self.edge_to_pseudo_xic[edge]
        except KeyError:
            f = reduce_feature(self.node.feature, edge.node_a.feature, False, False)
            if f:
                xic = PseudoXIC.from_feature(f)
            else:
                xic = None
            self.edge_to_pseudo_xic[edge] = xic
            return xic

    def get_xic_for_edge_child(self, edge: FeatureGraphEdge) -> PseudoXIC:
        try:
            return self.edge_to_pseudo_xic_children[edge]
        except KeyError:
            # Need to keep envelope and m/z as the product XIC's source will be used
            # to generate pseudo-peaks
            f = reduce_feature(self.node.feature, edge.node_b.feature, True, True)
            if f:
                xic = PseudoXIC.from_feature(f)
            else:
                xic = None
            self.edge_to_pseudo_xic_children[edge] = xic
            return xic

    def clear(self):
        self.edge_to_pseudo_xic.clear()
        self.edge_to_pseudo_xic_children.clear()

    def get_weighted_output(self, edge: FeatureGraphEdge, time: float, weight_scale_factor: float=1.0) -> float:
        weights = []
        index = None
        o: FeatureGraphEdge
        step_back = 0
        for i, o in enumerate(edge.node_b.edges):
            if o.node_a.index == self.node.index:
                index = i
            # If a parent node of a sharer of `edge` here doesn't span `time`,
            # evict it from the session, and if the index of `self.node` hasn't
            # been found yet, increment the backtrack counter by one before
            # skipping the interpolation.
            if not o.node_a.contains(time):
                self.edge_to_pseudo_xic.remove(o)
                if index is None:
                    step_back += 1
                continue
            f = self.get_xic_for_edge(o)
            if f is None:
                weights.append(0)
            else:
                v = f.interpolate(time)
                weights.append(v * o.transition ** 2)
        index -= step_back
        weights = np.array(weights)
        weights /= weights.max()
        if weight_scale_factor != 1.0:
            weights *= weight_scale_factor
        product = self.get_xic_for_edge_child(edge)
        return (weights[index] * product.interpolate(time) * edge.transition ** 2)

    def pseudopeak_for(self, edge: FeatureGraphEdge, time: float, weight_scale_factor: float=1.0) -> IonMobilityProfileDeconvolutedPeakSolution:
        weight = self.get_weighted_output(
            edge, time, weight_scale_factor=weight_scale_factor)
        xic = self.get_xic_for_edge_child(edge)
        f = xic.source
        peak = IonMobilityProfileDeconvolutedPeakSolution.from_feature(f)
        peak.intensity = weight
        return peak

    def pseudospectrum_for(self, time: float, weight_scale_factor: float = 1.0) -> DeconvolutedPeakSet:
        peaks = []
        for edge in self.node.edges:
            p = self.pseudopeak_for(
                edge, time, weight_scale_factor=weight_scale_factor)
            if p.intensity > self.minimum_intensity:
                peaks.append(p)
        peaks = DeconvolutedPeakSet(peaks)
        peaks.reindex()
        return peaks


class PrecursorProductCorrelatingIterator(LogUtilsMixin):
    precursor_graph: IonMobilityProfileDeconvolutedFeatureGraph
    product_graph: IonMobilityProfileDeconvolutedFeatureGraph
    weight_scale_factor: float

    edge_to_pseudo_xic: EdgeToXICCache
    edge_to_pseudo_xic_children: EdgeToXICCache

    def __init__(self, precursor_graph, product_graph, weight_scale_factor: float=1.0):
        self.precursor_graph = precursor_graph
        self.product_graph = product_graph
        self.scan_counter = 0
        self.ms1_times = None
        self.delta_ms1_time = None
        self.edge_to_pseudo_xic = EdgeToXICCache()
        self.edge_to_pseudo_xic_children = EdgeToXICCache()
        self._prepass_ms1()
        self.weight_scale_factor = weight_scale_factor

    def _prepass_ms1(self):
        times = set()
        for node in self.precursor_graph:
            node.generator = XICPseudoSpectrumGenerator(
                node, 1, self.edge_to_pseudo_xic.session(),
                self.edge_to_pseudo_xic_children.session())
            times.update(node.feature.times)
        self.ms1_times = np.array(sorted(times))
        self.delta_ms1_time = np.mean(np.diff(self.ms1_times))

    def __iter__(self):
        i = 0

        nodes_to_clear = []
        for i in range(len(self.ms1_times)):
            time = self.ms1_times[i]
            try:
                next_time = self.ms1_times[i + 1]
            except IndexError:
                next_time = time + self.delta_ms1_time

            peaks_for_scan = []
            msn_spectra = []

            for node_k in nodes_to_clear:
                node_k.generator.clear()
                node_k.generator = None
                for edge in list(node_k.edges):
                    del self.edge_to_pseudo_xic[edge]
                    edge.remove()

            nodes_to_clear = []
            self.log(
                f"... Processing time point {time:0.3f} with {self.scan_counter} scans produced, set size: \
{len(self.edge_to_pseudo_xic) + len(self.edge_to_pseudo_xic_children)} (\
{i / len(self.ms1_times) * 100.0:0.2f}%)")

            for node_k in self.precursor_graph.rt_tree.contains_point(time):
                cn, _ = node_k.feature.find_time(time)
                if not cn:
                    continue
                peaks_for_scan.extend(cn.members)
                generator: XICPseudoSpectrumGenerator = node_k.generator
                msn_for = generator.pseudospectrum_for(
                    time, weight_scale_factor=self.weight_scale_factor)
                if msn_for:
                    pinfo = generator.precursor_information.copy()
                    pinfo.intensity = cn.total_intensity()
                    msn_spectra.append(
                        (msn_for, pinfo)
                    )
                # Need to check if this is the last time point for this precursor node
                # so it can have its generator destroyed.
                if abs(node_k.end - time) < 1e-3:
                    nodes_to_clear.append(node_k)

            msn_spectra = [(t, ) + p for t, p in zip(np.linspace(time,
                                                                next_time, len(msn_spectra) + 2)[1:-1], msn_spectra)]
            peaks_for_scan = DeconvolutedPeakSet(peaks_for_scan)
            peaks_for_scan.reindex()

            ms1_scan = make_scan(
                id=f"merged={self.scan_counter}",
                ms_level=1, scan_time=time,
                index=self.scan_counter,
                deconvoluted_peak_set=peaks_for_scan)

            self.scan_counter += 1
            msn_scans = []
            pinfo: PrecursorInformation
            for msn_time, peaks, pinfo in msn_spectra:
                pinfo.precursor_scan_id = ms1_scan.id
                msn_scans.append(make_scan(
                    id=f"merged={self.scan_counter}",
                    ms_level=2, scan_time=msn_time,
                    index=self.scan_counter,
                    deconvoluted_peak_set=peaks,
                    precursor_information=pinfo))
                self.scan_counter += 1


            yield ScanBunch(ms1_scan, msn_scans)


PrecursorProductCorrelatingIterator.log_with_logger(logger)


class NaiveIonMobilityOverlapIterator(LogUtilsMixin):
    frame_pair_iterator: Union[Iterator[ScanBunch],
                               IonMobilitySourceRandomAccessFrameSource]

    def __init__(self, frame_pair_iterator):
        self.frame_pair_iterator = frame_pair_iterator
        self.scan_counter = 0

    def make_ms1_scan(self, precursor: IonMobilityFrame) -> WrappedScan:
        peaks_for_scan = [
            IonMobilityDeconvolutedPeak(
                f.neutral_mass, f.intensity, f.charge, f.intensity, 0,
                0, 0, 0, 0,
                score=f.score,
                envelope=f.sum_envelopes(),
                mz=f.mz,
                drift_time=f.apex_time)
            for f in
            precursor.deconvoluted_features
        ]

        ms1_scan = make_scan(
            id=precursor.id,
            ms_level=1, scan_time=precursor.time,
            index=self.scan_counter,
            deconvoluted_peak_set=DeconvolutedPeakSet(peaks_for_scan))

        self.scan_counter += 1
        return ms1_scan

    def __len__(self):
        return len(self.frame_pair_iterator)

    def split_products(self, ms1_scan: ScanBase, precursor: IonMobilityFrame, product: IonMobilityFrame) -> List[WrappedScan]:
        msn_scans = []
        for precursor_feature in precursor.deconvoluted_features:
            peaks = product.to_deconvoluted_peak_set(
                (precursor_feature.start_time,
                 precursor_feature.end_time),
                precursor_feature.neutral_mass)
            if not peaks:
                continue
            pinfo = PrecursorInformation(
                precursor_feature.mz,
                precursor_feature.intensity,
                precursor_feature.charge,
                ms1_scan.id,
                None,
                precursor_feature.neutral_mass,
                precursor_feature.charge,
                precursor_feature.intensity)
            pinfo._ion_mobility.add_ion_mobility(
                ion_mobility_drift_time,
                precursor_feature.apex_time
            )
            scan = make_scan(
                id=f"merged={self.scan_counter}",
                ms_level=2, scan_time=product.time,
                index=self.scan_counter,
                activation=ActivationInformation('cid', 2),
                deconvoluted_peak_set=peaks,
                precursor_information=pinfo)

            scan.annotations['source frame'] = product.id
            self.scan_counter += 1
            msn_scans.append(scan)
        return msn_scans

    def __iter__(self):
        return self

    def __next__(self):
        bunch: ScanBunch = next(self.frame_pair_iterator)
        while bunch.precursor is None:
            bunch: ScanBunch = next(self.frame_pair_iterator)
        ms1_scan = self.make_ms1_scan(bunch.precursor)
        msn_scans = []
        for product in bunch.products:
            msn_scans.extend(self.split_products(ms1_scan, bunch.precursor, product))
        return ScanBunch(ms1_scan, msn_scans)
