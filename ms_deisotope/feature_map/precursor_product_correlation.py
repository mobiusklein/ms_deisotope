import numpy as np

from ms_deisotope.task import LogUtilsMixin
from ms_deisotope.data_source.scan.base import ScanBunch

from ms_deisotope.peak_set import IonMobilityProfileDeconvolutedPeakSolution, DeconvolutedPeakSet
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


class PseudoXIC(object):
    __slots__ = ("x", "y", "source")

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
        return np.interp([t], self.x, self.y, left=0, right=0)[0]

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
                   match: IonMobilityProfileDeconvolutedLCMSFeature) -> IonMobilityProfileDeconvolutedLCMSFeature:
    submatch = match.between_ion_mobilities(
        feature.ion_mobility_interval.start,
        feature.ion_mobility_interval.end,
        feature.start_time,
        feature.end_time)
    return submatch


class PrecursorProductCorrelationGraph(LogUtilsMixin):
    edge_cls = FeatureGraphEdge

    precursor_graph: IonMobilityProfileDeconvolutedFeatureGraph
    product_graph: IonMobilityProfileDeconvolutedFeatureGraph
    edges: set
    _is_built: bool

    def __init__(self, precursor_graph, product_graph):
        self.precursor_graph = precursor_graph
        self.product_graph = product_graph
        self.edges = set()
        self._is_built = False

    def score_edge(self, feature: IonMobilityProfileDeconvolutedLCMSFeature,
                   match: IonMobilityProfileDeconvolutedLCMSFeature) -> float:
        submatch: IonMobilityProfileDeconvolutedLCMSFeature = match.between_ion_mobilities(
            feature.ion_mobility_interval.start,
            feature.ion_mobility_interval.end,
            feature.start_time, feature.end_time)
        if len(submatch) < 2:
            return -1.0
        _x, y1, y2 = features_to_aligned_arrays(feature, submatch)
        corr = np.corrcoef(
            sliding_mean(y1),
            sliding_mean(y2),
        )[0, 1]
        if np.isnan(corr):
            corr = -1.0
        return corr

    def find_edges(self, precursor_node: IonMobilityProfileDeconvolutedFeatureGraph.node_cls, minimum_score=0.1, **kwargs):
        nodes = self.product_graph.rt_tree.overlaps(
            [precursor_node.start_time, 0.0],
            [precursor_node.end_time, precursor_node.neutral_mass]
        )
        charge = precursor_node.charge
        for match in nodes:
            if match.charge > charge:
                continue
            if not match.ion_mobility_interval.overlaps(precursor_node.ion_mobility_interval):
                continue
            score = self.score_edge(precursor_node.feature, match.feature)
            if score < minimum_score:
                continue
            self.edges.add(self.edge_cls(
                precursor_node, match, score))

    def build(self, min_charge=2, min_size=3):
        n = len(self.precursor_graph)
        for i, precursor in enumerate(self.precursor_graph):
            if precursor.charge < min_charge:
                continue
            if len(precursor.feature) < min_size:
                continue
            if i % 1000 == 0:
                self.log(f"... Processed {i}/{n} precursors ({i / n * 100.0:0.2f}%)")
            self.find_edges(precursor)
        self._is_built = True

    def iterspectra(self, **kwargs):
        if not self._is_built:
            self.build(**kwargs)
        iterator = PrecursorProductCorrelatingIterator(self.precursor_graph, self.product_graph)
        for batch in iterator:
            yield batch


def edge_to_pseudopeak(edge: FeatureGraphEdge):
    precursor = edge.node_a.feature
    product = edge.node_b.feature
    sub_product = reduce_feature(precursor, product)
    peak = IonMobilityProfileDeconvolutedPeakSolution.from_feature(sub_product)
    peak.intensity *= edge.transition ** 2
    return peak


def get_weighted_output_xic(edge: FeatureGraphEdge, rt: float, return_weight=False):
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


class XICPseudoSpectrumGenerator(SpanningMixin):
    __slots__ = ('node', 'edge_to_pseudo_xic', 'edge_to_pseudo_xic_children', 'minimum_intensity', 'precursor_information')

    def __init__(self, node: IonMobilityProfileFeatureGraphNode, minimum_intensity: float=1):
        self.node = node
        self.edge_to_pseudo_xic = {}
        self.edge_to_pseudo_xic_children = {}
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
        if edge in self.edge_to_pseudo_xic:
            return self.edge_to_pseudo_xic[edge]
        f = reduce_feature(self.node.feature, edge.node_a.feature)
        if f:
            xic = PseudoXIC.from_feature(f)
        else:
            xic = None
        self.edge_to_pseudo_xic[edge] = xic
        return xic

    def get_xic_for_edge_child(self, edge: FeatureGraphEdge) -> PseudoXIC:
        if edge in self.edge_to_pseudo_xic_children:
            return self.edge_to_pseudo_xic_children[edge]
        f = reduce_feature(self.node.feature, edge.node_b.feature)
        if f:
            xic = PseudoXIC.from_feature(f)
        else:
            breakpoint()
            xic = None
        self.edge_to_pseudo_xic_children[edge] = xic
        return xic

    def get_weighted_output(self, edge: FeatureGraphEdge, time: float) -> float:
        weights = []
        index = None
        o: FeatureGraphEdge
        for i, o in enumerate(edge.node_b.edges):
            if o.node_a.index == self.node.index:
                index = i
            f = self.get_xic_for_edge(o)
            if not f:
                weights.append(0)
            else:
                v = f.interpolate(time)
                weights.append(v * o.transition ** 2)

        weights = np.array(weights)
        weights /= weights.max()
        product = self.get_xic_for_edge_child(edge)
        return (weights[index] * product.interpolate(time) * edge.transition ** 2)

    def pseudopeak_for(self, edge: FeatureGraphEdge, time: float) -> IonMobilityProfileDeconvolutedPeakSolution:
        weight = self.get_weighted_output(edge, time)
        xic = self.get_xic_for_edge_child(edge)
        f = xic.source
        peak = IonMobilityProfileDeconvolutedPeakSolution.from_feature(f)
        peak.intensity = weight
        return peak

    def pseudospectrum_for(self, time: float) -> DeconvolutedPeakSet:
        peaks = []
        for edge in self.node.edges:
            p = self.pseudopeak_for(edge, time)
            if p.intensity > self.minimum_intensity:
                peaks.append(p)
        peaks = DeconvolutedPeakSet(peaks)
        peaks.reindex()
        return peaks


class PrecursorProductCorrelatingIterator(LogUtilsMixin):
    precursor_graph: IonMobilityProfileDeconvolutedFeatureGraph
    product_graph: IonMobilityProfileDeconvolutedFeatureGraph

    def __init__(self, precursor_graph, product_graph):
        self.precursor_graph = precursor_graph
        self.product_graph = product_graph
        self.scan_counter = 0
        self.ms1_times = None
        self.delta_ms1_time = None
        self._prepass_ms1()

    def _prepass_ms1(self):
        times = set()
        for node in self.precursor_graph:
            node.generator = XICPseudoSpectrumGenerator(node)
            times.update(node.feature.times)
        self.ms1_times = np.array(sorted(times))
        self.delta_ms1_time = np.mean(np.diff(self.ms1_times))

    def __iter__(self):
        i = 0

        for i in range(len(self.ms1_times)):
            time = self.ms1_times[i]
            try:
                next_time = self.ms1_times[i + 1]
            except IndexError:
                next_time = time + self.delta_ms1_time

            peaks_for_scan = []
            msn_spectra = []

            self.log(
                f"... Processing time point {time} with {self.scan_counter} scans produced ({i / len(self.ms1_times) * 100.0:0.2f}%)")

            nodes_to_clear = []
            for node_k in self.precursor_graph.rt_tree.contains_point(time):
                cn, _ = node_k.feature.find_time(time)
                if not cn:
                    continue
                peaks_for_scan.extend(cn.members)
                msn_for = node_k.generator.pseudospectrum_for(time)
                if msn_for:
                    pinfo = node_k.generator.precursor_information.copy()
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

            for node_k in nodes_to_clear:
                node_k.generator = None

            yield ScanBunch(ms1_scan, msn_scans)
