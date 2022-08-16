
from collections import deque
from typing import Any, Deque, Iterator, List, Optional, Set, Type, Sequence

import numpy as np

from ms_deisotope.peak_dependency_network.intervals import (
    SpanningMixin,
    IntervalTreeNode,
    Interval,
    SimpleInterval,
    IntervalTreeNode2D
)

from ms_deisotope.feature_map.lcms_feature import LCMSFeature
from ._mz_feature_search import MZIndex, NeutralMassIndex


def mz_point_organizer_callback(contained_intervals):
    return IntervalTreeNode.build(
        [Interval(node.mz, node.mz, [node]) for node in contained_intervals])


def neutral_mass_point_organizer_callback(contained_intervals):
    return IntervalTreeNode.build(
        [Interval(node.neutral_mass, node.neutral_mass, [node]) for node in contained_intervals])


class FeatureGraphNode(SpanningMixin):
    feature: LCMSFeature
    index: int
    edges: Set['FeatureGraphEdge']
    center: float
    mz: float

    def __init__(self, feature: LCMSFeature, index: int, edges: Optional[Set]=None):
        self.start = feature.start_time
        self.end = feature.end_time
        if edges is None:
            edges = set()
        self.feature = feature
        self.index = index
        self.edges = edges

        total = 0
        abundance = 0
        for node in self.feature.nodes:
            intensity = node.total_intensity()
            total += node.time * intensity
            abundance += intensity
        self.center = total / abundance
        self.mz = feature.mz

    @property
    def start_time(self) -> float:
        return self.start

    @property
    def end_time(self) -> float:
        return self.end

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.feature,)

    def __index__(self):
        return self.index

    def __hash__(self):
        return hash(self.index)

    def __eq__(self, other):
        return self.feature == other.feature


class DeconvolutedFeatureGraphNode(FeatureGraphNode):
    neutral_mass: float
    charge: int

    def __init__(self, feature, index, edges=None):
        super().__init__(feature, index, edges)
        self.neutral_mass = feature.neutral_mass
        self.charge = feature.charge


class IonMobilityDeconvolutedFeatureGraphNode(DeconvolutedFeatureGraphNode):
    ion_mobility: float

    def __init__(self, feature, index, edges=None):
        super().__init__(feature, index, edges)
        self.ion_mobility = feature.drift_time


class IonMobilityProfileFeatureGraphNode(DeconvolutedFeatureGraphNode):
    ion_mobility_interval: SimpleInterval

    def __init__(self, feature, index, edges=None):
        super().__init__(feature, index, edges)
        self.ion_mobility_interval = feature.ion_mobility_interval


class FeatureGraphEdge(object):
    __slots__ = ('node_a', 'node_b', 'transition', 'weight', 'mass_error', 'rt_error', '_hash')

    node_a: FeatureGraphNode
    node_b: FeatureGraphNode
    transition: Any
    weight: float
    mass_error: float
    rt_error: float
    _hash: int

    def __init__(self, node_a, node_b, transition, weight=1, mass_error=0, rt_error=0):
        self.node_a = node_a
        self.node_b = node_b
        self.transition = transition
        self.weight = weight
        self.mass_error = mass_error
        self.rt_error = rt_error
        if self.node_a.index < self.node_b.index:
            self._hash = hash(((self.node_a.index, self.node_b.index)))
        else:
            self._hash = hash(((self.node_b.index, self.node_a.index)))
        self.node_a.edges.add(self)
        self.node_b.edges.add(self)

    def __repr__(self):
        return "%s(%s, %s, %s)" % (
            self.__class__.__name__, self.node_a.feature, self.node_b.feature,
            self.transition)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return self.node_a.index == other.node_a.index and self.node_b.index == other.node_b.index

    def __ne__(self, other):
        return not self == other

    def remove(self):
        self.node_a.edges.remove(self)
        self.node_b.edges.remove(self)


class TimeQuery(SpanningMixin):
    __slots__ = ()

    def __init__(self, feature: LCMSFeature, width: float = 0):
        self.start = feature.start_time - width
        self.end = feature.end_time + width

    def __repr__(self):
        return "TimeQuery(%f, %f)" % (self.start, self.end)


class PPMQuery(SpanningMixin):
    __slots__ = ()

    def __init__(self, center, error_tolerance):
        self.start = center - center * error_tolerance
        self.end = center + center * error_tolerance

    def __repr__(self):
        return "PPMQuery(%f, %f)" % (self.start, self.end)


try:
    from ms_deisotope._c.feature_map.feature_graph import (
        FeatureGraphNode, DeconvolutedFeatureGraphNode, IonMobilityProfileFeatureGraphNode,
        FeatureGraphEdge)
except ImportError:
    raise


class FeatureGraph(MZIndex):
    node_cls: Type = FeatureGraphNode
    edge_cls: Type = FeatureGraphEdge

    _features: Sequence[LCMSFeature]
    nodes: List[FeatureGraphNode]
    assigned_seed_queue: Deque[FeatureGraphNode]
    rt_tree: IntervalTreeNode2D
    edges: Set[FeatureGraphEdge]

    @property
    def features(self):
        return self.nodes

    def __init__(self, features):
        self._features = features
        self.assigned_seed_queue = deque()
        self.nodes = self._construct_graph_nodes(self._features)
        self.edges = set()
        self.rt_tree = None
        self._make_rt_tree()

    def _make_rt_tree(self):
        self.rt_tree = IntervalTreeNode2D.build(
            self.nodes, mz_point_organizer_callback)

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return iter(self.nodes)

    def __getitem__(self, i):
        return self.nodes[i]

    def _construct_graph_nodes(self, features):
        nodes = []
        for i, chroma in enumerate(features):
            node = (self.node_cls(chroma, i))
            nodes.append(node)
            self.enqueue_seed(node)
        return nodes

    def enqueue_seed(self, feature):
        self.assigned_seed_queue.append(feature)

    def pop_seed(self) -> FeatureGraphNode:
        return self.assigned_seed_queue.popleft()

    def iterseeds(self) -> Iterator[FeatureGraphNode]:
        while self.assigned_seed_queue:
            yield self.pop_seed()

    def find_edges(self, node: FeatureGraphNode, query_width: float = 2., error_tolerance=1.5e-5, **kwargs):
        query = TimeQuery(node.feature, query_width)
        mass_query = PPMQuery(node.mz, error_tolerance)
        nodes = self.rt_tree.overlaps_2d(
            np.array([query.start, mass_query.start]), np.array([query.end, mass_query.end]))

        for match in nodes:
            ppm_error = (node.mz - match.mz) / match.mz
            if abs(ppm_error) > error_tolerance:
                continue
            if match.index == node.index:
                continue
            rt_error = (node.center - match.center)
            self.edges.add(self.edge_cls(
                node, match, None, mass_error=ppm_error, rt_error=rt_error))

    def build(self, query_width: float = 2., error_tolerance=1.5e-5, **kwargs):
        for node in self.iterseeds():
            self.find_edges(node, query_width=query_width,
                            error_tolerance=error_tolerance, **kwargs)

    def adjacency_matrix(self):
        adjmat = np.zeros((len(self), len(self)))
        nodes = self.nodes
        for node in nodes:
            for edge in node.edges:
                adjmat[edge.node_a.index, edge.node_b.index] = 1
        return adjmat

    def connected_components(self) -> List[List[FeatureGraphNode]]:
        pool = set(self.nodes)
        components = []

        i = 0
        while pool:
            i += 1
            current_component = set()
            visited = set()
            node = pool.pop()
            current_component.add(node)
            j = 0
            while current_component:
                j += 1
                node = current_component.pop()
                visited.add(node)
                for edge in node.edges:
                    if edge.node_a not in visited and edge.node_a in pool:
                        current_component.add(edge.node_a)
                        pool.remove(edge.node_a)
                    if edge.node_b not in visited and edge.node_b in pool:
                        current_component.add(edge.node_b)
                        pool.remove(edge.node_b)
            components.append(list(visited))
        return components


class DeconvolutedFeatureGraph(FeatureGraph, NeutralMassIndex):
    node_cls = DeconvolutedFeatureGraphNode

    def _make_rt_tree(self):
        self.rt_tree = IntervalTreeNode2D.build(
            self.nodes, neutral_mass_point_organizer_callback)

    def find_all(self, mass: float, error_tolerance: float=0.00002):
        return NeutralMassIndex.find_all(self, mass, error_tolerance)

    def find_edges(self, node: DeconvolutedFeatureGraphNode, query_width: float = 2., error_tolerance: float=1.5e-5, **kwargs):
        query = TimeQuery(node.feature, query_width)
        mass_query = PPMQuery(node.neutral_mass, error_tolerance)
        nodes = self.rt_tree.overlaps_2d(
            np.array([query.start, mass_query.start]), np.array([query.end, mass_query.end]))
        charge = node.charge
        for match in nodes:
            if match.charge != charge:
                continue
            ppm_error = (node.neutral_mass - match.neutral_mass) / match.neutral_mass
            if abs(ppm_error) > error_tolerance:
                continue
            if match.index == node.index:
                continue
            rt_error = (node.center - match.center)
            self.edges.add(self.edge_cls(
                node, match, None, mass_error=ppm_error, rt_error=rt_error))


class IonMobilityDeconvolutedFeatureGraph(DeconvolutedFeatureGraph):
    node_cls = IonMobilityDeconvolutedFeatureGraphNode

    def find_edges(self, node: IonMobilityDeconvolutedFeatureGraphNode, query_width: float = 2., error_tolerance: float = 1.5e-5, ion_mobility_error_tolerance: float = 0.01, **kwargs):
        query = TimeQuery(node.feature, query_width)
        mass_query = PPMQuery(node.neutral_mass, error_tolerance)
        nodes = self.rt_tree.overlaps_2d(
            np.array([query.start, mass_query.start]), np.array([query.end, mass_query.end]))
        charge = node.charge
        for match in nodes:
            if match.charge != charge:
                continue
            # If there isn't an overlap in the IM dimension, then they aren't the same molecule even
            # if charge, mass, and elution time do match
            if abs(match.ion_mobility - node.ion_mobility) > ion_mobility_error_tolerance:
                continue
            ppm_error = (node.neutral_mass - match.neutral_mass) / \
                match.neutral_mass
            if abs(ppm_error) > error_tolerance:
                continue
            if match.index == node.index:
                continue
            rt_error = (node.center - match.center)
            self.edges.add(self.edge_cls(
                node, match, None, mass_error=ppm_error, rt_error=rt_error))


class IonMobilityProfileDeconvolutedFeatureGraph(DeconvolutedFeatureGraph):
    node_cls = IonMobilityProfileFeatureGraphNode

    def find_edges(self, node: IonMobilityProfileFeatureGraphNode, query_width: float = 2., error_tolerance: float=1.5e-5, **kwargs):
        query = TimeQuery(node.feature, query_width)
        mass_query = PPMQuery(node.neutral_mass, error_tolerance)
        nodes = self.rt_tree.overlaps_2d(
            np.array([query.start, mass_query.start]), np.array([query.end, mass_query.end]))
        charge = node.charge
        for match in nodes:
            if match.charge != charge:
                continue
            # If there isn't an overlap in the IM dimension, then they aren't the same molecule even
            # if charge, mass, and elution time do match
            if not match.ion_mobility_interval.overlaps(node.ion_mobility_interval):
                continue
            ppm_error = (node.neutral_mass - match.neutral_mass) / \
                match.neutral_mass
            if abs(ppm_error) > error_tolerance:
                continue
            if match.index == node.index:
                continue
            rt_error = (node.center - match.center)
            self.edges.add(self.edge_cls(
                node, match, None, mass_error=ppm_error, rt_error=rt_error))


try:
    from ms_deisotope._c.feature_map.feature_graph import (
        find_edges, find_edges_neutral_mass, find_edges_neutral_mass_ion_mobility, connected_components)

    FeatureGraph.find_edges = find_edges
    FeatureGraph.connected_components = connected_components
    DeconvolutedFeatureGraph.find_edges = find_edges_neutral_mass
    IonMobilityProfileDeconvolutedFeatureGraph.find_edges = find_edges_neutral_mass_ion_mobility
except ImportError:
    raise


class GapAwareFeatureSmoother(object):
    graph: FeatureGraph

    def __init__(self, features):
        self.graph = FeatureGraph(features)

    def _wrap_result(self, features) -> 'LCMSFeatureMap':
        from ms_deisotope.feature_map.feature_map import LCMSFeatureMap
        return LCMSFeatureMap(features)

    def connect_components(self, time_bridge: float=0.5, mass_error_tolerance: float=1.5e-5):
        self.graph.build(query_width=time_bridge, error_tolerance=mass_error_tolerance)
        features: List[LCMSFeature] = []
        for component in self.graph.connected_components():
            f = component[0].feature
            for f2 in component[1:]:
                f = f.merge(f2.feature)
            features.append(f)
        return self._wrap_result(features)

    @classmethod
    def smooth(cls, features, time_bridge: float=0.5, mass_error_tolerance: float=1.5e-5, **kwargs):
        return cls(features).connect_components(time_bridge=time_bridge, mass_error_tolerance=mass_error_tolerance, **kwargs)


class GapAwareDeconvolutedFeatureSmoother(GapAwareFeatureSmoother):
    graph: DeconvolutedFeatureGraph

    def __init__(self, features):
        self.graph = DeconvolutedFeatureGraph(features)

    def _wrap_result(self, features) -> 'DeconvolutedLCMSFeatureMap':
        from ms_deisotope.feature_map.feature_map import DeconvolutedLCMSFeatureMap
        return DeconvolutedLCMSFeatureMap(features)


class GapAwareIonMobilityDeconvolutedFeatureSmoother(GapAwareDeconvolutedFeatureSmoother):
    graph: IonMobilityDeconvolutedFeatureGraph

    def __init__(self, features):
        self.graph = IonMobilityDeconvolutedFeatureGraph(features)

    def connect_components(self, time_bridge: float = 0.5, mass_error_tolerance: float = 1.5e-5, ion_mobility_error_tolerance: float=0.1):
        self.graph.build(query_width=time_bridge,
                         error_tolerance=mass_error_tolerance,
                         ion_mobility_error_tolerance=ion_mobility_error_tolerance)
        features: List[LCMSFeature] = []
        for component in self.graph.connected_components():
            f = component[0].feature
            for f2 in component[1:]:
                f = f.merge(f2.feature)
            features.append(f)
        return self._wrap_result(features)

    def _wrap_result(self, features) -> 'IonMobilityDeconvolutedLCMSFeatureMap':
        from ms_deisotope.feature_map.feature_map import IonMobilityDeconvolutedLCMSFeatureForest
        return IonMobilityDeconvolutedLCMSFeatureForest(features)


class GapAwareIonMobilityProfileDeconvolutedFeatureSmoother(GapAwareDeconvolutedFeatureSmoother):
    graph: IonMobilityProfileDeconvolutedFeatureGraph

    def __init__(self, features):
        self.graph = IonMobilityProfileDeconvolutedFeatureGraph(features)
