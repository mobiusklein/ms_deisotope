
from collections import deque
from typing import Any, Deque, Iterator, List, Optional, Set, Type, Sequence

import numpy as np

from ms_deisotope.peak_dependency_network.intervals import SpanningMixin, IntervalTreeNode, Interval, SimpleInterval
from ms_deisotope.feature_map.lcms_feature import LCMSFeature


class IntervalTreeNode2D(IntervalTreeNode):
    def __init__(self, center, left, contained, right, level=0, parent=None, inner_organizer=None):
        super().__init__(center, left, contained, right, level=level, parent=parent)
        self.inner_organizer = inner_organizer
        self.organize()

    def set_organizer(self, organizer):
        self.inner_organizer = organizer
        if self.left is not None:
            self.left.set_organizer(organizer)
        if self.right is not None:
            self.right.set_organizer(organizer)
        self.organize()

    def organize(self):
        if self.inner_organizer is not None:
            self.organized = self.inner_organizer(self.contained)
        else:
            self.organized = None

    def _overlaps_interval(self, start, end):
        starts = start
        ends = end
        query = Interval(starts[0], ends[0])
        result = []
        if self.organized is None:
            return super()._overlaps_interval(starts[0], ends[0])
        else:
            for interv2 in self.organized.overlaps(starts[1], ends[1]):
                for interv in interv2.members:
                    if interv.overlaps(query):
                        result.append(interv)
        return result

    def overlaps(self, start, end):
        starts = start
        ends = end
        start = starts[0]
        end = ends[0]
        result = []
        if start > self.end:
            return result
        if end < self.start:
            return result
        elif start <= self.start:
            if end < self.start:
                return result
            else:
                if self.left is not None:
                    result.extend(self.left.overlaps(starts, ends))
                result.extend(self._overlaps_interval(starts, ends))
                if self.right is not None and end >= self.right.start:
                    result.extend(self.right.overlaps(starts, ends))
        elif start > self.start:
            if self.left is not None and self.left.end >= start:
                result.extend(self.left.overlaps(starts, ends))
            result.extend(self._overlaps_interval(starts, ends))
            if self.right is not None and end >= self.right.start:
                result.extend(self.right.overlaps(starts, ends))
        elif end > self.start:
            if self.left is not None:
                result.extend(self.left.overlaps(starts, ends))
            result.extend(self._overlaps_interval(starts, ends))
            if self.right is not None and end >= self.right.start:
                result.extend(self.right.overlaps(starts, ends))
        return result

    @classmethod
    def build(cls, intervals, organizer_callback=None):
        root: IntervalTreeNode2D = super(IntervalTreeNode2D, cls).build(intervals)
        root.set_organizer(organizer_callback)
        return root


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
        self.start_time = self.start = feature.start_time
        self.end_time = self.end = feature.end_time
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
        self._hash = None
        self.node_a.edges.add(self)
        self.node_b.edges.add(self)

    def __repr__(self):
        return "%s(%s, %s, %s)" % (
            self.__class__.__name__, self.node_a.feature, self.node_b.feature,
            self.transition)

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(
                frozenset((self.node_a.index, self.node_b.index)))
        return self._hash

    def __eq__(self, other):
        return (self.node_a.index, self.node_b.index) == (other.node_a.index, other.node_b.index)

    def __ne__(self, other):
        return not self == other


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


class FeatureGraph(object):
    node_cls: Type = FeatureGraphNode
    edge_cls: Type = FeatureGraphEdge

    features: Sequence[LCMSFeature]
    nodes: List[FeatureGraphNode]
    assigned_seed_queue: Deque[FeatureGraphNode]
    rt_tree: IntervalTreeNode2D
    edges: Set[FeatureGraphEdge]

    def __init__(self, features):
        self.features = features
        self.assigned_seed_queue = deque()
        self.nodes = self._construct_graph_nodes(self.features)
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
        nodes = self.rt_tree.overlaps(
            [query.start, mass_query.start], [query.end, mass_query.end])

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


class DeconvolutedFeatureGraph(FeatureGraph):
    node_cls = DeconvolutedFeatureGraphNode

    def _make_rt_tree(self):
        self.rt_tree = IntervalTreeNode2D.build(
            self.nodes, neutral_mass_point_organizer_callback)

    def find_edges(self, node: DeconvolutedFeatureGraphNode, query_width: float = 2., error_tolerance=1.5e-5, **kwargs):
        query = TimeQuery(node.feature, query_width)
        mass_query = PPMQuery(node.neutral_mass, error_tolerance)
        nodes = self.rt_tree.overlaps(
            [query.start, mass_query.start], [query.end, mass_query.end])
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


class IonMobilityProfileDeconvolutedFeatureGraph(DeconvolutedFeatureGraph):
    node_cls = IonMobilityProfileFeatureGraphNode

    def find_edges(self, node: IonMobilityProfileFeatureGraphNode, query_width: float = 2., error_tolerance=1.5e-5, **kwargs):
        query = TimeQuery(node.feature, query_width)
        mass_query = PPMQuery(node.neutral_mass, error_tolerance)
        nodes = self.rt_tree.overlaps(
            [query.start, mass_query.start], [query.end, mass_query.end])
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
    def smooth(cls, features, time_bridge: float=0.5, mass_error_tolerance: float=1.5e-5):
        return cls(features).connect_components(time_bridge=time_bridge, mass_error_tolerance=mass_error_tolerance)


class GapAwareDeconvolutedFeatureSmoother(GapAwareFeatureSmoother):
    graph: DeconvolutedFeatureGraph

    def __init__(self, features):
        self.graph = DeconvolutedFeatureGraph(features)

    def _wrap_result(self, features) -> 'DeconvolutedLCMSFeatureMap':
        from ms_deisotope.feature_map.feature_map import DeconvolutedLCMSFeatureMap
        return DeconvolutedLCMSFeatureMap(features)


class GapAwareIonMobilityProfileDeconvolutedFeatureSmoother(GapAwareDeconvolutedFeatureSmoother):
    graph: IonMobilityProfileDeconvolutedFeatureGraph

    def __init__(self, features):
        self.graph = IonMobilityProfileDeconvolutedFeatureGraph(features)
