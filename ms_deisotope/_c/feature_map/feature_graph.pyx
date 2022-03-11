cimport cython

from cpython.list cimport PyList_GetItem, PyList_Size

import numpy as np
cimport numpy as np

from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin, SimpleInterval, IntervalTreeNode2D

from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature, LCMSFeatureTreeNode


cdef class FeatureGraphNode(SpanningMixin):

    def __init__(self, feature, index, edges=None):
        self.start = feature.start_time
        self.end = feature.end_time
        if edges is None:
            edges = set()
        self.feature = feature
        self.index = index
        self.edges = edges
        self.center = self._average_center()
        self.mz = feature.mz

    def as_arrays(self, dtype=float):
        return self.feature.as_arrays(dtype)

    cpdef double _average_center(self):
        cdef:
            double total, abundance, intensity
            LCMSFeatureTreeNode node
            size_t i, n

        total = 0
        abundance = 0
        n = self.feature.get_size()
        for i in range(n):
            node = self.feature.getitem(i)
            intensity = node.total_intensity()
            total += node.time * intensity
            abundance += intensity
        return total / abundance

    @property
    def start_time(self):
        return self.start

    @property
    def end_time(self):
        return self.end

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.feature,)

    def __index__(self):
        return self.index

    def __hash__(self):
        return hash(self.index)

    def __eq__(self, other):
        return self.feature == other.feature


cdef class DeconvolutedFeatureGraphNode(FeatureGraphNode):

    def __init__(self, feature, index, edges=None):
        super().__init__(feature, index, edges)
        self.neutral_mass = feature.neutral_mass
        self.charge = feature.charge


cdef class IonMobilityProfileFeatureGraphNode(DeconvolutedFeatureGraphNode):

    def __init__(self, feature, index, edges=None):
        super().__init__(feature, index, edges)
        self.ion_mobility_interval = feature.ion_mobility_interval


cdef class FeatureGraphEdge(object):
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

    cpdef remove(self):
        self.node_a.edges.remove(self)
        self.node_b.edges.remove(self)



cdef class TimeQuery(SpanningMixin):
    def __init__(self, feature: LCMSFeature, width: float = 0):
        self.start = feature.get_start_time() - width
        self.end = feature.get_end_time() + width

    def __repr__(self):
        return "TimeQuery(%f, %f)" % (self.start, self.end)


cdef class PPMQuery(SpanningMixin):
    def __init__(self, center, error_tolerance):
        self.start = center - center * error_tolerance
        self.end = center + center * error_tolerance

    def __repr__(self):
        return "PPMQuery(%f, %f)" % (self.start, self.end)


cdef (double, double) time_query(LCMSFeature x, double width):
    cdef:
        double out1, out2

    out1 = x.get_start_time() - width
    out2 = x.get_end_time() + width
    return out1, out2


cdef (double, double) ppm_query(double x, double ppm):
    cdef:
        double out1, out2

    out1 = x - x * ppm
    out2 = x + x * ppm
    return out1, out2


@cython.binding(True)
@cython.cdivision(True)
def find_edges(self, FeatureGraphNode node, double query_width=2., double error_tolerance=1.5e-5, **kwargs):
    cdef:
        list nodes
        set edges
        FeatureGraphNode match
        size_t i, n
        IntervalTreeNode2D rt_tree
        double[2] starts, ends
        (double, double) query, mass_query

    query = time_query(node.feature, query_width)
    mass_query = ppm_query(node.mz, error_tolerance)
    starts = [query[0], mass_query[0]]
    ends = [query[1], mass_query[1]]
    rt_tree = self.rt_tree
    nodes = rt_tree.overlaps_2d(starts, ends)

    edge_cls = self.edge_cls
    edges = <set>self.edges

    n = PyList_Size(nodes)
    for i in range(n):
        match = <FeatureGraphNode>PyList_GetItem(nodes, i)
        ppm_error = (node.mz - match.mz) / match.mz
        if abs(ppm_error) > error_tolerance:
            continue
        if match.index == node.index:
            continue
        rt_error = (node.center - match.center)

        edges.add(
            edge_cls(
                node, match, None, mass_error=ppm_error, rt_error=rt_error))


@cython.binding(True)
@cython.cdivision(True)
def find_edges_neutral_mass(self, DeconvolutedFeatureGraphNode node, double query_width=2., double error_tolerance=1.5e-5, **kwargs):
    cdef:
        list nodes
        set edges
        IntervalTreeNode2D rt_tree
        DeconvolutedFeatureGraphNode match
        size_t i, n
        double[2] starts, ends
        (double, double) query, mass_query

    query = time_query(node.feature, query_width)
    mass_query = ppm_query(node.neutral_mass, error_tolerance)
    starts = [query[0], mass_query[0]]
    ends = [query[1], mass_query[1]]
    rt_tree = self.rt_tree
    nodes = rt_tree.overlaps_2d(starts, ends)
        # np.array([query.start, mass_query.start]), np.array([query.end, mass_query.end]))

    edge_cls = self.edge_cls
    edges = <set>self.edges

    n = PyList_Size(nodes)
    for i in range(n):
        match = <DeconvolutedFeatureGraphNode>PyList_GetItem(nodes, i)
        if match.charge != node.charge:
            continue
        ppm_error = (node.neutral_mass - match.neutral_mass) / match.neutral_mass
        if abs(ppm_error) > error_tolerance:
            continue
        if match.index == node.index:
            continue
        rt_error = (node.center - match.center)

        edges.add(
            edge_cls(
                node, match, None, mass_error=ppm_error, rt_error=rt_error))


@cython.binding(True)
@cython.cdivision(True)
def find_edges_neutral_mass_ion_mobility(self, IonMobilityProfileFeatureGraphNode node, double query_width=2., double error_tolerance=1.5e-5, **kwargs):
    cdef:
        list nodes
        set edges
        IntervalTreeNode2D rt_tree
        IonMobilityProfileFeatureGraphNode match
        size_t i, n
        double[2] starts, ends
        (double, double) query, mass_query

    query = time_query(node.feature, query_width)
    mass_query = ppm_query(node.neutral_mass, error_tolerance)
    starts = [query[0], mass_query[0]]
    ends = [query[1], mass_query[1]]
    rt_tree = self.rt_tree
    nodes = rt_tree.overlaps_2d(starts, ends)

    edge_cls = self.edge_cls
    edges = <set>self.edges

    n = PyList_Size(nodes)
    for i in range(n):
        match = <IonMobilityProfileFeatureGraphNode>PyList_GetItem(nodes, i)
        if match.charge != node.charge:
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
        self.edges.add(
            edge_cls(
                node, match, None, mass_error=ppm_error, rt_error=rt_error))


@cython.binding(True)
cpdef list connected_components(self):
    cdef:
        set pool, current_component, visited
        list components
        FeatureGraphNode node
        FeatureGraphEdge edge
        int i, j

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