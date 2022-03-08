from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin

from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature


cdef class FeatureGraphNode(SpanningMixin):

    def __init__(self, feature, index, edges=None):
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