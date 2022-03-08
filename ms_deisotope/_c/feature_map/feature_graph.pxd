from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin

from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature

cdef class FeatureGraphNode(SpanningMixin):
    cdef:
        public LCMSFeature feature
        public size_t index
        public set edges
        public double center
        public double mz


cdef class FeatureGraphEdge(object):

    cdef:
        public FeatureGraphNode node_a
        public FeatureGraphNode node_b
        public object transition
        public double weight
        public double mass_error
        public double rt_error
        public Py_hash_t _hash

    cpdef remove(self)