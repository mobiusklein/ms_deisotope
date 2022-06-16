from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin, SimpleInterval

from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature

cdef class FeatureGraphNode(SpanningMixin):
    cdef:
        public LCMSFeature feature
        public size_t index
        public set edges
        public double center
        public double mz
        public object generator # Factor this out somehow. Don't want to add a `__dict__`

    cpdef double _average_center(self)

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


cdef class DeconvolutedFeatureGraphNode(FeatureGraphNode):
    cdef:
        public double neutral_mass
        public int charge


cdef class IonMobilityProfileFeatureGraphNode(DeconvolutedFeatureGraphNode):
    cdef:
        public SimpleInterval ion_mobility_interval


cdef class TimeQuery(SpanningMixin):
    pass


cdef class PPMQuery(SpanningMixin):
    pass


cpdef list connected_components(self)