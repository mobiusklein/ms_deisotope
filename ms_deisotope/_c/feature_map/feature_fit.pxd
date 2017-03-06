from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature

cdef class map_coord:
    cdef:
        public double mz
        public double time


cdef class LCMSFeatureSetFit(object):
    cdef:
        public list features
        public list theoretical
        public double score
        public int charge
        public object data
        public size_t missing_features
        public LCMSFeature monoisotopic_feature

    cpdef bint _eq(self, LCMSFeatureSetFit other)
    cpdef bint _ne(self, LCMSFeatureSetFit other)
    cpdef bint _lt(self, LCMSFeatureSetFit other)
    cpdef bint _gt(self, LCMSFeatureSetFit other)