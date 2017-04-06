from ms_deisotope._c.feature_map.lcms_feature cimport FeatureBase, LCMSFeatureTreeNode
from ms_deisotope._c.peak_set cimport DeconvolutedPeak

cimport numpy as np
import numpy as np

cdef class map_coord:
    cdef:
        public double mz
        public double time


cdef class LCMSFeatureSetFit(object):
    cdef:
        public list features
        public list theoretical
        public list supporters
        public double score
        public double mz
        public double neutral_mass
        public int charge
        public object data
        public size_t missing_features
        public size_t n_points
        public FeatureBase monoisotopic_feature
        public np.ndarray scores

    cpdef bint _eq(self, LCMSFeatureSetFit other)
    cpdef bint _ne(self, LCMSFeatureSetFit other)
    cpdef bint _lt(self, LCMSFeatureSetFit other)
    cpdef bint _gt(self, LCMSFeatureSetFit other)


    @staticmethod
    cdef LCMSFeatureSetFit _create(list features, list theoretical, double score,
                                   int charge, size_t missing_features, list supporters,
                                   object data, double neutral_mass, size_t n_points,
                                   np.ndarray scores)
