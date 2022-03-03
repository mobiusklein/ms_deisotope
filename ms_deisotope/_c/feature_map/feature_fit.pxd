from ms_deisotope._c.feature_map.lcms_feature cimport FeatureBase, LCMSFeatureTreeNode, LCMSFeature
from ms_deisotope._c.peak_set cimport DeconvolutedPeak, Envelope
from ms_deisotope._c.averagine cimport TheoreticalIsotopicPattern

cimport numpy as np
import numpy as np

cdef class map_coord:
    cdef:
        public double mz
        public double time

    cdef bint ge(self, map_coord other)
    cdef bint le(self, map_coord other)

    @staticmethod
    cdef map_coord _create(double mz, double time)

    cpdef map_coord copy(self)


cdef class LCMSFeatureSetFit(object):
    cdef:
        public list features
        public TheoreticalIsotopicPattern theoretical
        public list supporters
        public double score
        public double mz
        public double neutral_mass
        public int charge
        public object data
        public size_t missing_features
        public size_t n_points
        public FeatureBase monoisotopic_feature
        public object scores
        public object times

    cpdef bint _eq(self, LCMSFeatureSetFit other)
    cpdef bint _ne(self, LCMSFeatureSetFit other)
    cpdef bint _lt(self, LCMSFeatureSetFit other)
    cpdef bint _gt(self, LCMSFeatureSetFit other)

    cpdef int count_null_features(self)
    cpdef bint has_multiple_real_features(self)

    cdef map_coord get_start(self)
    cdef map_coord get_end(self)

    @staticmethod
    cdef LCMSFeatureSetFit _create(list features, TheoreticalIsotopicPattern theoretical,
                                   double score, int charge, size_t missing_features,
                                   list supporters, object data, double neutral_mass,
                                   size_t n_points, object scores, object times)


cpdef Envelope _sum_envelopes(LCMSFeature self)