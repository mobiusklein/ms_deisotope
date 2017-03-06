from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM


cdef class LCMSFeatureMap(object):
    cdef:
        public list features

    cpdef list _find_all(self, double mz, double error_tolerance)
    cpdef LCMSFeature _search(self, double mz, double error_tolerance)

cpdef tuple binary_search_with_flag(list array, double mz, double error_tolerance)
cdef long binary_search(list array, double mz, double error_tolerance)
cdef void search_sweep(list array, double mz, double error_tolerance, long* loout, long* hiout)

