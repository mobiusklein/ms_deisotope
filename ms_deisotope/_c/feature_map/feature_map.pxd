from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM


cdef class LCMSFeatureMap(object):
    cdef:
        public list features

    cdef Py_ssize_t get_size(self)
    cdef LCMSFeature get(self, size_t i)
    cpdef list find_all(self, double mz, double error_tolerance=*)
    cpdef LCMSFeature search(self, double mz, double error_tolerance=*)
    cpdef list between(self, double lo, double hi, double error_tolerance=*)
    cpdef list spanning_time(self, double time_point)
    cpdef LCMSFeatureMap clone(self, bint deep=*)

cpdef tuple binary_search_with_flag(list array, double mz, double error_tolerance)
cdef long binary_search(list array, double mz, double error_tolerance)
cdef void search_sweep(list array, double mz, double error_tolerance, long* loout, long* hiout)

cpdef list split_sparse(LCMSFeatureMap self, double delta_rt=*, size_t min_size=*)
