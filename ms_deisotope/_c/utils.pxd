cimport cython

import numpy as np
cimport numpy as np

from ms_peak_picker._c.peak_index cimport PeakIndex
from ms_peak_picker._c.peak_set cimport PeakSet
from ms_deisotope._c.peak_set cimport DeconvolutedPeakSet, DeconvolutedPeakSetIndexed, PeakBase


cdef double ppm_error(double x, double y)
cdef double ppm2da(double mass, double error_tolerance)
cdef double da2ppm(double mass, double error_tolerance)


cpdef intensity_getter(peak)
cpdef mz_getter(peak)
cpdef mass_getter(peak)


cpdef list decode_envelopes(np.ndarray[np.float32_t, ndim=1] array)
cpdef DeconvolutedPeakSetIndexed deserialize_deconvoluted_peak_set(dict scan_dict, bint include_envelopes=*)
cpdef DeconvolutedPeakSetIndexed build_deconvoluted_peak_set_from_arrays(np.ndarray[double, ndim=1] mz_array,
                                                                         np.ndarray[double, ndim=1] intensity_array,
                                                                         np.ndarray[long, ndim=1] charge_array)


ctypedef fused peak_collection:
    PeakSet
    PeakIndex
    DeconvolutedPeakSet
    DeconvolutedPeakSetIndexed
    object


cpdef double _peak_sequence_tic(self, peak_collection peaks) except -1
cpdef PeakBase _peak_sequence_bp(self, peak_collection peaks)

cpdef double correlation(cython.floating[:] x, cython.floating[:] y)