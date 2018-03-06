
import numpy as np
cimport numpy as np

from ms_deisotope._c.peak_set cimport DeconvolutedPeakSetIndexed


cpdef list decode_envelopes(np.ndarray[np.float32_t, ndim=1] array)
cpdef DeconvolutedPeakSetIndexed deserialize_deconvoluted_peak_set(dict scan_dict)
cpdef DeconvolutedPeakSetIndexed build_deconvoluted_peak_set_from_arrays(np.ndarray[double, ndim=1] mz_array,
                                                                         np.ndarray[double, ndim=1] intensity_array,
                                                                         np.ndarray[long, ndim=1] charge_array)
