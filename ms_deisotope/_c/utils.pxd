
import numpy as np
cimport numpy as np

from ms_deisotope._c.peak_set cimport DeconvolutedPeakSet


cpdef list decode_envelopes(np.ndarray[np.float32_t, ndim=1] array)
cpdef DeconvolutedPeakSet marshal_deconvoluted_peak_set(dict scan_dict)
