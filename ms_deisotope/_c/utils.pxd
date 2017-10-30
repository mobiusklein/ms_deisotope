
import numpy as np
cimport numpy as np

from ms_deisotope._c.peak_set cimport DeconvolutedPeakSetIndexed


cpdef list decode_envelopes(np.ndarray[np.float32_t, ndim=1] array)
cpdef DeconvolutedPeakSetIndexed deserialize_deconvoluted_peak_set(dict scan_dict)
