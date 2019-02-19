cimport cython
from cython.parallel cimport prange
from libc.math cimport floor, sqrt
from libc.stdlib cimport calloc, free
from cpython.tuple cimport PyTuple_GET_ITEM
from ms_peak_picker._c.peak_index cimport PeakIndex
from ms_peak_picker._c.peak_set cimport PeakBase, FittedPeak, PeakSet
from ms_deisotope._c.peak_set cimport DeconvolutedPeak, DeconvolutedPeakSet

cimport numpy as cnp
import numpy as np

cdef object zeros = np.zeros


ctypedef fused peak_collection:
    PeakSet
    PeakIndex
    DeconvolutedPeakSet
    object


cpdef enum SimilarityMetrics:
    dot_product
    normalized_dot_product


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef double peak_set_similarity(peak_collection peak_set_a, peak_collection peak_set_b, int precision=0):
    """Computes the normalized dot product, also called cosine similarity between
    two peak sets, a similarity metric ranging between 0 (dissimilar) to 1.0 (similar).

    Parameters
    ----------
    peak_set_a : Iterable of Peak-like
    peak_set_b : Iterable of Peak-like
        The two peak collections to compare. It is usually only useful to
        compare the similarity of peaks of the same class, so the types
        of the elements of `peak_set_a` and `peak_set_b` should match.
    precision : int, optional
        The precision of rounding to use when binning spectra. Defaults to 0

    Returns
    -------
    float
        The similarity between peak_set_a and peak_set_b. Between 0.0 and 1.0
    """

    cdef:
        tuple binned_spectra
        double n_a, n_b, n_ab, z, a, b
        double hi, scaler
        PeakBase peak
        double *bin_a
        double *bin_b
        long i, n, k, index

    hi = 0
    if peak_collection is PeakIndex:
        peak = peak_set_a.peaks.getitem(peak_set_a.get_size() - 1)
        hi = peak.mz
        peak = peak_set_b.peaks.getitem(peak_set_b.get_size() - 1)
        hi = max(peak.mz, hi)
    elif peak_collection is object:
        peak = peak_set_a[-1]
        hi = peak.mz
        peak = peak_set_b[-1]
        hi = max(peak.mz, hi)
    else:
        peak = peak_set_a.getitem(peak_set_a.get_size() - 1)
        hi = peak.mz
        peak = peak_set_b.getitem(peak_set_b.get_size() - 1)
        hi = max(peak.mz, hi)
    scaler = 10 ** precision
    k = <long>((hi * scaler) + 1)
    bin_a = <double*>calloc(sizeof(double), k)
    bin_b = <double*>calloc(sizeof(double), k)

    if peak_collection is object:
        n = len(peak_set_a)
    else:
        n = peak_set_a.get_size()
    for i in range(n):
        if peak_collection is PeakIndex:
            peak = peak_set_a.peaks.getitem(i)
        elif peak_collection is object:
            peak = peak_set_a[i]
        else:
            peak = peak_set_a.getitem(i)
        index = int(peak.mz * scaler)
        bin_a[index] += peak.intensity

    if peak_collection is object:
        n = len(peak_set_b)
    else:
        n = peak_set_b.get_size()
    for i in range(n):
        if peak_collection is PeakIndex:
            peak = peak_set_b.peaks.getitem(i)
        elif peak_collection is object:
            peak = peak_set_b[i]
        else:
            peak = peak_set_b.getitem(i)
        index = int(peak.mz * scaler)
        bin_b[index] += peak.intensity

    n_a = 0
    n_b = 0
    z = 0
    n = k
    for i in prange(n, nogil=True, schedule='static', num_threads=4):
        a = bin_a[i]
        b = bin_b[i]
        n_a += a * a
        n_b += b * b
        z += a * b
    free(bin_a)
    free(bin_b)
    n_ab = sqrt(n_a) * sqrt(n_b)
    if n_ab == 0:
        return 0.0
    return z / n_ab


@cython.boundscheck(False)
cpdef tuple bin_peaks(peak_collection peak_set_a, peak_collection peak_set_b, int precision=0):
    cdef:
        double hi, scaler
        PeakBase peak
        cnp.ndarray[double, ndim=1] bin_a, bin_b
        size_t i, n, index

    if peak_collection is PeakIndex:
        hi = 0
        peak = peak_set_a.peaks.getitem(peak_set_a.get_size() - 1)
        hi = peak.mz
        peak = peak_set_b.peaks.getitem(peak_set_b.get_size() - 1)
        hi = max(peak.mz, hi)
    elif peak_collection is object:
        peak = peak_set_a[-1]
        hi = peak.mz
        peak = peak_set_b[-1]
        hi = max(peak.mz, hi)
    else:
        hi = 0
        peak = peak_set_a.getitem(peak_set_a.get_size() - 1)
        hi = peak.mz
        peak = peak_set_b.getitem(peak_set_b.get_size() - 1)
        hi = max(peak.mz, hi)
    scaler = 10 ** precision
    hi = (hi * scaler) + 1
    bin_a = zeros(<int>hi)
    bin_b = zeros(<int>hi)
    
    if peak_collection is object:
        n = len(peak_set_a)
    else:
        n = peak_set_a.get_size()
    for i in range(n):
        if peak_collection is PeakIndex:
            peak = peak_set_a.peaks.getitem(i)
        elif peak_collection is object:
            peak = peak_set_a[i]
        else:
            peak = peak_set_a.getitem(i)
        index = int(peak.mz * scaler)
        bin_a[index] += peak.intensity

    if peak_collection is object:
        n = len(peak_set_b)
    else:
        n = peak_set_b.get_size()
    for i in range(n):
        if peak_collection is PeakIndex:
            peak = peak_set_b.peaks.getitem(i)
        elif peak_collection is object:
            peak = peak_set_b[i]
        else:
            peak = peak_set_b.getitem(i)
        index = int(peak.mz * scaler)
        bin_b[index] += peak.intensity
    return (bin_a, bin_b)


@cython.boundscheck(False)
@cython.cdivision(True)
cdef double calculate_dot_product(cnp.ndarray[double, ndim=1] bin_a, cnp.ndarray[double, ndim=1] bin_b, bint normalize=True):
    cdef:
        double n_a, n_b, n_ab, z, a, b
        long i, n
    with nogil:
        n_a = 0
        n_b = 0
        z = 0
        n = bin_a.shape[0]
        for i in prange(n, schedule='static'):
            a = bin_a[i]
            b = bin_b[i]
            n_a += a * a
            n_b += b * b
            z += a * b
        if normalize:
            n_ab = sqrt(n_a) * sqrt(n_b)
            if n_ab == 0:
                return 0.0
            return z / n_ab
        return z

