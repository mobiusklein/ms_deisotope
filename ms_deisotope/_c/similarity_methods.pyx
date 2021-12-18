# cython: profile=True
cimport cython
from cython.parallel cimport prange
from libc.math cimport floor, sqrt
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
from cpython cimport PyErr_SetString

from cpython.object cimport PyObject
from cpython.sequence cimport PySequence_Fast, PySequence_Fast_ITEMS
from cpython.list cimport PyList_Size, PyList_GET_ITEM, PyList_Sort
from cpython.tuple cimport PyTuple_GET_ITEM, PyTuple_Size

from ms_peak_picker._c.peak_index cimport PeakIndex

from ms_peak_picker._c.peak_set cimport PeakBase, FittedPeak, PeakSet
from ms_deisotope._c.peak_set cimport DeconvolutedPeak, DeconvolutedPeakSet, DeconvolutedPeakSetIndexed

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


cpdef enum PeakType:
    Fitted
    Deconvoluted


@cython.freelist(10000)
@cython.final
cdef class PeakMatch(object):
    cdef:
        public PeakBase peak_a
        public PeakBase peak_b
        public double score

    def __init__(self, peak_a, peak_b, score=0.0):
        self.peak_a = peak_a
        self.peak_b = peak_b
        self.score = score

    cpdef bint _lt(self, PeakMatch other):
        return self.score < other.score

    cpdef bint _gt(self, PeakMatch other):
        return self.score > other.score

    cpdef bint _eq(self, PeakMatch other):
        return self.peak_a == other.peak_a and self.peak_b == other.peak_b

    def __richcmp__(self, PeakMatch other, int code):
        if other is None:
            if code == 3:
                return True
            else:
                return False

        if code == 0:
            return self._lt(other)
        elif code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))
        elif code == 4:
            return self._gt(other)

    @staticmethod
    cdef PeakMatch _create(PeakBase peak_a, PeakBase peak_b, double score=0.0):
        cdef PeakMatch self = PeakMatch.__new__(PeakMatch)
        self.peak_a = peak_a
        self.peak_b = peak_b
        self.score = score
        return self

    def __getitem__(self, i):
        if i == 0:
            return self.peak_a
        elif i == 1:
            return self.peak_b
        else:
            raise IndexError(i)

    def __len__(self):
        return 2

    def __repr__(self):
        template = "{self.__class__.__name__}({self.peak_a}, {self.peak_b}, {self.score})"
        return template.format(self=self)


cdef PyObject** _make_fast_array(peak_collection peak_set):
    cdef:
        PyObject** items

    if peak_collection is PeakIndex:
        items = PySequence_Fast_ITEMS(
            PySequence_Fast(peak_set.peaks.peaks, "Error, could not create fast array from PeakIndex"))
    elif peak_collection is PeakSet:
        items = PySequence_Fast_ITEMS(
            PySequence_Fast(peak_set.peaks, "Error, could not create fast array from PeakSet"))
    elif peak_collection is DeconvolutedPeakSet:
        items = PySequence_Fast_ITEMS(
            PySequence_Fast(peak_set.peaks, "Error, could not create fast array from DeconvolutedPeakSet"))
    else:
        # This block is unsafe for any type that is not a list or a tuple, risking segfaults.
        # The list/tuple must have a lifespan exceeding the duration of this function call.
        items = PySequence_Fast_ITEMS(
            PySequence_Fast(peak_set, "Error, could not create fast array from object"))
    return items


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef double peak_set_similarity(peak_collection peak_set_a, peak_collection peak_set_b, int precision=0) except -1:
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
        list _hold_a, _hold_b
        double n_a, n_b, n_ab, z, a, b
        double hi, scaler
        PeakBase peak
        double *bin_a
        double *bin_b
        long i, n, k, index, n_peaks_a, n_peaks_b
        PyObject** items_a
        PyObject** items_b
        PyObject* p_peak

    n_a = 0
    n_b = 0
    a = 0
    b = 0
    z = 0

    # Determine maximum value
    hi = 0
    if peak_collection is PeakIndex:
        peak = peak_set_a.peaks.getitem(peak_set_a.get_size() - 1)
        hi = peak.mz
        peak = peak_set_b.peaks.getitem(peak_set_b.get_size() - 1)
        hi = max(peak.mz, hi)
    elif peak_collection is object:
        hi = 0
        for obj in peak_set_a:
            peak = <PeakBase>obj
            if peak.mz > hi:
                hi = peak.mz
        for obj in peak_set_b:
            peak = <PeakBase>obj
            if peak.mz > hi:
                hi = peak.mz
    elif peak_collection is PeakSet:
        peak = peak_set_a.getitem(peak_set_a.get_size() - 1)
        hi = peak.mz
        peak = peak_set_b.getitem(peak_set_b.get_size() - 1)
        hi = max(peak.mz, hi)
    elif peak_collection is DeconvolutedPeakSet or peak_collection is DeconvolutedPeakSetIndexed:
        peak = peak_set_a._mz_ordered[-1]
        hi = peak.mz
        peak = peak_set_b._mz_ordered[-1]
        hi = max(peak.mz, hi)

    # Calculate bin array size and allocate memory
    scaler = 10 ** precision
    k = <long>((hi * scaler) + 1)
    bin_a = <double*>malloc(k * sizeof(double))
    bin_b = <double*>malloc(k * sizeof(double))
    if bin_a == NULL or bin_b == NULL:
        if bin_a != NULL:
            free(bin_a)
        if bin_b != NULL:
            free(bin_b)
        PyErr_SetString(MemoryError, "Could not allocate bins")
        return -1
    for i in range(k):
        bin_a[i] = bin_b[i] = 0

    # Fill bins
    if peak_collection is object:
        _hold_a = list(peak_set_a)
        n_peaks_a = len(_hold_a)
    else:
        n_peaks_a = peak_set_a.get_size()

    if peak_collection is object:
        _hold_b = list(peak_set_b)
        n_peaks_b = len(_hold_b)
    else:
        n_peaks_b = peak_set_b.get_size()

    if peak_collection is object:
        items_a = _make_fast_array(_hold_a)
        items_b = _make_fast_array(_hold_b)
    else:
        items_a = _make_fast_array(peak_set_a)
        items_b = _make_fast_array(peak_set_b)
    with nogil:
        # Fill bins
        for i in range(n_peaks_a):
            p_peak = items_a[i]
            index = int((<PeakBase>p_peak).mz * scaler)
            bin_a[index] += (<PeakBase>p_peak).intensity

        for i in range(n_peaks_b):
            p_peak = items_b[i]
            index = int((<PeakBase>p_peak).mz * scaler)
            bin_b[index] += (<PeakBase>p_peak).intensity

        # Calculate dot product and normalizers in parallel
        n = k
        for i in prange(n, schedule='static', num_threads=4):
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
        hi = 0
        for obj in peak_set_a:
            peak = <PeakBase>obj
            if peak.mz > hi:
                hi = peak.mz
        for obj in peak_set_b:
            peak = <PeakBase>obj
            if peak.mz > hi:
                hi = peak.mz
    elif peak_collection is DeconvolutedPeakSet or peak_collection is DeconvolutedPeakSetIndexed:
        peak = peak_set_a._mz_ordered[-1]
        hi = peak.mz
        peak = peak_set_b._mz_ordered[-1]
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



cdef class SpectrumAlignment(object):
    cdef:
        public list peak_pairs
        public double shift
        public double score
        public object peak_set_a
        public object peak_set_b
        public PeakType peak_type
        public double norm_aa_cosine
        public double norm_bb_cosine
        public bint normalize
        public bint sqrt_transform

    def __init__(self, peak_set_a, peak_set_b, double error_tolerance=2e-5, double shift=0.0,
                 bint normalize=True, bint sqrt_transform=False):
        self.peak_set_a = peak_set_a
        self.peak_set_b = peak_set_b
        self.peak_pairs = []
        self.shift = shift

        self.score = 0
        self.norm_aa_cosine = 1.0
        self.norm_bb_cosine = 1.0

        self.normalize = normalize
        self.sqrt_transform = sqrt_transform

        self._determine_peak_type()
        self.calculate_normalization()
        self.align(error_tolerance, shift)

    cpdef calculate_normalization(self):
        cdef:
            size_t i, n
            PeakSet peak_set
            DeconvolutedPeakSet dpeak_set
            double acc, intensity, acc_cosine
            bint sqrt_transform

        sqrt_transform = self.sqrt_transform
        if self.peak_type == PeakType.Fitted:
            peak_set = <PeakSet>self.peak_set_a
            n = peak_set.get_size()
            acc_cosine = 0.0
            for i in range(n):
                intensity = peak_set.getitem(i).intensity
                if sqrt_transform:
                    intensity = sqrt(intensity)
                acc_cosine += intensity ** 2
            self.norm_aa_cosine = sqrt(acc_cosine)
            peak_set = <PeakSet>self.peak_set_b
            n = peak_set.get_size()
            acc_cosine = 0.0
            for i in range(n):
                intensity = peak_set.getitem(i).intensity
                if sqrt_transform:
                    intensity = sqrt(intensity)
                acc_cosine += intensity ** 2
            self.norm_bb_cosine = sqrt(acc_cosine)
        elif self.peak_type == PeakType.Deconvoluted:
            dpeak_set = <DeconvolutedPeakSet>self.peak_set_a
            n = dpeak_set.get_size()
            acc_cosine = 0.0
            for i in range(n):
                intensity = dpeak_set.getitem(i).intensity
                if sqrt_transform:
                    intensity = sqrt(intensity)
                acc_cosine += intensity ** 2
            self.norm_aa_cosine = sqrt(acc_cosine)
            dpeak_set = <DeconvolutedPeakSet>self.peak_set_b
            n = dpeak_set.get_size()
            acc_cosine = 0.0
            for i in range(n):
                intensity = dpeak_set.getitem(i).intensity
                if sqrt_transform:
                    intensity = sqrt(intensity)
                acc_cosine += intensity ** 2
            self.norm_bb_cosine = sqrt(acc_cosine)

    cpdef _determine_peak_type(self):
        if isinstance(self.peak_set_a, DeconvolutedPeakSet):
            if isinstance(self.peak_set_b, DeconvolutedPeakSet):
                self.peak_type = PeakType.Deconvoluted
            else:
                raise TypeError(
                    "Peak set types must match, got %s and %s" % (type(self.peak_set_a), type(self.peak_set_b)))
        elif isinstance(self.peak_set_a, PeakSet):
            if isinstance(self.peak_set_b, PeakSet):
                self.peak_type = PeakType.Fitted
            else:
                raise TypeError(
                    "Peak set types must match, got %s and %s" % (type(self.peak_set_a), type(self.peak_set_b)))
        else:
            raise TypeError(
                "Peak set types must match, got %s and %s" % (type(self.peak_set_a), type(self.peak_set_b)))

    cpdef align(self, double error_tolerance=2e-5, double shift=0.0):
        cdef:
            list pairs
        pairs = []
        if self.peak_type == PeakType.Fitted:
            pairs = convolve_peak_sets(<PeakSet>self.peak_set_a, <PeakSet>self.peak_set_b, error_tolerance)
            if shift != 0:
                pairs.extend(
                    convolve_peak_sets(
                        <PeakSet>self.peak_set_a,
                        <PeakSet>self.peak_set_b, error_tolerance, shift))
        elif self.peak_type == PeakType.Deconvoluted:
            pairs = convolve_deconvoluted_peak_sets(
                <DeconvolutedPeakSet>self.peak_set_a,
                <DeconvolutedPeakSet>self.peak_set_b, error_tolerance)
            if shift != 0:
                pairs.extend(
                    convolve_deconvoluted_peak_sets(
                        <DeconvolutedPeakSet>self.peak_set_a,
                        <DeconvolutedPeakSet>self.peak_set_b, error_tolerance, shift))
        self.peak_pairs = pairs
        self._calculate_score()

    cpdef _calculate_score(self):
        self.score = convolved_dot_product(
            self.peak_pairs, self.peak_type, sqrt_transform=self.sqrt_transform)
        if self.normalize:
            self.score /= (self.norm_aa_cosine * self.norm_bb_cosine)



cdef list convolve_peak_sets(PeakSet peak_set_a, PeakSet peak_set_b, double error_tolerance=2e-5, double shift=0.0):
    cdef:
        size_t i, n, j, m
        FittedPeak peak
        FittedPeak other
        tuple peaks_slice
        list peak_pairs
        list pairs_for_peak
        double pt

    peak_pairs = []
    n = peak_set_a.get_size()
    for i in range(n):
        pairs_for_peak = []
        peak = peak_set_a.getitem(i)
        pt = peak.mz + shift
        peaks_slice = peak_set_b.all_peaks_for(pt, error_tolerance)
        m = PyTuple_Size(peaks_slice)
        for j in range(m):
            other = <FittedPeak>PyTuple_GET_ITEM(peaks_slice, j)
            peak_pairs.append(PeakMatch._create(peak, other))
    return peak_pairs


cdef list convolve_deconvoluted_peak_sets(DeconvolutedPeakSet peak_set_a, DeconvolutedPeakSet peak_set_b,
                                          double error_tolerance=2e-5, double shift=0.0):
    cdef:
        size_t i, n, j, m
        DeconvolutedPeak peak
        DeconvolutedPeak other
        tuple peaks_slice
        list peak_pairs
        list pairs_for_peak
        double pt

    peak_pairs = []
    n = peak_set_a.get_size()
    for i in range(n):
        peak = peak_set_a.getitem(i)
        pt = peak.neutral_mass + shift
        peaks_slice = peak_set_b.all_peaks_for(pt, error_tolerance)
        m = PyTuple_Size(peaks_slice)
        for j in range(m):
            other = <DeconvolutedPeak>PyTuple_GET_ITEM(peaks_slice, j)
            peak_pairs.append(PeakMatch._create(peak, other))
    return peak_pairs


cpdef double convolved_dot_product(list peak_pairs, PeakType peak_type=PeakType.Fitted, double norm_a=1.0,
                                   double norm_b=1.0, bint sqrt_transform=False):
    cdef:
        size_t i, n
        PeakBase peak, other
        list pairs
        PeakMatch pair
        double d
        set seen1, seen2
        object key1, key2
    seen1 = set()
    seen2 = set()
    d = 0.0
    n = PyList_Size(peak_pairs)
    for i in range(n):
        pair = <PeakMatch>PyList_GET_ITEM(peak_pairs, i)
        if sqrt_transform:
            pair.score = (sqrt(pair.peak_a.intensity) / norm_a) * (sqrt(pair.peak_b.intensity) / norm_b)
        else:
            pair.score = (pair.peak_a.intensity / norm_a) * (pair.peak_b.intensity / norm_b)
    PyList_Sort(peak_pairs)

    for i in range(n - 1, -1, -1):
        pair = <PeakMatch>PyList_GET_ITEM(peak_pairs, i)
        if peak_type == PeakType.Fitted:
            key1 = (<FittedPeak>pair.peak_a).peak_count
            if key1 in seen1:
                continue
            key2 = (<FittedPeak>pair.peak_b).peak_count
            if key2 in seen2:
                continue
        elif peak_type == PeakType.Deconvoluted:
            key1 = (<DeconvolutedPeak>pair.peak_a)._index.neutral_mass
            if key1 in seen1:
                continue
            key2 = (<DeconvolutedPeak>pair.peak_b)._index.neutral_mass
            if key2 in seen2:
                continue
        seen1.add(key1)
        seen2.add(key2)
        d += pair.score
    return d


cdef double convolve_peak_sets2_fitted(PeakSet peak_set_a, PeakSet peak_set_b, double error_tolerance=2e-5):
    cdef:
        size_t i, n, j, m
        FittedPeak peak
        FittedPeak other
        PeakSet peaks_slice
        double acc

    acc = 0.0
    n = peak_set_a.get_size()
    for i in range(n):
        peak = peak_set_a.getitem(i)
        peaks_slice = peak_set_b._between(
            peak.mz - peak.mz * error_tolerance, peak.mz + peak.mz * error_tolerance)
        m = peaks_slice.get_size()
        for j in range(m):
            other = peaks_slice.getitem(j)
            acc += peak.intensity * other.intensity
    return acc


cdef double convolve_peak_sets2_deconvoluted(DeconvolutedPeakSet peak_set_a, DeconvolutedPeakSet peak_set_b,
                                             double error_tolerance=2e-5):
    cdef:
        size_t i, n, j, m
        DeconvolutedPeak peak
        DeconvolutedPeak other
        tuple peaks_slice
        double acc

    acc = 0.0
    n = peak_set_a.get_size()
    for i in range(n):
        peak = peak_set_a.getitem(i)
        peaks_slice = peak_set_b.all_peaks_for(peak.neutral_mass, error_tolerance)
        m = PyTuple_Size(peaks_slice)
        for j in range(m):
            other = <DeconvolutedPeak>PyTuple_GET_ITEM(peaks_slice, j)
            acc += peak.intensity * other.intensity
    return acc


cpdef double ppm_peak_set_similarity(peak_collection peak_set_a, peak_collection peak_set_b, double error_tolerance=2e-5):
    cdef:
        double ab, aa, bb
    if peak_collection is PeakSet:
        ab = convolve_peak_sets2_fitted(peak_set_a, peak_set_b, error_tolerance)
        aa = convolve_peak_sets2_fitted(peak_set_a, peak_set_a, error_tolerance)
        bb = convolve_peak_sets2_fitted(peak_set_b, peak_set_b, error_tolerance)
    elif peak_collection is PeakIndex:
        ab = convolve_peak_sets2_fitted(peak_set_a.peaks, peak_set_b.peaks, error_tolerance)
        aa = convolve_peak_sets2_fitted(peak_set_a.peaks, peak_set_a.peaks, error_tolerance)
        bb = convolve_peak_sets2_fitted(peak_set_b.peaks, peak_set_b.peaks, error_tolerance)
    elif peak_collection is DeconvolutedPeakSet:
        ab = convolve_peak_sets2_deconvoluted(peak_set_a, peak_set_b, error_tolerance)
        aa = convolve_peak_sets2_deconvoluted(peak_set_a, peak_set_a, error_tolerance)
        bb = convolve_peak_sets2_deconvoluted(peak_set_b, peak_set_b, error_tolerance)
    elif peak_collection is object:
        raise TypeError("Cannot handle objects of type %s" % type(peak_set_a))

    return (ab) / (sqrt((aa)) * sqrt((bb)))


cpdef SpectrumAlignment align_peaks(peak_collection peak_set_a, peak_collection peak_set_b, double error_tolerance=2e-5,
                                    double shift=0.0, bint normalize=True):
    if peak_collection is PeakSet or peak_collection is DeconvolutedPeakSet or peak_collection is DeconvolutedPeakSetIndexed:
        return SpectrumAlignment(peak_set_a, peak_set_b, error_tolerance, shift, normalize, sqrt_transform=False)
    elif peak_collection is PeakIndex:
        return SpectrumAlignment(peak_set_a.peaks, peak_set_b.peaks, error_tolerance, shift, normalize, sqrt_transform=False)
    elif peak_collection is object:
        raise TypeError("Cannot handle objects of type %s" % type(peak_set_a))
