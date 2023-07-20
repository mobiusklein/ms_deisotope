# distutils: language = c++
cimport cython
from cython.operator cimport dereference as deref, preincrement as inc
from cython.parallel cimport prange
from libc.math cimport floor, sqrt
from libc.stdlib cimport malloc, free
from libc.stdio cimport printf

from libcpp.vector cimport vector
from libcpp.unordered_set cimport unordered_set
from libcpp.algorithm cimport sort

from cpython cimport PyErr_SetString

from cpython.object cimport PyObject
from cpython.sequence cimport PySequence_Fast, PySequence_Fast_ITEMS
from cpython.list cimport PyList_Size, PyList_GET_ITEM, PyList_Sort
from cpython.tuple cimport PyTuple_GET_ITEM, PyTuple_Size

from ms_peak_picker._c.peak_index cimport PeakIndex

from ms_peak_picker._c.peak_set cimport PeakBase, FittedPeak, PeakSet
from ms_deisotope._c.peak_set cimport (
    DeconvolutedPeak,
    DeconvolutedPeakSet,
    DeconvolutedPeakSetIndexed,
    deconvoluted_peak_set_t,
    deconvoluted_peak_t,
    _CPeakSet,
    CPeaksFlags,
    create_deconvoluted_peak_set_t,
    free_deconvoluted_peak_set_t,
    deconvoluted_peak_set_all_peaks_for,
    deconvoluted_peak_set_has_peak,
)

cimport numpy as cnp
import numpy as np


cdef struct peak_match_t:
    size_t peak1_index
    size_t peak2_index
    double score


cdef int _sort_score_descending(peak_match_t a, peak_match_t b) nogil:
    return b.score <= a.score


cdef class _AlignableSpectrum:
    cdef:
        public object scan
        public _CPeakSet cpeaks

    def __init__(self, scan):
        self.scan = scan
        if self.scan.deconvoluted_peak_set is None:
            raise ValueError("Must provide a deconvoluted spectrum!")
        self.cpeaks = _CPeakSet.from_peak_list(self.scan.deconvoluted_peak_set)
        if self.cpeaks.ptr == NULL:
            raise MemoryError()
        with nogil:
            self._transform_peak_set(self.cpeaks.ptr)

    def __getitem__(self, i):
        return self.scan.deconvoluted_peak_set[i]

    def __len__(self):
        return len(self.scan.deconvoluted_peak_set)

    @cython.cdivision(True)
    cdef int _transform_peak_set(self, deconvoluted_peak_set_t* peaks) nogil:
        cdef:
            Py_ssize_t i, n
            double total
        n = peaks.size
        total = 0.0
        for i in prange(n):
            total += peaks.peaks[i].intensity
            peaks.peaks[i].intensity = sqrt(peaks.peaks[i].intensity)
        total = sqrt(total)
        if total == 0:
            return 0
        for i in prange(n):
            peaks.peaks[i].intensity /= total
        return 0


cdef enum _PeaksOwnership:
    owning
    borrowing_1
    borrowing_2
    borrowing_1_2


cdef class SpectrumAlignment:
    cdef:
        public object peak_set1
        public object peak_set2
        public double shift
        public double score

        deconvoluted_peak_set_t cpeaks1
        deconvoluted_peak_set_t cpeaks2
        vector[peak_match_t] peak_pairs
        _PeaksOwnership flags


    def __init__(self, peak_set1, peak_set2, double shift, double error_tolerance=2e-5):
        self.peak_set1 = peak_set1
        self.peak_set2 = peak_set2
        self.shift = shift
        self.peak_pairs = vector[peak_match_t]()
        self._init_peak_sets()
        with nogil:
            self.align_peaks(error_tolerance, 0.0)
            if shift != 0:
                self.align_peaks(error_tolerance, shift)
            self.score = self.calculate_score()

    cpdef _init_peak_sets(self):
        self.flags = _PeaksOwnership.owning
        if isinstance(self.peak_set1, DeconvolutedPeakSet):
            create_deconvoluted_peak_set_t(self.peak_set1, &self.cpeaks1)
            with nogil:
                self._transform_peak_set(&self.cpeaks1)
        elif isinstance(self.peak_set1, _AlignableSpectrum):
            self.cpeaks1 = deref((<_AlignableSpectrum>(self.peak_set1)).cpeaks.ptr)
            self.flags = <_PeaksOwnership>(self.flags | _PeaksOwnership.borrowing_1)
        else:
            raise TypeError(f"Cannot align object of type {type(self.peak_set1)}")

        if isinstance(self.peak_set2, DeconvolutedPeakSet):
            create_deconvoluted_peak_set_t(self.peak_set2, &self.cpeaks2)
            with nogil:
                self._transform_peak_set(&self.cpeaks2)
        elif isinstance(self.peak_set2, _AlignableSpectrum):
            self.cpeaks2 = deref((<_AlignableSpectrum>(self.peak_set2)).cpeaks.ptr)
            self.flags = <_PeaksOwnership>(self.flags | _PeaksOwnership.borrowing_2)
        else:
            raise TypeError(f"Cannot align object of type {type(self.peak_set2)}")

    @cython.cdivision(True)
    cdef int _transform_peak_set(self, deconvoluted_peak_set_t* peaks) nogil:
        cdef:
            Py_ssize_t i, n
            double total
        n = peaks.size
        total = 0.0
        for i in prange(n):
            total += peaks.peaks[i].intensity
            peaks.peaks[i].intensity = sqrt(peaks.peaks[i].intensity)
        total = sqrt(total)
        if total == 0:
            return 0
        for i in prange(n):
            peaks.peaks[i].intensity /= total
        return 0

    def __getitem__(self, Py_ssize_t i):
        cdef:
            peak_match_t pair
        if i < self.peak_pairs.size():
            pair = self.peak_pairs[i]
            return self.peak_set1[pair.peak1_index], self.peak_set2[pair.peak2_index], pair.score
        else:
            raise IndexError(i)

    def __len__(self):
        return self.peak_pairs.size()

    cdef int align_peaks(self, double error_tolerance, double shift) nogil:
        cdef:
            size_t i, n
            size_t j, m
            deconvoluted_peak_set_t subset
            peak_match_t pair
            int z1

        pair.score = 0
        n = self.cpeaks1.size
        for i in range(n):
            z1 = self.cpeaks1.peaks[i].charge
            subset = deconvoluted_peak_set_all_peaks_for(
                &self.cpeaks2,
                self.cpeaks1.peaks[i].neutral_mass + shift,
                error_tolerance
            )
            pair.peak1_index = i
            m = subset.size
            for j in range(m):
                if z1 == subset.peaks[j].charge:
                    pair.peak2_index = subset.peaks[j].index
                    self.peak_pairs.push_back(pair)
        return 0

    cdef double calculate_score(self) nogil:
        cdef:
            size_t i, n
            unordered_set[size_t] used_peaks1, used_peaks2
            double total_score
            peak_match_t* pair


        total_score = 0
        n = self.peak_pairs.size()
        for i in range(n):
            pair = &self.peak_pairs[i]
            pair.score = self.cpeaks1.peaks[pair.peak1_index].intensity * self.cpeaks2.peaks[pair.peak2_index].intensity

        sort(self.peak_pairs.begin(), self.peak_pairs.end(), _sort_score_descending)
        for i in range(n):
            pair = &self.peak_pairs[i]
            if used_peaks1.find(pair.peak1_index) == used_peaks1.end() and used_peaks2.find(pair.peak2_index) == used_peaks2.end():
                used_peaks1.insert(pair.peak1_index)
                used_peaks2.insert(pair.peak2_index)
                total_score += pair.score
        return total_score

    def __dealloc__(self):
        if not (self.flags & _PeaksOwnership.borrowing_1):
            free_deconvoluted_peak_set_t(&self.cpeaks1)
        if not (self.flags & _PeaksOwnership.borrowing_2):
            free_deconvoluted_peak_set_t(&self.cpeaks2)