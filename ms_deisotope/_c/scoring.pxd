from brainpy._c.isotopic_distribution cimport TheoreticalPeak
from ms_peak_picker._c.peak_set cimport PeakSet, FittedPeak

from ms_deisotope._c.deconvoluter_base cimport DeconvoluterBase
from ms_deisotope._c.averagine cimport TheoreticalIsotopicPattern

ctypedef fused fit_collection:
    list
    set
    object


cdef double INFINITY = float('inf')


cdef class IsotopicFitRecord(object):
    cdef:
        public FittedPeak seed_peak
        public double score
        public int charge
        public list experimental
        public TheoreticalIsotopicPattern theoretical
        public FittedPeak monoisotopic_peak
        public int missed_peaks
        public object data
        public Py_hash_t _hash

    @staticmethod
    cdef IsotopicFitRecord _create(FittedPeak seed_peak, double score, int charge, TheoreticalIsotopicPattern theoretical,
                                   list experimental, object data=*, int missed_peaks=*)

    cpdef bint _eq(self, IsotopicFitRecord other)
    cpdef bint _ne(self, IsotopicFitRecord other)
    cpdef bint _lt(self, IsotopicFitRecord other)
    cpdef bint _gt(self, IsotopicFitRecord other)


cdef class FitSelectorBase(object):
    cdef:
        public double minimum_score

    cpdef IsotopicFitRecord best(self, object results)
    cdef IsotopicFitRecord _best_from_set(self, set results)

    cpdef bint reject(self, IsotopicFitRecord result)
    cpdef bint reject_score(self, double score)
    cpdef bint is_maximizing(self)


cdef class IsotopicFitterBase(object):
    cdef:
        public FitSelectorBase select

    cpdef double _evaluate(self, PeakSet peaklist, list observed, list expected)
    cpdef bint reject(self, IsotopicFitRecord fit)
    cpdef bint reject_score(self, double score)
    cpdef bint is_maximizing(self)

    cpdef IsotopicFitterBase _configure(self, DeconvoluterBase deconvoluter, dict kwargs)


cdef class LeastSquaresFitter(IsotopicFitterBase):
    pass


cdef class GTestFitter(IsotopicFitterBase):
    pass

cdef class ScaledGTestFitter(IsotopicFitterBase):
    pass

cdef class MSDeconVFitter(IsotopicFitterBase):
    cdef public double mass_error_tolerance

    cdef double score_peak(self, FittedPeak obs, TheoreticalPeak theo, double mass_error_tolerance=*, double minimum_signal_to_noise=*) nogil


cdef class PenalizedMSDeconVFitter(IsotopicFitterBase):
    cdef:
        public double mass_error_tolerance
        public double penalty_factor


cdef class FunctionScorer(IsotopicFitterBase):
    cdef:
        object function


cdef class InterferenceDetection(object):
    cdef public PeakSet peaklist

    cdef double detect_interference(self, list experimental_peaks, double lower=*, double upper=*)


cdef class DistinctPatternFitter(IsotopicFitterBase):
    cdef:
        public InterferenceDetection interference_detector
        public ScaledGTestFitter g_test_scaled
        public double peak_count_scale
        public double domain_scale

cdef class ScaledPenalizedMSDeconvFitter(IsotopicFitterBase):
    cdef:
        public double scale_factor
        public PenalizedMSDeconVFitter scorer

    cpdef double _calculate_scale_factor(self, PeakSet peaklist)
    cdef void scale_fitted_peaks(self, list experimental, double factor)
    cdef void scale_theoretical_peaks(self, list theoretical, double factor)


cdef class DotProductFitter(IsotopicFitterBase):
    pass
