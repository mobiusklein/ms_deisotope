from brainpy._c.isotopic_distribution cimport TheoreticalPeak
from ms_peak_picker._c.peak_set cimport PeakSet, FittedPeak
from ms_peak_picker._c.peak_index cimport PeakIndex

cdef class IsotopicFitRecord(object):
    cdef:
        public FittedPeak seed_peak
        public double score
        public int charge
        public list experimental
        public list theoretical
        public FittedPeak monoisotopic_peak
        public int missed_peaks
        public object data

    cpdef bint _eq(self, IsotopicFitRecord other)
    cpdef bint _ne(self, IsotopicFitRecord other)
    cpdef bint _lt(self, IsotopicFitRecord other)
    cpdef bint _gt(self, IsotopicFitRecord other)


cdef class FitSelectorBase(object):
    cdef:
        public double minimum_score

    cpdef IsotopicFitRecord best(self, object results)
    cpdef bint reject(self, IsotopicFitRecord result)
    cpdef bint is_maximizing(self)


cdef class IsotopicFitterBase(object):
    cdef:
        public FitSelectorBase select

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected)
    cpdef bint reject(self, IsotopicFitRecord fit)
    cpdef bint is_maximizing(self)

cdef class LeastSquaresFitter(IsotopicFitterBase):
    pass


cdef class ScaledGTestFitter(IsotopicFitterBase):
    pass

cdef class MSDeconVFitter(IsotopicFitterBase):
    cdef double score_peak(self, FittedPeak obs, TheoreticalPeak theo, double mass_error_tolerance=*, double minimum_signal_to_noise=*) nogil


cdef class PenalizedMSDeconVFitter(IsotopicFitterBase):
    cdef:
        MSDeconVFitter msdeconv
        ScaledGTestFitter penalizer
        double penalty_factor


cdef class FunctionScorer(IsotopicFitterBase):
    cdef:
        object function


cdef class InterferenceDetection(object):
    cdef public PeakIndex peaklist

    cdef double detect_interference(self, list experimental_peaks)


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

    cpdef double _calculate_scale_factor(self, PeakIndex peaklist)
    cpdef void scale_fitted_peaks(self, list experimental, double factor)
    cpdef void scale_theoretical_peaks(self, list theoretical, double factor)
