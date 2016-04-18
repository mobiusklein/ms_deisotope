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

    cpdef bint _eq(self, IsotopicFitRecord other)
    cpdef bint _ne(self, IsotopicFitRecord other)
    cpdef bint _lt(self, IsotopicFitRecord other)
    cpdef bint _gt(self, IsotopicFitRecord other)


cdef class FitSelectorBase(object):
    cdef:
        public double minimum_score

    cpdef IsotopicFitRecord best(self, object results)
    cpdef bint reject(self, IsotopicFitRecord result)


cdef class IsotopicFitterBase(object):
    cdef:
        public FitSelectorBase select

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected)
    cpdef bint reject(self, IsotopicFitRecord fit)


cdef class LeastSquaresFitter(IsotopicFitterBase):
    pass


cdef class ScaledGTestFitter(IsotopicFitterBase):
    pass

cdef class MSDeconVFitter(IsotopicFitterBase):

    cdef double _reweight(self, FittedPeak obs, TheoreticalPeak theo, double scale_obs, double scale_theo) nogil
    cpdef double reweight(self, FittedPeak obs, TheoreticalPeak theo, double scale_obs, double scale_theo)
    cdef double score_peak(self, FittedPeak obs, TheoreticalPeak theo, double mass_error_tolerance=*, double minimum_signal_to_noise=*) nogil


cdef class PenalizedMSDeconVFitter(IsotopicFitterBase):
    cdef:
        MSDeconVFitter msdeconv
        ScaledGTestFitter penalizer


cdef class FunctionScorer(IsotopicFitterBase):
    cdef:
        object function
