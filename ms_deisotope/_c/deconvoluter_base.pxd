from ms_peak_picker._c.peak_index cimport PeakIndex
from ms_peak_picker._c.peak_set cimport PeakSet, FittedPeak
from brainpy._c.isotopic_distribution cimport TheoreticalPeak

from ms_deisotope._c.scoring cimport IsotopicFitterBase, IsotopicFitRecord
from ms_deisotope._c.averagine cimport AveragineCache


cdef class DeconvoluterBase(object):
    cdef:
        public bint use_subtraction
        public str scale_method
        public bint merge_isobaric_peaks
        public double minimum_intensity

        public PeakIndex peaklist

        public IsotopicFitterBase scorer
        public bint verbose
        public dict _slice_cache

    cpdef PeakSet between(self, double m1, double m2)
    cpdef FittedPeak has_peak(self, double mz, double error_tolerance)
    cpdef list match_theoretical_isotopic_distribution(self, list theoretical_distribution, double error_tolerance=*)
    cpdef scale_theoretical_distribution(self, list theoretical_distribution, list experimental_distribution)
    cpdef subtraction(self, list isotopic_cluster, double error_tolerance=*)
    cpdef list _find_next_putative_peak(self, double mz, int charge, int step=*, double tolerance=*)
    cpdef list _find_previous_putative_peak(self, double mz, int charge, int step=*, double tolerance=*)


cdef class AveragineDeconvoluterBase(DeconvoluterBase):
    cdef:
        public AveragineCache averagine

    cpdef IsotopicFitRecord fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge,
                                                         double charge_carrier=*)

    cpdef set _fit_peaks_at_charges(self, set peak_charge_set, double error_tolerance, double charge_carrier=*)


cdef class MultiAveragineDeconvoluterBase(DeconvoluterBase):
    cdef:
        public list averagines

    cpdef IsotopicFitRecord fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge,
                                                         AveragineCache  averagine, double charge_carrier=*)
    cpdef set _fit_peaks_at_charges(self, set peak_charge_set, double error_tolerance, double charge_carrier=*)
