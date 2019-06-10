cimport cython

from ms_peak_picker._c.peak_set cimport PeakSet, PeakSetIndexed, FittedPeak
from ms_peak_picker._c.peak_index cimport PeakIndex

from brainpy._c.isotopic_distribution cimport TheoreticalPeak

from ms_deisotope._c.scoring cimport IsotopicFitterBase, IsotopicFitRecord
from ms_deisotope._c.averagine cimport AveragineCache, TheoreticalIsotopicPattern

cimport numpy as np


ctypedef fused FittedPeakCollection:
    PeakSet
    PeakIndex


ctypedef fused FittedPeakOrPosition:
    FittedPeak
    size_t
    long


cdef class DeconvoluterBase(object):
    cdef:
        public bint use_subtraction
        public str scale_method
        public bint merge_isobaric_peaks
        public double minimum_intensity

        public PeakSet peaklist

        public IsotopicFitterBase scorer
        public bint verbose
        public dict _slice_cache

    cpdef PeakSet between(self, double m1, double m2)
    cpdef FittedPeak has_peak(self, double mz, double error_tolerance)
    cdef FittedPeak _has_peak(self, double mz, double error_tolerance)

    cpdef list match_theoretical_isotopic_distribution(self, list theoretical_distribution, double error_tolerance=*)
    cpdef scale_theoretical_distribution(self, TheoreticalIsotopicPattern theoretical_distribution, list experimental_distribution)
    cpdef IsotopicFitRecord _evaluate_theoretical_distribution(self, list experimental, TheoreticalIsotopicPattern theoretical, FittedPeak peak, int charge)
    cpdef subtraction(self, TheoreticalIsotopicPattern isotopic_cluster, double error_tolerance=*)
    cpdef list fit_incremental_truncation(self, IsotopicFitRecord seed_fit, double lower_bound)

    cpdef bint _check_fit(self, IsotopicFitRecord fit)

    cpdef list _find_next_putative_peak(self, double mz, int charge, int step=*, double tolerance=*)
    cpdef list _find_previous_putative_peak(self, double mz, int charge, int step=*, double tolerance=*)

    cdef int _find_next_putative_peak_inplace(self, double mz, int charge, set result, int step=*, double tolerance=*)
    cdef int _find_previous_putative_peak_inplace(self, double mz, int charge, set result, int step=*, double tolerance=*)


cdef class AveragineDeconvoluterBase(DeconvoluterBase):
    cdef:
        public AveragineCache averagine

    cpdef IsotopicFitRecord fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge,
                                                         double charge_carrier=*, double truncate_after=*,
                                                         double ignore_below=*)

    cpdef set _fit_peaks_at_charges(self, set peak_charge_set, double error_tolerance, double charge_carrier=*, double truncate_after=*,
                                    double ignore_below=*)


cdef class MultiAveragineDeconvoluterBase(DeconvoluterBase):
    cdef:
        public list averagines

    cpdef IsotopicFitRecord fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge,
                                                         AveragineCache  averagine, double charge_carrier=*, double truncate_after=*,
                                                         double ignore_below=*)
    cpdef set _fit_peaks_at_charges(self, set peak_charge_set, double error_tolerance, double charge_carrier=*, double truncate_after=*,
                                    double ignore_below=*)


cdef bint has_multiple_real_peaks(list peaklist)


cpdef set _get_all_peak_charge_pairs(DeconvoluterBase self, FittedPeak peak, double error_tolerance=*,
                                     object charge_range=*,
                                     int left_search_limit=*, int right_search_limit=*,
                                     bint recalculate_starting_peak=*, bint use_quick_charge=*)

cpdef np.ndarray[np.int32_t, ndim=1] quick_charge(FittedPeakCollection peak_set, FittedPeakOrPosition position_spec, int min_charge, int max_charge)


cdef class ChargeIterator(object):
    cdef:
        public int lower
        public int upper
        public int sign
        int* values
        public size_t size
        public size_t index


    @staticmethod
    cdef ChargeIterator _create(int lo, int hi)

    @staticmethod
    cdef ChargeIterator from_tuple(tuple charge_range)

    @staticmethod
    cdef ChargeIterator from_quickcharge(tuple charge_range, PeakSet peaks, FittedPeak peak)

    cpdef sequence_from_quickcharge(self, PeakSet peaks, FittedPeak peak)
    cpdef make_sequence(self)

    cdef void release_sequence(self)

    cdef void set_bounds(self, int lo, int hi)

    cpdef bint has_more(self)
    cpdef reset(self)
    cpdef int get_next_value(self)

    cdef size_t get_size(self)
