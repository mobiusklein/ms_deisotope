# cython: profile=True
# cython: embedsignature=True


cimport cython

from ms_peak_picker._c.peak_index cimport PeakIndex
from ms_peak_picker._c.peak_set cimport PeakSet, FittedPeak
from brainpy._c.isotopic_distribution cimport TheoreticalPeak

from ms_deisotope._c.scoring cimport IsotopicFitterBase, IsotopicFitRecord
from ms_deisotope._c.averagine cimport AveragineCache, isotopic_shift

from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE
from cpython.tuple cimport PyTuple_GET_ITEM
from cpython.int cimport PyInt_AsLong, PyInt_Check
from cpython.long cimport PyLong_Check
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.object cimport PyObject

import operator


cdef double sum_intensity(list peaklist):
    cdef:
        size_t i
        double total
        FittedPeak peak
    total = 0
    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <FittedPeak>PyList_GET_ITEM(peaklist, i)
        total += peak.intensity
    return total


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

    def __cinit__(self):
        self.use_subtraction = False
        self.scale_method = 'sum'
        self.merge_isobaric_peaks = True
        self.minimum_intensity = 5.
        self._slice_cache = {}

    cpdef PeakSet between(self, double m1, double m2):
        cdef:
            tuple key
            PyObject* p
            PeakSet region
        key = (m1, m2)
        p = PyDict_GetItem(self._slice_cache, key)
        if p == NULL:
            region = self.peaklist._between(m1, m2)
            PyDict_SetItem(self._slice_cache, key, region)
            return region
        else:
            region = <PeakSet>p
            return region

    cpdef FittedPeak has_peak(self, double mz, double error_tolerance):
        peak = self.peaklist._has_peak(mz, error_tolerance)
        if peak is None or peak.intensity < self.minimum_intensity:
            return FittedPeak(mz, 1.0, 0, 0, 0, 0, 0)
        return peak

    cpdef list match_theoretical_isotopic_distribution(self, list theoretical_distribution, double error_tolerance=2e-5):
        cdef:
            list experimental_distribution
            size_t i
            TheoreticalPeak theo_peak

        experimental_distribution = []

        for i in range(PyList_GET_SIZE(theoretical_distribution)):
            theo_peak = <TheoreticalPeak>PyList_GET_ITEM(theoretical_distribution, i)
            experimental_distribution.append(self.has_peak(theo_peak.mz, error_tolerance))


        return experimental_distribution

    cpdef scale_theoretical_distribution(self, list theoretical_distribution, list experimental_distribution):
        cdef:
            size_t i
            TheoreticalPeak peak
            double total_abundance

        if self.scale_method == 'sum':
            total_abundance = sum_intensity(experimental_distribution)
            for i in range(PyList_GET_SIZE(theoretical_distribution)):
                peak = <TheoreticalPeak>PyList_GET_ITEM(theoretical_distribution, i)            
                peak.intensity *= total_abundance
            return theoretical_distribution
        # elif self.scale_method == 'max':
        #     i, peak = max(enumerate(theoretical_distribution), key=lambda x: x[1].intensity)
        #     scale_factor = experimental_distribution[i].intensity / peak.intensity
        #     for peak in theoretical_distribution:
        #         peak.intensity *= scale_factor
        #     return theoretical_distribution

    cpdef subtraction(self, list isotopic_cluster, double error_tolerance=2e-5):
        cdef:
            size_t i
            TheoreticalPeak peak
            FittedPeak match
        if self.use_subtraction:
            for i in range(PyList_GET_SIZE(isotopic_cluster)):
                peak = <TheoreticalPeak>PyList_GET_ITEM(isotopic_cluster, i)
                match = self.peaklist._has_peak(peak.mz, error_tolerance)
                if match is not None:
                    match.intensity -= peak.intensity
                    if match.intensity < 0:
                        match.intensity = 1.

    def _merge_peaks(self, peak_list):
        peak_list = sorted(peak_list, key=operator.attrgetter("neutral_mass"))
        current_peak = peak_list[0]
        merged_peaks = []
        for peak in peak_list[1:]:
            if current_peak.neutral_mass == peak.neutral_mass and current_peak.charge == peak.charge:
                current_peak.intensity += peak.intensity
            else:
                merged_peaks.append(current_peak)
                current_peak = peak
        return merged_peaks

    cpdef list _find_next_putative_peak(self, double mz, int charge, int step=1, double tolerance=2e-5):
        """
        Recalibrates the current peak location given the position of the next putative peak
        in a theoretical isotopic cluster.

        Suppose that the peak at `mz` is roughly in the neighborhood of a real isotopic peak,
        but the alignment is bad, so it won't make a good starting point for the search for the
        rest of the peaks in its cluster under a stringent error tolerance.

        However, if we're willing to search for the next putative peak with a more permissive error
        tolerance, which we expect will be properly aligned with the rest of its isotopic cluster,
        we can recalibrate the proper starting peak's mz and use that for isotopic cluster fitting.

        Parameters
        ----------
        mz : float
            Starting m/z value to search from
        charge : int
            Charge state to use when calculating the step size in m/z
        step : int, optional
            The number of steps into the putative isotopic cluster to take. Defaults to 1
        tolerance : float, optional
            The error tolerance to accept for finding supporting peaks.

        Returns
        -------
        list
        """
        cdef:
            double shift, next_peak, prev_peak_mz
            PeakSet peaklist_slice
            list candidates
            FittedPeak dummy_peak, forward
            size_t i, n

        shift = isotopic_shift(charge)
        next_peak = mz + (shift * step)
        peaklist_slice = self.between(
            next_peak - (next_peak * tolerance),
            next_peak + (next_peak * tolerance))
        candidates = []

        n = peaklist_slice._get_size()
        for i in range(n):
            forward = peaklist_slice.getitem(i)
            prev_peak_mz = forward.mz - (shift * step)
            dummy_peak = FittedPeak(prev_peak_mz, 1.0, 0, 0, 0, 0, 0)
            candidates.append((dummy_peak, charge))
        return candidates

    cpdef list _find_previous_putative_peak(self, double mz, int charge, int step=1, double tolerance=2e-5):
        cdef:
            double shift, prev_peak, prev_peak_mz
            PeakSet peaklist_slice
            list candidates
            FittedPeak backward
            size_t i, n

        shift = isotopic_shift(charge)
        prev_peak = mz - (shift)
        peaklist_slice = self.between(
            prev_peak - (prev_peak * tolerance),
            prev_peak + (prev_peak * tolerance))
        candidates = []
        n = peaklist_slice._get_size()
        for i in range(n):
            backward = peaklist_slice.getitem(i)
            prev_peak_mz = backward.mz
            if step == 1:
                candidates.extend(self._find_next_putative_peak(prev_peak_mz, charge, 1, tolerance))
            else:
                candidates.extend(
                    self._find_previous_putative_peak(prev_peak_mz, charge, step - 1, tolerance))
        return candidates


cdef bint has_multiple_real_peaks(list peaklist):
    cdef:
        size_t i
        int n
        FittedPeak peak

    n = 0
    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <FittedPeak>PyList_GET_ITEM(peaklist, i)
        if peak.mz > 1 and peak.intensity > 1:
            n += 1
    return n > 1


cdef class AveragineDeconvoluterBase(DeconvoluterBase):
    cdef:
        public AveragineCache averagine

    def __init__(self, *args, **kwargs):
        super(AveragineDeconvoluterBase, self).__init__()

    cpdef IsotopicFitRecord fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge):
        cdef:
            list tid, eid
            double score
        tid = self.averagine.isotopic_cluster(peak.mz, charge)
        eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance=error_tolerance)
        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer._evaluate(self.peaklist, eid, tid)
        return IsotopicFitRecord(peak, score, charge, tid, eid)

    cpdef set _fit_peaks_at_charges(self, set peak_charge_set, double error_tolerance):
        cdef:
            list results
            tuple peak_charge
            IsotopicFitRecord fit
            size_t i
            int charge
            peak_charge_list
        results = []
        peak_charge_list = list(peak_charge_set)
        for i in range(PyList_GET_SIZE(peak_charge_list)):
            peak_charge = <tuple>PyList_GET_ITEM(peak_charge_list, i)
            peak = <FittedPeak>PyTuple_GET_ITEM(peak_charge, 0)
            charge = PyInt_AsLong(<object>PyTuple_GET_ITEM(peak_charge, 1))

            if peak.mz < 1:
                continue

            fit = self.fit_theoretical_distribution(
                     peak, error_tolerance, charge)
            if not has_multiple_real_peaks(fit.experimental) and fit.charge > 1:
                continue
            if self.scorer.reject(fit):
                continue
            results.append(fit)
        return set(results)


cdef class MultiAveragineDeconvoluterBase(DeconvoluterBase):
    cdef:
        public list averagines

    def __init__(self, *args, **kwargs):
        super(MultiAveragineDeconvoluterBase, self).__init__()

    cpdef list fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge):
        cdef:
            list fits
            AveragineCache averagine
            list tid, eid
            double score
        fits = []
        for averagine in self.averagines:
            tid = averagine.isotopic_cluster(peak.mz, charge)
            eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance=error_tolerance)
            self.scale_theoretical_distribution(tid, eid)
            score = self.scorer._evaluate(self.peaklist, eid, tid)
            fits.append(IsotopicFitRecord(score, charge, tid, eid))
        return fits
