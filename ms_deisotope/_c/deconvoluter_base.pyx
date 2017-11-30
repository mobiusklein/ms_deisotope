# cython: embedsignature=True

cimport cython
from libc.stdlib cimport malloc, free

from ms_peak_picker._c.peak_set cimport PeakSet, FittedPeak, PeakSetIndexed
from brainpy._c.isotopic_distribution cimport TheoreticalPeak

from ms_deisotope.constants import ERROR_TOLERANCE as _ERROR_TOLERANCE
from ms_deisotope._c.scoring cimport IsotopicFitterBase, IsotopicFitRecord
from ms_deisotope._c.averagine cimport (AveragineCache, isotopic_shift, PROTON,
                                        TheoreticalIsotopicPattern)

from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE
from cpython.tuple cimport PyTuple_GET_ITEM
from cpython.int cimport PyInt_AsLong, PyInt_Check
from cpython.long cimport PyLong_Check
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.object cimport PyObject

import operator


cdef double ERROR_TOLERANCE = _ERROR_TOLERANCE


cdef size_t count_missed_peaks(list peaklist):
    cdef:
        size_t i
        int n, t
        # FittedPeak peak
        void* peak

    t = n = PyList_GET_SIZE(peaklist)
    for i in range(t):
        # peak = <FittedPeak>PyList_GET_ITEM(peaklist, i)
        peak = PyList_GET_ITEM(peaklist, i)
        if (<FittedPeak>peak).mz > 1 and (<FittedPeak>peak).intensity > 1:
            n -= 1
    return n


cdef FittedPeak make_placeholder_peak(double mz):
    cdef FittedPeak peak = FittedPeak._create(
        mz, intensity=1.0, signal_to_noise=1.0, peak_count=-1, index=0, full_width_at_half_max=0.0,
        area=1.0, left_width=0.0, right_width=0.0)
    return peak


cdef class DeconvoluterBase(object):

    def __init__(self, use_subtraction=False, scale_method="sum", merge_isobaric_peaks=True,
                  minimum_intensity=5., *args, **kwargs):
        self.use_subtraction = use_subtraction
        self.scale_method = scale_method
        self.merge_isobaric_peaks = merge_isobaric_peaks
        self.minimum_intensity = minimum_intensity
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
        return self._has_peak(mz, error_tolerance)

    @cython.final
    cdef FittedPeak _has_peak(self, double mz, double error_tolerance):
        peak = self.peaklist._has_peak(mz, error_tolerance)
        if peak is None or peak.intensity < self.minimum_intensity:
            return make_placeholder_peak(mz)
        return peak

    cpdef list match_theoretical_isotopic_distribution(self, list theoretical_distribution, double error_tolerance=2e-5):
        cdef:
            list experimental_distribution
            size_t i
            double mz

        experimental_distribution = []

        for i in range(PyList_GET_SIZE(theoretical_distribution)):
            mz = (<TheoreticalPeak>PyList_GET_ITEM(theoretical_distribution, i)).mz
            experimental_distribution.append(self._has_peak(mz, error_tolerance))


        return experimental_distribution

    cpdef scale_theoretical_distribution(self, TheoreticalIsotopicPattern theoretical_distribution,
                                         list experimental_distribution):
        cdef:
            size_t i
            TheoreticalPeak peak
            double total_abundance
        return theoretical_distribution.scale(experimental_distribution, self.scale_method)    

    cpdef subtraction(self, TheoreticalIsotopicPattern isotopic_cluster, double error_tolerance=2e-5):
        cdef:
            size_t i
            double existing
            TheoreticalPeak peak
            FittedPeak match
        for i in range(isotopic_cluster.get_size()):
            peak = isotopic_cluster.get(i)
            match = self.peaklist._has_peak(peak.mz, error_tolerance)
            if match is not None:
                existing = match.intensity
                match.intensity -= peak.intensity
                if (match.intensity < 0) or (peak.intensity > (existing * 0.7)):
                    match.intensity = 1.

    def _merge_peaks(self, peak_list):
        peak_list = sorted(peak_list, key=operator.attrgetter("neutral_mass"))
        if not peak_list:
            return []
        current_peak = peak_list[0]
        merged_peaks = []
        for peak in peak_list[1:]:
            if current_peak.neutral_mass == peak.neutral_mass and current_peak.charge == peak.charge:
                current_peak.intensity += peak.intensity
            else:
                merged_peaks.append(current_peak)
                current_peak = peak
        merged_peaks.append(current_peak)
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

        n = peaklist_slice.get_size()
        for i in range(n):
            forward = peaklist_slice.getitem(i)
            prev_peak_mz = forward.mz - (shift * step)
            dummy_peak = make_placeholder_peak(prev_peak_mz)
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
        n = peaklist_slice.get_size()
        for i in range(n):
            backward = peaklist_slice.getitem(i)
            prev_peak_mz = backward.mz
            if step == 1:
                candidates.extend(self._find_next_putative_peak(prev_peak_mz, charge, 1, tolerance))
            else:
                candidates.extend(
                    self._find_previous_putative_peak(prev_peak_mz, charge, step - 1, tolerance))
        return candidates

    def __repr__(self):
        type_name = self.__class__.__name__
        return "%s(peaklist=%s, scorer=%s)" % (type_name, self.peaklist, self.scorer)


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

    def __init__(self, bint use_subtraction=False, str scale_method="sum", bint merge_isobaric_peaks=True,
                 double minimum_intensity=5., *args, **kwargs):
        super(AveragineDeconvoluterBase, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)

    cpdef IsotopicFitRecord fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge,
                                                         double charge_carrier=PROTON, double truncate_after=0.95,
                                                         double ignore_below=0):
        cdef:
            list eid
            TheoreticalIsotopicPattern tid
            double score
        tid = self.averagine.isotopic_cluster(
            peak.mz, charge, charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)
        eid = self.match_theoretical_isotopic_distribution(tid.peaklist, error_tolerance=error_tolerance)
        # self.scale_theoretical_distribution(tid, eid)
        tid._scale(eid, self.scale_method)
        score = self.scorer._evaluate(self.peaklist, eid, tid.peaklist)
        return IsotopicFitRecord._create(peak, score, charge, tid, eid, None, 0)

    cpdef set _fit_peaks_at_charges(self, set peak_charge_set, double error_tolerance, double charge_carrier=PROTON,
                                    double truncate_after=0.95, double ignore_below=0):
        cdef:
            set results
            tuple peak_charge
            IsotopicFitRecord fit
            size_t i

            int charge
            list peak_charge_list
        results = set()
        peak_charge_list = list(peak_charge_set)
        for i in range(PyList_GET_SIZE(peak_charge_list)):
            peak_charge = <tuple>PyList_GET_ITEM(peak_charge_list, i)
            peak = <FittedPeak>PyTuple_GET_ITEM(peak_charge, 0)
            charge = PyInt_AsLong(<object>PyTuple_GET_ITEM(peak_charge, 1))

            if peak.mz < 1:
                continue

            fit = self.fit_theoretical_distribution(
                     peak, error_tolerance, charge,
                     charge_carrier, truncate_after=truncate_after,
                     ignore_below=ignore_below)
            fit.missed_peaks = count_missed_peaks(fit.experimental)
            if not has_multiple_real_peaks(fit.experimental) and fit.charge > 1:
                continue
            if self.scorer.reject(fit):
                continue
            results.add(fit)
        return results


cdef class MultiAveragineDeconvoluterBase(DeconvoluterBase):

    def __init__(self, bint use_subtraction=False, str scale_method="sum", bint merge_isobaric_peaks=True,
                 double minimum_intensity=5., *args, **kwargs):
        super(MultiAveragineDeconvoluterBase, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)

    cpdef IsotopicFitRecord fit_theoretical_distribution(self, FittedPeak peak, double error_tolerance, int charge,
                                                         AveragineCache  averagine, double charge_carrier=PROTON,
                                                         double truncate_after=0.95, double ignore_below=0):
        cdef:
            list eid
            TheoreticalIsotopicPattern tid
            double score
        tid = averagine.isotopic_cluster(
            peak.mz, charge, charge_carrier=charge_carrier,
            truncate_after=truncate_after, ignore_below=ignore_below)
        eid = self.match_theoretical_isotopic_distribution(tid.peaklist, error_tolerance=error_tolerance)
        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer._evaluate(self.peaklist, eid, tid.peaklist)
        return IsotopicFitRecord._create(peak, score, charge, tid, eid, None, 0)

    cpdef set _fit_peaks_at_charges(self, set peak_charge_set, double error_tolerance, double charge_carrier=PROTON,
                                    double truncate_after=0.95, double ignore_below=0):
        cdef:
            list results
            tuple peak_charge
            IsotopicFitRecord fit
            size_t i, j, n_averagine
            int charge
            list peak_charge_list
        results = []
        n_averagine = PyList_GET_SIZE(self.averagines)
        peak_charge_list = list(peak_charge_set)
        for i in range(PyList_GET_SIZE(peak_charge_list)):
            peak_charge = <tuple>PyList_GET_ITEM(peak_charge_list, i)
            peak = <FittedPeak>PyTuple_GET_ITEM(peak_charge, 0)
            charge = PyInt_AsLong(<object>PyTuple_GET_ITEM(peak_charge, 1))

            if peak.mz < 1:
                continue
            for j in range(n_averagine):
                averagine = <AveragineCache>PyList_GET_ITEM(self.averagines, j)
                fit = self.fit_theoretical_distribution(
                    peak, error_tolerance, charge, averagine, charge_carrier,
                    truncate_after=truncate_after, ignore_below=ignore_below)
                fit.missed_peaks = count_missed_peaks(fit.experimental)
                fit.data = averagine
                if not has_multiple_real_peaks(fit.experimental) and fit.charge > 1:
                    continue
                if self.scorer.reject(fit):
                    continue
                # should we track the best fit for each hypothetical peak charge pair
                # and only add the best one to the result set? This would save time
                # later.
                results.append(fit)

        return set(results)


cdef FittedPeak has_previous_peak_at_charge(DeconvoluterBase peak_collection, FittedPeak peak, int charge, int step, double error_tolerance):
    """Get the `step`th *preceding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_collection : DeconvoluterBase
        Peak collection to look up peaks in. Calls :meth:`has_peak` with default accuracy
    peak : ms_peak_picker.FittedPeak
        The peak to use as a point of reference
    charge : int, optional
        The charge state to interpolate from. Defaults to `2`.
    step : int, optional
        The number of peaks along the isotopic pattern to search.

    Returns
    -------
    FittedPeak
    """
    prev = peak.mz - isotopic_shift(charge) * step
    return peak_collection._has_peak(prev, error_tolerance)


cdef FittedPeak has_successor_peak_at_charge(DeconvoluterBase peak_collection, FittedPeak peak, int charge, int step, double error_tolerance):
    """Get the `step`th *succeeding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_collection : DeconvoluterBase
        Peak collection to look up peaks in. Calls :meth:`has_peak` with default accuracy
    peak : ms_peak_picker.FittedPeak
        The peak to use as a point of reference
    charge : int, optional
        The charge state to interpolate from. Defaults to `2`.
    step : int, optional
        The number of peaks along the isotopic pattern to search.

    Returns
    -------
    FittedPeak
    """
    nxt = peak.mz + isotopic_shift(charge) * step
    return peak_collection._has_peak(nxt, error_tolerance)


@cython.final
cdef class ChargeIterator(object):
    cdef:
        public int lower
        public int upper
        public int sign
        int* values
        public size_t size
        public size_t index

    def __init__(self, int lo, int hi):
        self.set_bounds(lo, hi)
        self.make_sequence()

    def __dealloc__(self):
        free(self.values)

    @staticmethod
    cdef ChargeIterator _create(int lo, int hi):
        cdef:
            ChargeIterator inst
        inst = ChargeIterator.__new__(ChargeIterator)
        inst.set_bounds(lo, hi)
        inst.make_sequence()
        return inst

    @staticmethod
    cdef ChargeIterator _from_tuple(tuple charge_range):
        cdef:
            ChargeIterator inst
            int a, b

        a = PyInt_AsLong(<object>PyTuple_GET_ITEM(charge_range, 0))
        b = PyInt_AsLong(<object>PyTuple_GET_ITEM(charge_range, 1))
        return ChargeIterator._create(a, b)

    cdef void set_bounds(self, int lo, int hi):
        cdef:
            int abs_lo, abs_hi
        self.sign = -1 if lo < 0 else 1
        abs_lo, abs_hi = abs(lo), abs(hi)
        if abs_lo < abs_hi:
            self.lower = abs_lo
            self.upper = abs_hi
        else:
            self.lower = abs_hi
            self.upper = abs_lo
        self.size = self.upper

    cdef void make_sequence(self):
        cdef:
            int v
            size_t i, n
        self.index = 0
        n = self.get_size()
        if n == 0:
            return
        self.values = <int*>malloc(sizeof(int) * n)

        for i in range(n):
            self.values[i] = (self.upper - i) * self.sign

    cpdef bint has_more(self):
        return self.index < self.get_size()

    cpdef int get_next_value(self):
        cdef:
            int value
        value = self.values[self.index]
        self.index += 1
        return value

    cdef size_t get_size(self):
        return self.size


@cython.binding(True)
cpdef set _get_all_peak_charge_pairs(DeconvoluterBase self, FittedPeak peak, double error_tolerance=ERROR_TOLERANCE,
                                     object charge_range=(1, 8), int left_search_limit=3, int right_search_limit=3,
                                     bint recalculate_starting_peak=True):
        """Construct the set of all unique candidate (monoisotopic peak, charge state) pairs using
        the provided search parameters.

        The search is performed using :func:`has_previous_peak_at_charge`, :func:`has_successor_peak_at_charge`,
        :meth:`_find_previous_putative_peak`, and :meth:`_find_next_putative_peak`.

        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to ERROR_TOLERANCE
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        recalculate_starting_peak : bool, optional
            Whether or not to re-calculate the putative starting peak m/z based upon nearby
            peaks close to where isotopic peaks for `peak` should be. Defaults to True

        Returns
        -------
        set
            The set of all unique candidate (monoisotopic peak, charge state)
        """
        cdef:
            ChargeIterator charge_iterator
            int charge
            size_t i
            set target_peaks
            FittedPeak prev_peak, nxt_peak
            object add_, update_

        charge_iterator = ChargeIterator._from_tuple(tuple(charge_range))

        target_peaks = set()
        add_ = target_peaks.add
        update_ = target_peaks.update

        while charge_iterator.has_more():
            charge = charge_iterator.get_next_value()
            add_((peak, charge))

            # Look Left
            for i in range(1, left_search_limit):
                prev_peak = has_previous_peak_at_charge(
                    self, peak, charge, i, error_tolerance)
                if prev_peak is None:
                    continue
                add_((prev_peak, charge))

                if recalculate_starting_peak:
                    update_(self._find_previous_putative_peak(
                        peak.mz, charge, i, 2 * error_tolerance))

            # Look Right
            for i in range(1, right_search_limit):
                nxt_peak = has_successor_peak_at_charge(
                    self, peak, charge, i, error_tolerance)
                if nxt_peak is None:
                    continue
                add_((nxt_peak, charge))


                if recalculate_starting_peak:
                    update_(self._find_next_putative_peak(
                        peak.mz, charge, i, 2 * error_tolerance))

            if recalculate_starting_peak:
                for i in range(min(left_search_limit, 2)):
                    update_(self._find_next_putative_peak(
                        peak.mz, charge, step=i, tolerance=2 * error_tolerance))

        return target_peaks
