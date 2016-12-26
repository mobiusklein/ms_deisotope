# -*- coding: utf-8 -*-
import operator

from .averagine import (
    AveragineCache, peptide, glycopeptide, glycan, neutral_mass, isotopic_variants,
    isotopic_shift, PROTON, shift_isotopic_pattern)
from .peak_set import DeconvolutedPeak, DeconvolutedPeakSolution, DeconvolutedPeakSet
from .scoring import IsotopicFitRecord, penalized_msdeconv, distinct_pattern_fitter
from .utils import range, Base, TrivialTargetedDeconvolutionResult, DeconvolutionProcessResult
from .envelope_statistics import a_to_a2_ratio, average_mz, most_abundant_mz
from .peak_dependency_network import PeakDependenceGraph, NetworkedTargetedDeconvolutionResult

from ms_peak_picker import FittedPeak

import logging

logger = logging.getLogger("deconvolution")
info = logger.info
debug = logger.debug


def mean(numbers):
    n = 0.
    total = 0
    for x in numbers:
        n += 1.
        total += x
    return total / n


def has_previous_peak_at_charge(peak_index, peak, charge=2, step=1, error_tolerance=2e-5):
    """Get the `step`th *preceding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_index : ms_peak_picker.PeakIndex
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
    return peak_index.has_peak(prev, error_tolerance)


def has_successor_peak_at_charge(peak_index, peak, charge=2, step=1, error_tolerance=2e-5):
    """Get the `step`th *succeeding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_index : ms_peak_picker.PeakIndex
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
    return peak_index.has_peak(nxt, error_tolerance)


def first_peak(peaks):
    """Get the first non-placeholder peak in a list of peaks

    Parameters
    ----------
    peaks : Iterable of FittedPeak

    Returns
    -------
    FittedPeak
    """
    for peak in peaks:
        if peak.intensity > 1 and peak.mz > 1:
            return peak


def drop_placeholders(peaks):
    """Removes all placeholder peaks from an iterable of peaks

    Parameters
    ----------
    peaks : Iterable of FittedPeak

    Returns
    -------
    list
    """
    return [peak for peak in peaks if peak.mz > 1 and peak.intensity > 1]


def count_placeholders(peaks):
    """Counts the number of placeholder peaks in an iterable
    of FittedPeaks

    Parameters
    ----------
    peaks : Iterable of FittedPeak

    Returns
    -------
    int
        Number of placeholder peaks
    """
    i = 0
    for peak in peaks:
        if peak.intensity <= 1:
            i += 1
    return i


def drop_placeholders_parallel(peaks, otherpeaks):
    """Given two parallel iterables of Peak objects, `peaks` and `otherpeaks`,
    for each position that is not a placeholder in `peaks`, include that Peak object
    and its counterpart in `otherpeaks` in a pair of output lists.

    Parameters
    ----------
    peaks : Iterable of FittedPeak
        Peak collection to filter against
    otherpeaks : Iterable
        Collection of objects (Peak-like) to include based upon
        contents of `peaks`

    Returns
    -------
    list
        Filtered form of `peaks`
    list
        Filtered form of `otherpeaks`
    """
    new_peaks = []
    new_otherpeaks = []
    for i in range(len(peaks)):
        peak = peaks[i]
        if peak.intensity > 1:
            new_peaks.append(peak)
            new_otherpeaks.append(otherpeaks[i])
    return new_peaks, new_otherpeaks


class DeconvoluterBase(Base):
    """Base class for all Deconvoluter types. Provides basic configuration for common operations,
    regardless of implementation. Because these methods form the backbone of all deconvolution algorithms,
    this class has a C-extension implementation as well.

    Attributes
    ----------
    merge_isobaric_peaks : bool
        If multiple passes produce peaks with identical mass values,
        should those peaks be summed
    minimum_intensity : float
        Experimental peaks whose intensity is below this level will be ignored
        by peak querying methods
    scale_method : str
        The name of the method to use to scale theoretical isotopic pattern intensities
        to match the experimental isotopic pattern
    use_subtraction : bool
        Whether or not to apply a subtraction procedure to experimental peaks after they
        have been fitted. This is only necessary if the same signal may be examined multiple
        times as in a multi-pass method or when peak dependence is not considered
    verbose : bool
        Produce extra logging information
    """
    use_subtraction = False
    scale_method = 'sum'
    merge_isobaric_peaks = True
    minimum_intensity = 5.
    verbose = False

    def __init__(self, use_subtraction=False, scale_method="sum", merge_isobaric_peaks=True,
                 minimum_intensity=5., *args, **kwargs):
        self.use_subtraction = use_subtraction
        self.scale_method = scale_method
        self.merge_isobaric_peaks = merge_isobaric_peaks
        self.minimum_intensity = minimum_intensity
        self._slice_cache = {}

    def has_peak(self, mz, error_tolerance):
        """Query :attr:`peaklist` for a peak at `mz` within `error_tolerance` ppm. If a peak
        is not found, this method returns a placeholder peak.

        Parameters
        ----------
        mz : float
            The m/z to search for a peak at
        error_tolerance : float
            The parts-per-million error tolerance to search with

        Returns
        -------
        FittedPeak
            A peak from :attr:`peaklist` if present, else a placeholder peak.
        """
        peak = self.peaklist.has_peak(mz, error_tolerance)
        if peak is None or peak.intensity < self.minimum_intensity:
            return FittedPeak(mz, 1.0, 0, 0, 0, 0, 0)
        return peak

    def between(self, m1, m2):
        """Take a slice from :attr:`peaklist` using it's :meth:`between` method. Caches
        repeated queries in :attr:`_scan_cache`.

        Parameters
        ----------
        m1 : float
            The low m/z to slice from
        m2 : float
            The high m/z to slice from

        Returns
        -------
        PeakSet
            The subset of peaks from :attr:`peaklist` between `m1` and `m2`
        """
        key = (m1, m2)
        if key in self._slice_cache:
            return self._slice_cache[key]
        else:
            region = self.peaklist.between(m1, m2)
            self._slice_cache[key] = region
            return region

    def match_theoretical_isotopic_distribution(self, theoretical_distribution, error_tolerance=2e-5):
        """Given a list of theoretical peaks, find their counterparts in :attr:`peaklist` within `error_tolerance`
        ppm error. If no experimental peak is found, a placeholder will be used in its stead.

        Parameters
        ----------
        theoretical_distribution : list of TheoreticalPeak
            The theoretical isotopic pattern to match
        error_tolerance : float, optional
            Parts-per-million error tolerance to permit in searching for matches

        Returns
        -------
        list of FittedPeak
            The list of matched peaks
        """
        experimental_distribution = [self.has_peak(
            p.mz, error_tolerance) for p in theoretical_distribution]
        return experimental_distribution

    def scale_theoretical_distribution(self, theoretical_distribution, experimental_distribution):
        """Scale up a theoretical isotopic pattern such that its total intensity matches the experimental
        isotopic pattern. This mutates `theoretical_distribution`.

        This assumes that the sum of intensities in `theoretical_distribution` is 1.0. This method is also
        not particularly kind to isotopic patterns where there are missing peaks.

        The scaling algorithm used is controlled by :attr:`scale_method`

        Parameters
        ----------
        theoretical_distribution : list of TheoreticalPeak
            The theoretical isotopic pattern to scale up
        experimental_distribution : list of FittedPeak
            The experimental isotopic pattern to use as a reference

        Returns
        -------
        list of TheoreticalPeak
        """
        if self.scale_method == 'sum':
            total_abundance = sum(
                p.intensity for p in experimental_distribution)
            for peak in theoretical_distribution:
                peak.intensity *= total_abundance
            return theoretical_distribution
        elif self.scale_method == 'max':
            i, peak = max(enumerate(theoretical_distribution),
                          key=lambda x: x[1].intensity)
            scale_factor = experimental_distribution[
                i].intensity / peak.intensity
            for peak in theoretical_distribution:
                peak.intensity *= scale_factor
            return theoretical_distribution

    def subtraction(self, isotopic_cluster, error_tolerance=2e-5):
        """Subtract signal attributed to `isotopic_cluster` from the equivalent
        peaks in :attr:`peaklist`, mutating the peaks within.

        This will change the intensity value of the matched `FittedPeak` instances,
        and this is reflected in all of there references.

        Parameters
        ----------
        isotopic_cluster : list of TheoreticalPeak
            The isotopic cluster to subtract
        error_tolerance : float, optional
            Parts-per-million mass accuracy error tolerance to permit when
            finding matches for `isotopic_cluster`
        """
        for peak in isotopic_cluster:
            match = self.peaklist.has_peak(peak.mz, error_tolerance)
            if match is not None:
                match.intensity -= peak.intensity
                if match.intensity < 0:
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

    def _find_next_putative_peak(self, mz, charge, step=1, tolerance=2e-5):
        """
        Recalibrates the current peak location given the position of the **next** putative peak
        in a theoretical isotopic cluster.

        Suppose that the peak at `mz` is roughly in the neighborhood of a real isotopic peak,
        but the alignment is bad, so it won't make a good starting point for the search for the
        rest of the peaks in its cluster under a stringent error tolerance.

        However, if we're willing to search for the **next** putative peak with a more permissive error
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
        shift = isotopic_shift(charge)
        next_peak = mz + (shift * step)
        peaklist_slice = self.between(
            next_peak - (next_peak * tolerance),
            next_peak + (next_peak * tolerance))
        candidates = []
        for forward in peaklist_slice:
            prev_peak_mz = forward.mz - (shift * step)
            dummy_peak = FittedPeak(prev_peak_mz, 1.0, 0, 0, 0, 0, 0)
            candidates.append((dummy_peak, charge))
        return candidates

    def _find_previous_putative_peak(self, mz, charge, step=1, tolerance=2e-5):
        """
        Recalibrates the current peak location given the position of the **previous** putative peak
        in a theoretical isotopic cluster.

        Suppose that the peak at `mz` is roughly in the neighborhood of a real isotopic peak,
        but the alignment is bad, so it won't make a good starting point for the search for the
        rest of the peaks in its cluster under a stringent error tolerance.

        However, if we're willing to search for the **previous** putative peak with a more permissive error
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
        shift = isotopic_shift(charge)
        prev_peak = mz - (shift)
        peaklist_slice = self.between(
            prev_peak - (prev_peak * tolerance),
            prev_peak + (prev_peak * tolerance))
        candidates = []
        for backward in peaklist_slice:
            prev_peak_mz = backward.mz
            if step == 1:
                candidates.extend(self._find_next_putative_peak(
                    prev_peak_mz, charge, 1, tolerance))
            else:
                candidates.extend(
                    self._find_previous_putative_peak(prev_peak_mz, charge, step - 1, tolerance))
        return candidates

    def __repr__(self):
        type_name = self.__class__.__name__
        return "%s(peaklist=%s)" % (type_name, self.peaklist)


class AveragineDeconvoluterBase(DeconvoluterBase):
    """A base class derived from :class:`DeconvoluterBase` which provides some common methods
    for fitting isotopic patterns using an Averagine model.

    Because these methods form the backbone of all deconvolution algorithms, this class has a
    C-extension implementation as well.
    """
    def __init__(self, use_subtraction=False, scale_method="sum", merge_isobaric_peaks=True,
                 minimum_intensity=5., *args, **kwargs):
        super(AveragineDeconvoluterBase, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)

    def fit_theoretical_distribution(self, peak, error_tolerance, charge, charge_carrier=PROTON, truncate_after=0.8):
        """Fit an isotopic pattern seeded at `peak` at `charge` charge.

        Generates a theoretical isotopic pattern using :attr:`averagine`, calls
        :meth:`match_theoretical_isotopic_distribution`
        to extract experimental peaks matching this theoretical pattern, scales the theoretical distribution using
        :meth:`scale_theoretical_distribution`, and evaluates the quality of the fit using :attr:`scorer`.

        Parameters
        ----------
        peak : FittedPeak
            The putative monoisotopic peak to use for interpolating an isotopic pattern
        error_tolerance : float
            Parts-per-million error tolerance for isotopic pattern matching
        charge : int
            The charge state to produce an isotopic pattern for
        charge_carrier : float, optional
            The charge carrier mass, defaults to `PROTON`

        Returns
        -------
        IsotopicFitRecord
            The fitted isotopic pattern
        """
        tid = self.averagine.isotopic_cluster(
            peak.mz, charge, charge_carrier=charge_carrier,
            truncate_after=truncate_after)
        eid = self.match_theoretical_isotopic_distribution(
            tid, error_tolerance=error_tolerance)
        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer(self.peaklist, eid, tid)
        return IsotopicFitRecord(peak, score, charge, tid, eid)

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON, truncate_after=0.8):
        """Given a set of candidate monoisotopic peaks and charge states, and a PPM error tolerance,
        fit each putative isotopic pattern.

        Calls :meth:`fit_theoretical_distribution` on each candidate.

        If a fit does not satisfy :attr:`scorer` `.reject`, it is discarded. If a fit has only one real peak
        and has a charge state greater than 1, it will also be discarded.

        Parameters
        ----------
        peak_charge_set : set
            The set of candidate (FittedPeak, charge) tuples to try to fit
        error_tolerance : float
            Matching error tolerance
        charge_carrier : float, optional
            The charge carrier to use. Defaults to `PROTON`

        Returns
        -------
        set
            The set of IsotopicFitRecord instances produced
        """
        results = []
        for peak, charge in peak_charge_set:
            if peak.mz < 1:
                continue
            fit = self.fit_theoretical_distribution(
                peak, error_tolerance, charge, charge_carrier, truncate_after)
            fit.missed_peaks = count_placeholders(fit.experimental)
            if len(drop_placeholders(fit.experimental)) == 1 and fit.charge > 1:
                continue
            if self.scorer.reject(fit):
                continue
            results.append(fit)
        return set(results)


try:
    from ms_deisotope._c.deconvoluter_base import DeconvoluterBase, AveragineDeconvoluterBase
except ImportError:
    pass


def charge_range_(lo, hi, step=None):
    sign = -1 if lo < 0 else 1
    abs_lo, abs_hi = abs(lo), abs(hi)
    upper = max(abs_lo, abs_hi)
    lower = min(abs_lo, abs_hi)

    for c in range(upper, lower - 1, -1):
        yield c * sign


class ExhaustivePeakSearchDeconvoluterBase(object):
    """Provides common methods for algorithms which attempt to find a deconvolution for every peak
    in a spectrum. This assumes no dependence between different peaks, instead it relies on subtraction,
    breadth of search, and order of encounter to avoid artefactual fits. This is usually not reasonable,
    so instead please use this class's extension, :class:`PeakDependenceGraphDeconvoluterBase` which can
    express dependence of fits on common resources.

    This class is not meant to be instantiated, but instead used as a mixin for classes that also
    inherit from :class:`DeconvoluterBase` and provide methods `fit_theoretical_distribution`
    and `_fit_peaks_at_charges`

    """
    def _update_charge_bounds_with_prediction(self, peak, charge_range):
        """Update the charge range upper limit in `charge_range` based upon the
        Fourier-Patterson charge state estimate for `peak`

        Parameters
        ----------
        peak : FittedPeak
            The peak to estimate the charge state for
        charge_range : tuple
            The charge state bounds to be updated

        Returns
        -------
        tuple
            The updated charge state bounds
        """
        polarity = 1 if min(charge_range) > 0 else -1
        upper_charge_limit = self.peaklist.predict_charge_state(peak)
        if abs(upper_charge_limit) <= 1:
            return charge_range
        else:
            upper_charge_limit += (2)
            upper_charge_limit *= polarity
        low = min(charge_range, key=abs)
        return (low, upper_charge_limit)

    def _get_all_peak_charge_pairs(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=3,
                                   right_search_limit=3, use_charge_state_hint=False,
                                   recalculate_starting_peak=True):
        """Construct the set of all unique candidate (monoisotopic peak, charge state) pairs using
        the provided search parameters.

        The search is performed using :func:`has_previous_peak_at_charge`, :func:`has_successor_peak_at_charge`,
        :meth:`_find_previous_putative_peak`, and :meth:`_find_next_putative_peak`.

        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        recalculate_starting_peak : bool, optional
            Whether or not to re-calculate the putative starting peak m/z based upon nearby
            peaks close to where isotopic peaks for `peak` should be. Defaults to True

        Returns
        -------
        set
            The set of all unique candidate (monoisotopic peak, charge state)
        """
        if use_charge_state_hint:
            charge_range = self._update_charge_bounds_with_prediction(
                peak, charge_range)

        target_peaks = set()
        if self.verbose:
            info("Considering charge range %r for %r" %
                 (list(charge_range_(*charge_range)), peak))
        for charge in charge_range_(*charge_range):
            target_peaks.add((peak, charge))

            # Look Left
            for i in range(1, left_search_limit):
                prev_peak = has_previous_peak_at_charge(
                    self, peak, charge, i)
                if prev_peak is None:
                    continue
                target_peaks.add((prev_peak, charge))

                if recalculate_starting_peak:
                    target_peaks.update(self._find_previous_putative_peak(
                        peak.mz, charge, i, 2 * error_tolerance))

            # Look Right
            for i in range(1, right_search_limit):
                nxt_peak = has_successor_peak_at_charge(
                    self, peak, charge, i)
                if nxt_peak is None:
                    continue
                target_peaks.add((nxt_peak, charge))

                if recalculate_starting_peak:
                    target_peaks.update(self._find_next_putative_peak(
                        peak.mz, charge, i, 2 * error_tolerance))

            if recalculate_starting_peak:
                for i in range(min(left_search_limit, 2)):
                    target_peaks.update(self._find_next_putative_peak(
                        peak.mz, charge, step=i, tolerance=2 * error_tolerance))

        return target_peaks

    def _fit_all_charge_states(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=3,
                               right_search_limit=3, use_charge_state_hint=False,
                               recalculate_starting_peak=True, charge_carrier=PROTON,
                               truncate_after=0.8):
        """Carry out the fitting process for `peak`.

        This method calls :meth:`_get_all_peak_charge_pairs` to collect all hypothetical solutions
        for `peak`, and invokes :meth:`_fit_peaks_at_charges` to evaluate them.

        The method :meth:`_fit_peaks_at_charges` is required by this interface, but is not defined by
        it, as it depends upon the underlying isotopic pattern fitting algorithm. See one of the
        Averagine-based algorithms for an implementation, such as :class:`AveragineDeconvoluterBase`,
        a complementary ancestor with this class.

        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        recalculate_starting_peak : bool, optional
            Whether or not to re-calculate the putative starting peak m/z based upon nearby
            peaks close to where isotopic peaks for `peak` should be. Defaults to True
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        set
            The set of IsotopicFitRecord instances produced
        """
        target_peaks = self._get_all_peak_charge_pairs(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            use_charge_state_hint=use_charge_state_hint, recalculate_starting_peak=True)

        results = self._fit_peaks_at_charges(
            target_peaks, error_tolerance, charge_carrier=charge_carrier, truncate_after=truncate_after)
        return (results)

    def charge_state_determination(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=3,
                                   right_search_limit=3, use_charge_state_hint=False, charge_carrier=PROTON,
                                   truncate_after=0.8):
        """Determine the optimal isotopic fit for `peak`, extracting it's charge state and monoisotopic peak.

        This method invokes :meth:`_fit_all_charge_states`, and then uses :attr:`scorer`'s `select` method to
        choose the optimal solution.

        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        IsotopicFitRecord
            The best scoring isotopic fit
        """
        results = self._fit_all_charge_states(
            peak, error_tolerance=error_tolerance, charge_range=charge_range, left_search_limit=left_search_limit,
            right_search_limit=right_search_limit, use_charge_state_hint=use_charge_state_hint,
            charge_carrier=charge_carrier, truncate_after=truncate_after)

        if self.verbose:
            info("Fits for %r" % peak)
            for rec in sorted(results)[-10:]:
                info(rec)

        try:
            result = self.scorer.select(results)
            return result
        except ValueError:
            return None

    def deconvolute_peak(self, peak, error_tolerance=2e-5, charge_range=(1, 8), use_charge_state_hint=False,
                         left_search_limit=3, right_search_limit=3, charge_carrier=PROTON, truncate_after=0.8):
        """Perform a deconvolution for `peak`, generating a new :class:`ms_deisotope.peak_set.DeconvolutedPeak` instance
        corresponding to the optimal solution.

        This new peak has an m/z matching the monoisotopic peak of the pattern containing `peak`, and its intensity
        is the sum of all the matched peaks in its isotopic pattern. Its charge, isotopic fit, and other qualities
        are derived from the :class:`ms_deisotope.scoring.IsotopicFitRecord` instance corresponding to its best
        solution.

        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        DeconvolutedPeak
        """
        fit = self.charge_state_determination(
            peak, error_tolerance=error_tolerance,
            charge_range=charge_range, use_charge_state_hint=use_charge_state_hint,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after)
        if fit is None:
            return
        score, charge, eid, tid = fit
        rep_eid = drop_placeholders(eid)
        total_abundance = sum(p.intensity for p in rep_eid)
        monoisotopic_mass = neutral_mass(
            eid[0].mz, charge, charge_carrier=charge_carrier)
        reference_peak = first_peak(eid)

        dpeak = DeconvolutedPeak(
            neutral_mass=monoisotopic_mass, intensity=total_abundance, charge=charge,
            signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
            index=reference_peak.index,
            full_width_at_half_max=mean(
                p.full_width_at_half_max for p in rep_eid),
            a_to_a2_ratio=a_to_a2_ratio(tid),
            most_abundant_mass=neutral_mass(most_abundant_mz(eid), charge),
            average_mass=neutral_mass(average_mz(eid), charge),
            score=score,
            envelope=[(p.mz, p.intensity) for p in rep_eid],
            mz=eid[0].mz,
            fit=fit,
            area=sum(e.area for e in eid))

        self._deconvoluted_peaks.append(dpeak)
        if self.use_subtraction:
            self.subtraction(tid, error_tolerance)
        return dpeak

    def targeted_deconvolution(self, peak, error_tolerance=2e-5, charge_range=(1, 8), use_charge_state_hint=False,
                               left_search_limit=3, right_search_limit=3, charge_carrier=PROTON, truncate_after=0.8):
        """Express the intent that this peak's deconvolution solution will be retrieved at a later point in the process
        and that it should be deconvoluted, and return a handle to retrieve the results with.

        This algorithm's implementation is simple enough that this is equivalent to just performing the deconvolution
        now and storing the result in a :class:`ms_deisotope.utils.TrivialTargetedDeconvolutionResult` instance.

        Otherwise identical to :meth:`deconvolute_peak`.

        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        TrivialTargetedDeconvolutionResult
        """
        dpeak = self.deconvolute_peak(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            use_charge_state_hint=use_charge_state_hint, left_search_limit=left_search_limit,
            right_search_limit=right_search_limit, charge_carrier=charge_carrier, truncate_after=truncate_after)
        result = TrivialTargetedDeconvolutionResult(self, dpeak, peak)
        return result

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8),
                    order_chooser=operator.attrgetter('index'),
                    use_charge_state_hint=False, left_search_limit=3,
                    right_search_limit=3, charge_carrier=PROTON, truncate_after=0.8):
        """Completely deconvolute the spectrum.

        Visit each peak in the order chosen by `order_chooser`, and call :meth:`deconvolute_peak`
        on it with the provided arguments. This assumes all overlaps in isotopic pattern are captured
        by the search limits. This is usually not the case. For an alternative see
        :class:`PeakDependenceGraphDeconvoluterBase`

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        order_chooser : callable, optional:
            A callable used as a key function for sorting peaks into the order they will
            be visited during deconvolution. Defaults to `operator.attrgetter("index")`
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        DeconvolutedPeakSet
        """
        i = 0
        for peak in sorted(self.peaklist, key=order_chooser, reverse=True):
            if peak.mz < 2 or peak.intensity < self.minimum_intensity:
                continue
            self.deconvolute_peak(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                use_charge_state_hint=use_charge_state_hint, left_search_limit=left_search_limit,
                right_search_limit=right_search_limit, charge_carrier=charge_carrier, truncate_after=truncate_after)
            i += 1

        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(
                self._deconvoluted_peaks)

        return DeconvolutedPeakSet(self._deconvoluted_peaks)._reindex()


class AveragineDeconvoluter(AveragineDeconvoluterBase, ExhaustivePeakSearchDeconvoluterBase):
    """A Deconvoluter which uses an :title-reference:`averagine` [1] model to generate theoretical isotopic patterns
    for each peak to consider. Combines :class:`AveragineDeconvoluterBase` and
    :class:`ExhaustivePeakSearchDeconvoluterBase` to create a working Deconvoluter type.


    Attributes
    ----------
    averagine : ms_deisotope.averagine.AveragineCache
        The averagine model and associated theoretical isotopic pattern cache to use
        to build theoretical isotopic patterns.
    peaklist : ms_peak_picker.PeakIndex
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : ms_deisotope.scoring.IsotopicFitterBase
        An object derived from IsotopicFitterBase which can evaluate isotopic fits
    verbose : bool
        How much diagnostic information to provide

    References
    ----------
    [1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
        for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
        Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8
    """
    def __init__(self, peaklist, averagine=None, scorer=penalized_msdeconv,
                 use_subtraction=True, scale_method='sum',
                 verbose=False, **kwargs):
        if averagine is None:
            averagine = AveragineCache(peptide, dict())
        else:
            if not isinstance(averagine, AveragineCache):
                averagine = AveragineCache(averagine, dict())
        self.peaklist = peaklist.clone()
        self.averagine = averagine
        self.scorer = scorer
        self._deconvoluted_peaks = []
        self.verbose = verbose

        super(AveragineDeconvoluter, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks=True)

    def config(self):
        return {
            "scale_method": self.scale_method,
            "use_subtraction": self.use_subtraction,
            "verbose": self.verbose,
            "scorer": self.scorer,
            "averagine": self.averagine
        }


class MultiAveragineDeconvoluterBase(DeconvoluterBase):

    def fit_theoretical_distribution(self, peak, error_tolerance, charge, averagine, charge_carrier=PROTON,
                                     truncate_after=0.8):
        tid = averagine.isotopic_cluster(
            peak.mz, charge, charge_carrier=charge_carrier, truncate_after=truncate_after)
        eid = self.match_theoretical_isotopic_distribution(
            tid, error_tolerance=error_tolerance)
        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer(self.peaklist, eid, tid)
        return IsotopicFitRecord(peak, score, charge, tid, eid)

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON, truncate_after=0.8):
        results = []
        for peak, charge in peak_charge_set:
            for averagine in self.averagine:
                if peak.mz < 1:
                    continue
                fit = self.fit_theoretical_distribution(
                    peak, error_tolerance, charge, averagine, charge_carrier=charge_carrier,
                    truncate_after=truncate_after)
                fit.missed_peaks = count_placeholders(fit.experimental)
                fit.data = averagine
                if len(drop_placeholders(fit.experimental)) == 1 and fit.charge > 1:
                    continue
                if self.scorer.reject(fit):
                    continue
                results.append(fit)
        return set(results)


try:
    from ms_deisotope._c.deconvoluter_base import MultiAveragineDeconvoluterBase
except ImportError:
    pass


class MultiAveragineDeconvoluter(MultiAveragineDeconvoluterBase, ExhaustivePeakSearchDeconvoluterBase):
    """A Deconvoluter which uses multiple :title-reference:`averagine` [1] model to generate theoretical
    isotopic patterns for each peak to consider. Combines :class:`MultiAveragineDeconvoluterBase` and
    :class:`ExhaustivePeakSearchDeconvoluterBase` to create a working Deconvoluter type.

    This differs from :class:`AveragineDeconvoluter`, in that it will produce multiple isotopic fits for
    each (peak, charge) pair. This is advantageous when the isotopic patterns produced by different models
    are sufficiently different enough that they will favor different peak sets.

    Attributes
    ----------
    averagine : list of ms_deisotope.averagine.AveragineCache
        The averagine models and associated theoretical isotopic pattern caches to use
        to build theoretical isotopic patterns.
    peaklist : ms_peak_picker.PeakIndex
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : ms_deisotope.scoring.IsotopicFitterBase
        An object derived from IsotopicFitterBase which can evaluate isotopic fits
    verbose : bool
        How much diagnostic information to provide

    References
    ----------
    [1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
        for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
        Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8
    """
    def __init__(self, peaklist, averagine=None, scorer=penalized_msdeconv,
                 use_subtraction=True, scale_method='sum',
                 merge_isobaric_peaks=True, minimum_intensity=5.,
                 verbose=False, *args, **kwargs):
        self.peaklist = peaklist.clone()
        self.scorer = scorer
        self.use_subtraction = use_subtraction
        self.scale_method = scale_method

        cache_backend = dict
        if averagine is None:
            averagine = [peptide, glycopeptide, glycan]
        averagine = [
            AveragineCache(avg, backend=cache_backend()) if not isinstance(
                avg, AveragineCache) else avg
            for avg in averagine]
        self.averagine = averagine
        self.verbose = verbose

        self._deconvoluted_peaks = []

        super(MultiAveragineDeconvoluter, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)


class PeakDependenceGraphDeconvoluterBase(ExhaustivePeakSearchDeconvoluterBase):
    """Extends the concept of :class:`ExhaustivePeakSearchDeconvoluterBase` to include a way to handle
    conflicting solutions which claim the same experimental peak.

    This lets the Deconvoluter assign a single peak only once, and to the "best" solution to use it. To
    do this, the Deconvoluter constructs a graph where peaks are nodes, and isotopic fits are hyperedges
    connecting multiple nodes. Rather than deconvoluting the spectrum step by step, assigning signal as
    it explores the spectrum, the Deconvoluter instead inserts each considered isotopic fit into the graph.
    After completely traversing the spectrum, the Deconvoluter solves the dependence graph attempting to
    maximize some criterion and produce a set of disjoint isotopic fits. These fits are then assigned signal
    and added to the deconvoluted spectrum as normal.

    The criterion used is currently a greedy maximization (or minimization) of each connected component of
    the peak dependence graph.

    Attributes
    ----------
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : PeakDependenceGraph
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
    """
    def __init__(self, peaklist, *args, **kwargs):
        max_missed_peaks = kwargs.pop("max_missed_peaks", 1)
        ExhaustivePeakSearchDeconvoluterBase.__init__(self)
        self.peak_dependency_network = PeakDependenceGraph(
            self.peaklist, maximize=self.scorer.is_maximizing())
        self.max_missed_peaks = max_missed_peaks
        self._priority_map = {}

    @property
    def max_missed_peaks(self):
        return self.peak_dependency_network.max_missed_peaks

    @max_missed_peaks.setter
    def max_missed_peaks(self, value):
        self.peak_dependency_network.max_missed_peaks = value

    def _explore_local(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, use_charge_state_hint=False, charge_carrier=PROTON,
                       truncate_after=0.8):
        """Given a peak, explore the local neighborhood for candidate isotopic fits and add each
        fit above a threshold to the peak dependence graph.

        The threshold assumes that a single peak's neighborhood will contain many, many fits, but
        that only the top `n` scoring fits are worth considering. For now, `n` is fixed at 100 or
        the half number of fits returned, whichever is larger. This is to prevent the fit graph
        from growing out of control and wasting time storing impractical fits. Any fit added to
        the graph will have to pass :attr:`scorer.select` as well, so weak fits will never be added,
        regardless of how many fits are allowed to be inserted.
    
        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        
        Returns
        -------
        int
            The number of fits added to the graph
        """
        results = self._fit_all_charge_states(
            peak, error_tolerance=error_tolerance, charge_range=charge_range, left_search_limit=left_search_limit,
            right_search_limit=right_search_limit,
            use_charge_state_hint=use_charge_state_hint, charge_carrier=charge_carrier, truncate_after=truncate_after)

        hold = set()
        for fit in results:
            if fit.charge > 1 and len(drop_placeholders(fit.experimental)) == 1:
                continue
            hold.add(fit)

        results = hold

        n = len(results)
        stop = max(min(n / 2, 100), 10)
        if n == 0:
            return 0

        if self.verbose:
            info("\nFits for %r (%f)" % (peak, peak.mz))

        for i in range(stop):
            if len(results) == 0:
                break
            candidate = self.scorer.select(results)
            if self.verbose:
                info("Candidate: %r", candidate)
            self.peak_dependency_network.add_fit_dependence(candidate)
            results.discard(candidate)

        return i

    def populate_graph(self, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, use_charge_state_hint=False, charge_carrier=PROTON,
                       truncate_after=0.8):
        """Visit each experimental peak and execute :meth:`_explore_local` on it with the provided
        parameters, populating the peak dependence graph with all viable candidates.

        Parameters
        ----------
        peak : FittedPeak
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        """
        for peak in self.peaklist:
            if peak in self._priority_map:
                continue
            out = self._explore_local(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                use_charge_state_hint=use_charge_state_hint, charge_carrier=charge_carrier,
                truncate_after=truncate_after)

    def postprocess_fits(self, fit_postprocessor=None, error_tolerance=2e-5, charge_range=(1, 8),
                         charge_carrier=PROTON, *args, **kwargs):
        if fit_postprocessor is None:
            return

    def select_best_disjoint_subgraphs(self, error_tolerance=2e-5, charge_carrier=PROTON):
        """Construct connected envelope graphs from :attr:`peak_dependency_network` and
        extract the best disjoint isotopic pattern fits in each envelope graph. This in turn
        produces one or more :class:`DeconvolutedPeak` instances from each disjoint fit,
        which are processed and added to the results set.
        
        Parameters
        ----------
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result
        """
        disjoint_envelopes = self.peak_dependency_network.find_non_overlapping_intervals()
        i = 0
        for cluster in disjoint_envelopes:
            disjoint_best_fits = cluster.disjoint_best_fits()
            for fit in disjoint_best_fits:
                score, charge, eid, tid = fit
                rep_eid = drop_placeholders(eid)
                if len(rep_eid) == 0:
                    continue
                total_abundance = sum(p.intensity for p in rep_eid)
                monoisotopic_mass = neutral_mass(
                    eid[0].mz, charge, charge_carrier)
                reference_peak = first_peak(eid)

                dpeak = DeconvolutedPeak(
                    neutral_mass=monoisotopic_mass, intensity=total_abundance, charge=charge,
                    signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
                    index=reference_peak.index,
                    full_width_at_half_max=mean(
                        p.full_width_at_half_max for p in rep_eid),
                    a_to_a2_ratio=a_to_a2_ratio(tid),
                    most_abundant_mass=neutral_mass(
                        most_abundant_mz(eid), charge),
                    average_mass=neutral_mass(average_mz(eid), charge),
                    score=score,
                    envelope=[(p.mz, p.intensity) for p in rep_eid],
                    mz=eid[0].mz, fit=fit,
                    area=sum(e.area for e in eid))

                self.peak_dependency_network.add_solution(fit, dpeak)

                self._deconvoluted_peaks.append(dpeak)
                i += 1
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)

    def targeted_deconvolution(self, peak, error_tolerance=2e-5, charge_range=(1, 8), use_charge_state_hint=False,
                               left_search_limit=3, right_search_limit=3, charge_carrier=PROTON, truncate_after=0.8):
        """Express the intent that this peak's deconvolution solution will be retrieved at a later point in the process
        and that it should be deconvoluted, and return a handle to retrieve the results with.

        As the operation does not immediately result in a deconvoluted peak but just adds the resulting fits to
        :attr:`peak_dependency_network`, this method constructs an instance of
        :class:`ms_deisotope.peak_dependency_network.peak_network.NetworkedTargetedDeconvolutionResult` which holds all
        the required information for recovering the best fit containing `peak`.

        Parameters
        ----------
        peak : FittedPeak
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`

        Returns
        -------
        ms_deisotope.peak_dependency_network.peak_network.NetworkedTargetedDeconvolutionResult
        """
        self._explore_local(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            use_charge_state_hint=use_charge_state_hint, left_search_limit=left_search_limit,
            right_search_limit=right_search_limit, charge_carrier=charge_carrier, truncate_after=truncate_after)
        result = NetworkedTargetedDeconvolutionResult(self, peak)
        self._priority_map[peak] = result
        return result

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8),
                    use_charge_state_hint=False, left_search_limit=1,
                    right_search_limit=0, iterations=1, charge_carrier=PROTON,
                    truncate_after=0.8, fit_postprocessor=None):
        """Completely deconvolute the spectrum.

        For each iteration, clear :attr:`peak_depencency_network`, then invoke :meth:`populate_graph`
        followed by :meth:`select_best_disjoint_subgraphs` to populate the resulting :class:`DeconvolutedPeakSet`

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to 2e-5
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        order_chooser : callable, optional:
            A callable used as a key function for sorting peaks into the order they will
            be visited during deconvolution. Defaults to `operator.attrgetter("index")`
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        use_charge_state_hint : bool, optional
            Whether or not to try to estimate the upper limit of the charge states to consider
            using :meth:`_update_charge_bounds_with_prediction`. Defaults to False
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to `PROTON`
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        DeconvolutedPeakSet
        """
        for i in range(iterations):
            self.peak_dependency_network.reset()
            self.populate_graph(
                error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                use_charge_state_hint=use_charge_state_hint, charge_carrier=charge_carrier,
                truncate_after=truncate_after)
            self.postprocess_fits(
                fit_postprocessor=fit_postprocessor, charge_range=charge_range, charge_carrier=charge_carrier,
                error_tolerance=error_tolerance)
            self.select_best_disjoint_subgraphs(
                error_tolerance, charge_carrier)

        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(
                self._deconvoluted_peaks)

        return DeconvolutedPeakSet(list(self._deconvoluted_peaks))._reindex()


class AveraginePeakDependenceGraphDeconvoluter(AveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    """Extends :class:`AveragineDeconvoluter` to include features from
    :class:`PeakDependenceGraphDeconvoluterBase` making it suitable for deconvoluting complex spectra where
    peak overlaps are common.
    """
    def __init__(self, peaklist, *args, **kwargs):
        AveragineDeconvoluter.__init__(self, peaklist, *args, **kwargs)
        PeakDependenceGraphDeconvoluterBase.__init__(self, peaklist, **kwargs)


class MultiAveraginePeakDependenceGraphDeconvoluter(MultiAveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    """Extends :class:`MultiAveragineDeconvoluter` to include features from
    :class:`PeakDependenceGraphDeconvoluterBase` making it suitable for deconvoluting complex spectra where
    peak overlaps are common.
    """
    def __init__(self, peaklist, *args, **kwargs):
        MultiAveragineDeconvoluter.__init__(self, peaklist, *args, **kwargs)
        PeakDependenceGraphDeconvoluterBase.__init__(self, peaklist, **kwargs)


class CompositionListDeconvoluterBase(object):
    """A mixin class to provide common features for deconvoluters which process spectra
    using a list of targeted compositions.
    
    Attributes
    ----------
    composition_list : Sequence of Mapping
        A series of objects which represent elemental compositions and support
        the Mapping interface to access their individual elements.
    """
    def __init__(self, composition_list):
        self.composition_list = composition_list

    def generate_theoretical_isotopic_cluster(self, composition, charge, truncate_after=0.8,
                                              mass_shift=None, charge_carrier=PROTON):
        """Generate a theoretical isotopic pattern for `compos
        
        Parameters
        ----------
        composition : Mapping
            An object representing an elemental composition
        charge : int
            The charge state to generate the isotopic pattern for
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift : float, optional
            An arbitrary mass shift to apply to the generated theoretical isotopic pattern,
            moving all peaks forward by that mass charge ratio transformed mass.
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to `PROTON`
        
        Returns
        -------
        list of TheoreticalPeak
            The theoretical isotopic pattern generated
        """
        cumsum = 0
        result = []
        for peak in isotopic_variants(composition, charge=charge, charge_carrier=charge_carrier):
            cumsum += peak.intensity
            result.append(peak)
            if cumsum >= truncate_after:
                break
        if mass_shift is not None:
            shift_isotopic_pattern(mass_shift / abs(charge), result)
        return result

    def fit_composition_at_charge(self, composition, charge, error_tolerance=2e-5, charge_carrier=PROTON,
                                  truncate_after=0.8, mass_shift=None):
        """Produce an isotopic fit for `composition` at `charge` against the experimental peak set.

        This method requires that the instance also possess a method named `match_theoretical_isotopic_distribution`
        such as the one implemented in :class:`DeconvoluterBase`.

        Parameters
        ----------
        composition : Mapping
            An object representing an elemental composition
        charge : int
            The charge state to generate the isotopic pattern for
        error_tolerance : float
            The mass accuracy required to for peak matches
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift : float, optional
            An arbitrary mass shift to apply to the generated theoretical isotopic pattern,
            moving all peaks forward by that mass charge ratio transformed mass.
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to `PROTON`

        Returns
        -------
        ms_deisotope.scoring.IsotopicFitRecord
        """
        tid = self.generate_theoretical_isotopic_cluster(composition, charge=charge, truncate_after=truncate_after,
                                                         charge_carrier=charge_carrier, mass_shift=mass_shift)
        eid = self.match_theoretical_isotopic_distribution(
            tid, error_tolerance)

        missed_peaks = count_placeholders(eid)

        if missed_peaks > len(eid) / 2:
            return None

        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer.evaluate(self.peaklist, eid, tid)
        fit = IsotopicFitRecord(None, score, charge, tid, eid)
        fit.missed_peaks = missed_peaks
        return fit

    def deconvolute_composition(self, composition, error_tolerance=2e-5, charge_range=(1, 8), charge_carrier=PROTON,
                                truncate_after=0.8, mass_shift=None):
        """For each charge state under consideration, fit the theoretical isotopic pattern for this composition,
        and if the fit is satisfactory, add it to the results set.

        Parameters
        ----------
        composition : Mapping
            An object representing an elemental composition
        error_tolerance : float
            The mass accuracy required to for peak matches
        charge_range : tuple
            The charge state range to generate the isotopic patterns for
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift : float, optional
            An arbitrary mass shift to apply to the generated theoretical isotopic pattern,
            moving all peaks forward by that mass charge ratio transformed mass.
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to `PROTON`
        """
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(composition, charge=charge, error_tolerance=error_tolerance,
                                                 truncate_after=truncate_after, charge_carrier=charge_carrier,
                                                 mass_shift=mass_shift)
            if fit is None:
                continue
            if not self.scorer.reject(fit):
                eid = fit.experimental
                tid = fit.theoretical

                total_abundance = sum(
                    p.intensity for p in eid if p.intensity > 1)

                monoisotopic_mass = neutral_mass(
                    eid[0].mz, charge, charge_carrier)
                monoisotopic_mz = eid[0].mz

                reference_peak = first_peak(eid)
                rep_eid = drop_placeholders(eid)
                if len(rep_eid) < 2 or len(rep_eid) < len(tid) / 2. or len(rep_eid) == 1 and fit.charge > 1:
                    continue
                peak = DeconvolutedPeakSolution(
                    composition, fit,
                    monoisotopic_mass, total_abundance, charge,
                    signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
                    index=reference_peak.index,
                    full_width_at_half_max=mean(
                        p.full_width_at_half_max for p in rep_eid),
                    a_to_a2_ratio=a_to_a2_ratio(tid),
                    most_abundant_mass=neutral_mass(
                        most_abundant_mz(eid), charge),
                    average_mass=neutral_mass(average_mz(eid), charge),
                    score=fit.score,
                    envelope=[(p.mz, p.intensity) for p in rep_eid],
                    mz=monoisotopic_mz, area=sum(e.area for e in eid))

                self._deconvoluted_peaks.append(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)


class CompositionListDeconvoluter(CompositionListDeconvoluterBase, DeconvoluterBase):

    def __init__(self, peaklist, composition_list, scorer,
                 use_subtraction=False, scale_method='sum',
                 verbose=False):
        self.peaklist = peaklist.clone()
        self.scorer = scorer
        self.verbose = verbose
        self._deconvoluted_peaks = []
        CompositionListDeconvoluterBase.__init__(self, composition_list)
        DeconvoluterBase.__init__(
            self,
            use_subtraction=use_subtraction, scale_method=scale_method, merge_isobaric_peaks=True)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), charge_carrier=PROTON, truncate_after=0.8,
                    mass_shift=None, **kwargs):
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance=error_tolerance,
                                         charge_range=charge_range, charge_carrier=charge_carrier,
                                         truncate_after=truncate_after, mass_shift=mass_shift)
        return DeconvolutedPeakSet(self._deconvoluted_peaks)._reindex()


class CompositionListPeakDependenceGraphDeconvoluter(CompositionListDeconvoluter):

    def __init__(self, peaklist, composition_list, scorer,
                 use_subtraction=False, scale_method='sum',
                 verbose=False, **kwargs):
        max_missed_peaks = kwargs.get("max_missed_peaks", 1)
        super(CompositionListPeakDependenceGraphDeconvoluter, self).__init__(
            peaklist, composition_list, scorer, use_subtraction, scale_method,
            verbose)

        self.peak_dependency_network = PeakDependenceGraph(
            self.peaklist, maximize=self.scorer.is_maximizing(), **kwargs)
        self.max_missed_peaks = max_missed_peaks

    @property
    def max_missed_peaks(self):
        return self.peak_dependency_network.max_missed_peaks

    @max_missed_peaks.setter
    def max_missed_peaks(self, value):
        self.peak_dependency_network.max_missed_peaks = value

    def _save_peak_solution(self, solution):
        self._deconvoluted_peaks.append(solution)

    def deconvolute_composition(self, composition, error_tolerance=2e-5, charge_range=(1, 8), truncate_after=0.8,
                                charge_carrier=PROTON, mass_shift=None):
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(
                composition, charge, error_tolerance, charge_carrier=charge_carrier,
                truncate_after=truncate_after, mass_shift=mass_shift)
            if fit is None:
                continue
            rep_eid = drop_placeholders(fit.experimental)
            if len(rep_eid) == 1 and fit.charge > 1:
                continue
            if not self.scorer.reject(fit):
                self.peak_dependency_network.add_fit_dependence(fit)

    def populate_graph(self, error_tolerance=2e-5, charge_range=(1, 8), truncate_after=0.8, charge_carrier=PROTON,
                       mass_shift=None):
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance, charge_range,
                                         truncate_after=truncate_after, charge_carrier=charge_carrier,
                                         mass_shift=mass_shift)

    def select_best_disjoint_subgraphs(self, error_tolerance=2e-5, charge_carrier=PROTON):
        disjoint_envelopes = self.peak_dependency_network.find_non_overlapping_intervals()

        for cluster in disjoint_envelopes:
            for fit in cluster.disjoint_best_fits():
                eid = fit.experimental
                tid = fit.theoretical
                composition = fit.data
                charge = fit.charge
                score = fit.score

                total_abundance = sum(
                    p.intensity for p in eid if p.intensity > 1)
                monoisotopic_mass = neutral_mass(
                    eid[0].mz, charge, charge_carrier)
                monoisotopic_mz = eid[0].mz
                reference_peak = first_peak(eid)
                rep_eid = drop_placeholders(eid)
                if len(rep_eid) < 2 or len(rep_eid) < len(tid) / 2.:
                    continue
                peak = DeconvolutedPeakSolution(
                    composition, fit,
                    monoisotopic_mass, total_abundance, charge,
                    signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
                    index=reference_peak.index,
                    full_width_at_half_max=mean(
                        p.full_width_at_half_max for p in rep_eid),
                    a_to_a2_ratio=a_to_a2_ratio(tid),
                    most_abundant_mass=neutral_mass(
                        most_abundant_mz(eid), charge),
                    average_mass=neutral_mass(average_mz(eid), charge),
                    score=score,
                    envelope=[(p.mz, p.intensity) for p in rep_eid],
                    mz=monoisotopic_mz, area=sum(e.area for e in eid))

                self._save_peak_solution(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), iterations=1, truncate_after=0.8,
                    charge_carrier=PROTON, mass_shift=None, **kwargs):
        for i in range(iterations):
            self.populate_graph(error_tolerance, charge_range, charge_carrier=charge_carrier,
                                truncate_after=truncate_after, mass_shift=mass_shift)
            self.select_best_disjoint_subgraphs(error_tolerance)
        return DeconvolutedPeakSet(self._deconvoluted_peaks)._reindex()


class HybridAveragineCompositionListPeakDependenceGraphDeconvoluter(
        AveraginePeakDependenceGraphDeconvoluter, CompositionListDeconvoluterBase):
    def __init__(self, peaklist, composition_list, *args, **kwargs):
        AveraginePeakDependenceGraphDeconvoluter.__init__(self, peaklist, *args, **kwargs)
        CompositionListDeconvoluterBase.__init__(self, composition_list)

    def deconvolute_composition(self, composition, error_tolerance=2e-5, charge_range=(1, 8), truncate_after=0.8,
                                charge_carrier=PROTON, mass_shift=None):
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(
                composition, charge, error_tolerance, charge_carrier=charge_carrier,
                truncate_after=truncate_after, mass_shift=mass_shift)
            if fit is None:
                continue
            rep_eid = drop_placeholders(fit.experimental)
            if len(rep_eid) == 1 and fit.charge > 1:
                continue
            if not self.scorer.reject(fit):
                fit.data = (composition, charge)
                self.peak_dependency_network.add_fit_dependence(fit)

    def populate_graph(self, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, use_charge_state_hint=False, charge_carrier=PROTON,
                       truncate_after=0.8, mass_shift=None):
        for composition in self.composition_list:
            self.deconvolute_composition(
                composition, error_tolerance, charge_range=charge_range,
                truncate_after=truncate_after, charge_carrier=charge_carrier,
                mass_shift=mass_shift)
        AveraginePeakDependenceGraphDeconvoluter.populate_graph(
            self, error_tolerance,
            charge_range=charge_range,
            left_search_limit=left_search_limit,
            right_search_limit=right_search_limit,
            use_charge_state_hint=use_charge_state_hint,
            charge_carrier=charge_carrier, truncate_after=truncate_after)

    def subtraction(self, isotopic_cluster, error_tolerance):
        super(HybridAveragineCompositionListPeakDependenceGraphDeconvoluter, self).subtraction(
            isotopic_cluster, error_tolerance)


def deconvolute_peaks(peaklist, deconvoluter_type=AveraginePeakDependenceGraphDeconvoluter,
                      decon_config=None, charge_range=(1, 8), error_tolerance=2e-5, priority_list=None,
                      use_charge_state_hint_for_priorities=False, left_search_limit=3, right_search_limit=3,
                      left_search_limit_for_priorities=None, right_search_limit_for_priorities=None,
                      verbose_priorities=False, verbose=False, charge_carrier=PROTON, truncate_after=0.8, **kwargs):
    if priority_list is None:
        priority_list = []
    if left_search_limit_for_priorities is None:
        left_search_limit_for_priorities = left_search_limit
    if right_search_limit_for_priorities is None:
        right_search_limit_for_priorities = right_search_limit

    decon_config = decon_config or {}
    decon_config.update(kwargs)
    decon_config.setdefault("use_subtraction", False)
    decon = deconvoluter_type(peaklist=peaklist, **decon_config)

    if verbose_priorities or verbose:
        decon.verbose = True

    priority_list_results = []
    for p in priority_list:
        try:
            target_info = p
            p = target_info.peak
            hinted_charge_range = target_info.charge_range_hint(charge_range)
        except AttributeError:
            hinted_charge_range = charge_range
        if not isinstance(p, FittedPeak):
            p = decon.peaklist.has_peak(p, error_tolerance)
        priority_result = decon.targeted_deconvolution(
            p, error_tolerance=error_tolerance, charge_range=hinted_charge_range,
            use_charge_state_hint=use_charge_state_hint_for_priorities,
            left_search_limit=left_search_limit_for_priorities,
            right_search_limit=right_search_limit_for_priorities,
            charge_carrier=charge_carrier,
            truncate_after=truncate_after)
        priority_list_results.append(priority_result)

    if verbose_priorities and not verbose:
        decon.verbose = False

    deconvoluted_peaks = decon.deconvolute(
        error_tolerance=error_tolerance, charge_range=charge_range, left_search_limit=left_search_limit,
        right_search_limit=right_search_limit, charge_carrier=charge_carrier, truncate_after=truncate_after)

    acc = []
    for pr in priority_list_results:
        try:
            result = pr.get()
        except ValueError:
            result = None
            logger.info("Could not extract a solution for %r", pr.query_peak, exc_info=True)
        acc.append(result)

    priority_list_results = acc

    return DeconvolutionProcessResult(decon, deconvoluted_peaks, priority_list_results)
