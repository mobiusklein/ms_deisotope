# -*- coding: utf-8 -*-
'''A base class for deconvolution strategies.
'''
import operator
import logging

from ms_peak_picker import FittedPeak


from ms_deisotope.scoring import IsotopicFitRecord
from ms_deisotope.constants import (
    ERROR_TOLERANCE)

from ms_deisotope.utils import (
    Base)

from .utils import (
    isotopic_shift, drop_placeholders)


logger = logging.getLogger("deconvolution")
info = logger.info
debug = logger.debug
error = logger.error


class DeconvoluterBase(Base):
    """Base class for all Deconvoluter types. Provides basic configuration for common operations,
    regardless of implementation. Because these methods form the backbone of all deconvolution algorithms,
    this class has a C-extension implementation as well.

    Attributes
    ----------
    peaklist : ms_peak_picker.PeakSet
        The centroided mass spectrum to deconvolute
    scorer : IsotopicFitterBase
        The criterion for evaluating individual isotopic pattern fits
    merge_isobaric_peaks : bool
        If multiple passes produce peaks with identical mass values,
        should those peaks be summed
    minimum_intensity : float
        Experimental peaks whose intensity is below this level will be ignored
        by peak querying methods
    scale_method : str
        The name of the method to use to scale theoretical isotopic pattern intensities
        to match the experimental isotopic pattern. For a description of options, see
        :meth:`~.TheoreticalIsotopicPattern.scale`. The default method is `"sum"`.
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
    peaklist = None
    scorer = None

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
            return FittedPeak(mz, 1.0, 1.0, -1, 0, 0, 0)
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

    def match_theoretical_isotopic_distribution(self, theoretical_distribution, error_tolerance=ERROR_TOLERANCE):
        """Given a list of theoretical peaks, find their counterparts in :attr:`peaklist` within `error_tolerance`
        ppm error. If no experimental peak is found, a placeholder will be used in its stead.

        Parameters
        ----------
        theoretical_distribution : list of :class:`~.TheoreticalPeak`
            The theoretical isotopic pattern to match
        error_tolerance : float, optional
            Parts-per-million error tolerance to permit in searching for matches

        Returns
        -------
        list of :class:`~.FittedPeak`
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

        The scaling algorithm used is controlled by :attr:`scale_method`. For a description of options, see
        :meth:`~.TheoreticalIsotopicPattern.scale`.

        Parameters
        ----------
        theoretical_distribution : :class:`~.TheoreticalIsotopicPattern`
            The theoretical isotopic pattern to scale up
        experimental_distribution : list of :class:`~.FittedPeak`
            The experimental isotopic pattern to use as a reference

        Returns
        -------
        list of :class:`~.TheoreticalPeak`
        """
        theoretical_distribution.scale(
            experimental_distribution, self.scale_method)

    def _evaluate_theoretical_distribution(self, experimental, theoretical, peak, charge):
        """Evaluate a provided theoretical isotopic pattern fit against a
        set of matched experimental peaks.

        Parameters
        ----------
        experimental : list
            A list of experimental fitted peaks to test
        theoretical : :class:`~.TheoreticalIsotopicPattern`
            A list of theoretical isotopic peaks to test
        peak : :class:`~.FittedPeak`
            The seed peak for the isotopic pattern fit
        charge : int
            The target charge state

        Returns
        -------
        :class:`~.IsotopicFitRecord`
        """
        self.scale_theoretical_distribution(theoretical, experimental)
        score = self.scorer(self.peaklist, experimental, theoretical)  # pylint: disable=not-callable
        return IsotopicFitRecord(peak, score, charge, theoretical, experimental)

    def fit_incremental_truncation(self, seed_fit, lower_bound):
        """Fit incrementally truncated versions of the seed fit to check to see if a narrower
        theoretical fit matches the data better.

        Parameters
        ----------
        seed_fit : :class:`~.IsotopicFitRecord`
            The original fit to explore truncations of
        lower_bound : float
            The percentage of total signal remaining to stop truncating at.

        Returns
        -------
        :class:`list` of :class:`~.IsotopicFitRecord`
        """
        patterns = seed_fit.theoretical.incremental_truncation(lower_bound)
        fits = [seed_fit]
        for pattern in patterns[1:]:
            k = len(pattern)
            eid = seed_fit.experimental[:k]
            fit = self._evaluate_theoretical_distribution(
                eid, pattern, seed_fit.seed_peak, seed_fit.charge)
            fits.append(fit)
        return fits

    def subtraction(self, isotopic_cluster, error_tolerance=ERROR_TOLERANCE):
        """Subtract signal attributed to `isotopic_cluster` from the equivalent
        peaks in :attr:`peaklist`, mutating the peaks within.

        This will change the intensity value of the matched `FittedPeak` instances,
        and this is reflected in all of there references.

        Parameters
        ----------
        isotopic_cluster : list of :class:`~.TheoreticalPeak`
            The isotopic cluster to subtract
        error_tolerance : float, optional
            Parts-per-million mass accuracy error tolerance to permit when
            finding matches for `isotopic_cluster`
        """
        for peak in isotopic_cluster:
            match = self.peaklist.has_peak(peak.mz, error_tolerance)
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
            if abs(current_peak.neutral_mass - peak.neutral_mass) < 1e-3 and current_peak.charge == peak.charge:
                current_peak.intensity += peak.intensity
            else:
                merged_peaks.append(current_peak)
                current_peak = peak
        merged_peaks.append(current_peak)
        return merged_peaks

    def _find_next_putative_peak(self, mz, charge, step=1, tolerance=ERROR_TOLERANCE):
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
            dummy_peak = FittedPeak(prev_peak_mz, 1.0, 1.0, -1, 0, 0, 0)
            candidates.append((dummy_peak, charge))
        return candidates

    def _find_previous_putative_peak(self, mz, charge, step=1, tolerance=ERROR_TOLERANCE):
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

    def _check_fit(self, fit):
        if len(drop_placeholders(fit.experimental)) == 1 and fit.charge > 1:
            return False
        if self.scorer.reject(fit):
            return False
        return True

    def __repr__(self):
        type_name = self.__class__.__name__
        return "%s(peaklist=%s)" % (type_name, self.peaklist)


try:
    _has_c = True
    _DeconvoluterBase = DeconvoluterBase
    from ms_deisotope._c.deconvoluter_base import DeconvoluterBase
except ImportError as e:
    _has_c = False
