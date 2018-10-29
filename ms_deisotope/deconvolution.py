# -*- coding: utf-8 -*-
import operator
import logging
import numpy as np

from ms_peak_picker import (
    FittedPeak, PeakSet, PeakIndex, simple_peak, is_peak)

from .averagine import (
    AveragineCache, peptide, glycopeptide, glycan, neutral_mass, isotopic_variants,
    isotopic_shift, PROTON, TheoreticalIsotopicPattern)
from .peak_set import DeconvolutedPeak, DeconvolutedPeakSolution, DeconvolutedPeakSet, Envelope
from .scoring import IsotopicFitRecord, penalized_msdeconv, PenalizedMSDeconVFitter, DotProductFitter
from .utils import range, Base, TrivialTargetedDeconvolutionResult, DeconvolutionProcessResult
from .envelope_statistics import a_to_a2_ratio, average_mz, most_abundant_mz
from .peak_dependency_network import PeakDependenceGraph, NetworkedTargetedDeconvolutionResult
from .constants import (
    TRUNCATE_AFTER,
    MAX_ITERATION,
    ERROR_TOLERANCE,
    IGNORE_BELOW,
    CONVERGENCE,
    SCALE_METHOD)

logger = logging.getLogger("deconvolution")
info = logger.info
debug = logger.debug
error = logger.error


def prepare_peaklist(peaks):
    if isinstance(peaks, PeakIndex):
        peaks = PeakSet(peaks.peaks).clone()
    else:
        peaks = tuple(peaks)
        if len(peaks) == 0:
            return PeakSet([])
        if not isinstance(peaks[0], FittedPeak):
            if is_peak(peaks[0]):
                peaks = [simple_peak(p.mz, p.intensity, 0.01) for p in peaks]
            elif isinstance(peaks[0], (list, tuple)):
                peaks = [simple_peak(p[0], p[1], 0.01) for p in peaks]
            else:
                raise TypeError("Cannot convert peaks into a PeakSet")

        peaks = PeakSet(peaks).clone()
    peaks.reindex()
    return peaks


def mean(numbers):
    n = 0.
    total = 0
    for x in numbers:
        n += 1.
        total += x
    return total / n


def has_previous_peak_at_charge(peak_collection, peak, charge=2, step=1, error_tolerance=ERROR_TOLERANCE):
    """Get the `step`th *preceding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_collection : DeconvoluterBase
        Peak collection to look up peaks in. Calls :meth:`has_peak`
    peak : ms_peak_picker.FittedPeak
        The peak to use as a point of reference
    charge : int, optional
        The charge state to interpolate from. Defaults to 2.
    step : int, optional
        The number of peaks along the isotopic pattern to search.

    Returns
    -------
    FittedPeak
    """
    prev = peak.mz - isotopic_shift(charge) * step
    return peak_collection.has_peak(prev, error_tolerance)


def has_successor_peak_at_charge(peak_collection, peak, charge=2, step=1, error_tolerance=ERROR_TOLERANCE):
    """Get the `step`th *succeeding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_collection : DeconvoluterBase
        Peak collection to look up peaks in. Calls :meth:`has_peak`
    peak : ms_peak_picker.FittedPeak
        The peak to use as a point of reference
    charge : int, optional
        The charge state to interpolate from. Defaults to 2.
    step : int, optional
        The number of peaks along the isotopic pattern to search.

    Returns
    -------
    FittedPeak
    """
    nxt = peak.mz + isotopic_shift(charge) * step
    return peak_collection.has_peak(nxt, error_tolerance)


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


def from_fitted_peak(peak, charge=1):
    """Convert a :class:`~.FittedPeak` into a :class:`~.DeconvolutedPeak`
    at the specified charge state.

    Parameters
    ----------
    peak : :class:`~.FittedPeak`
        The fitted peak to use as the template
    charge : int, optional
        The charge state to use, defaults to 1+

    Returns
    -------
    :class:`~.DeconvolutedPeak`
    """
    mass = neutral_mass(peak.mz, charge)
    dpeak = DeconvolutedPeak(
        mass, peak.intensity, charge,
        peak.signal_to_noise, -1, peak.full_width_at_half_max,
        0, mass, mass, 0, Envelope([(peak.mz, peak.intensity)]),
        peak.mz, area=peak.area)
    return dpeak


def quick_charge(peak_set, index, min_charge, max_charge):
    """An implementation of Hoopman's QuickCharge [1] algorithm for quickly capping charge
    state queries

    Parameters
    ----------
    peak_set : :class:`ms_peak_picker.PeakSet
        The centroided peak set to search
    index : int
        The index of the peak to start the search from
    min_charge : int
        The minimum charge state to consider
    max_charge : int
        The maximum charge state to consider

    Returns
    -------
    np.ndarray
        The list of feasible charge states

    References
    ----------
    [1] Hoopmann, M. R., Finney, G. L., MacCoss, M. J., Michael R. Hoopmann, Gregory L. Finney,
        and, MacCoss*, M. J., … MacCoss, M. J. (2007). "High-speed data reduction, feature detection
        and MS/MS spectrum quality assessment of shotgun proteomics data sets using high-resolution
        Mass Spectrometry". Analytical Chemistry, 79(15), 5620–5632. https://doi.org/10.1021/ac0700833
    """
    min_intensity = peak_set[index].intensity / 4.
    charges = np.zeros(max_charge, dtype=int)
    for j in range(index + 1, len(peak_set)):
        if peak_set[j].intensity < min_intensity:
            continue
        diff = peak_set[j].mz - peak_set[index].mz
        if diff > 1.1:
            break
        raw_charge = 1 / diff
        charge = int(raw_charge + 0.5)
        remain = raw_charge - int(raw_charge)
        if 0.2 < remain < 0.8:
            continue
        if charge < min_charge or charge > max_charge:
            continue
        charges[charge] = 1
    if not np.any(charges):
        return np.array([], dtype=int)
    for j in range(index - 1, -1, -1):
        diff = peak_set[index].mz - peak_set[j].mz
        if diff > 1.1:
            break
        raw_charge = 1 / diff
        charge = int(raw_charge + 0.5)
        remain = raw_charge - int(raw_charge)
        if 0.2 < remain < 0.8:
            continue
        if charge < min_charge or charge > max_charge:
            continue
        charges[charge] = 1
    return np.where(charges)[0]


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

        The scaling algorithm used is controlled by :attr:`scale_method`

        Parameters
        ----------
        theoretical_distribution : list of :class:`~.TheoreticalPeak`
            The theoretical isotopic pattern to scale up
        experimental_distribution : list of :class:`~.FittedPeak`
            The experimental isotopic pattern to use as a reference

        Returns
        -------
        list of :class:`~.TheoreticalPeak`
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
        score = self.scorer(self.peaklist, experimental, theoretical)
        return IsotopicFitRecord(peak, score, charge, theoretical, experimental)

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


class AveragineDeconvoluterBase(DeconvoluterBase):
    """A base class derived from :class:`DeconvoluterBase` which provides some common methods
    for fitting isotopic patterns using an Averagine model.
    """
    def __init__(self, use_subtraction=False, scale_method="sum", merge_isobaric_peaks=True,
                 minimum_intensity=5., *args, **kwargs):
        super(AveragineDeconvoluterBase, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)

    def fit_theoretical_distribution(self, peak, error_tolerance, charge, charge_carrier=PROTON, truncate_after=0.8,
                                     ignore_below=IGNORE_BELOW):
        """Fit an isotopic pattern seeded at `peak` at `charge` charge.

        Generates a theoretical isotopic pattern using :attr:`averagine`, calls
        :meth:`match_theoretical_isotopic_distribution`
        to extract experimental peaks matching this theoretical pattern, scales the theoretical distribution using
        :meth:`scale_theoretical_distribution`, and evaluates the quality of the fit using :attr:`scorer`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The putative monoisotopic peak to use for interpolating an isotopic pattern
        error_tolerance : float
            Parts-per-million error tolerance for isotopic pattern matching
        charge : int
            The charge state to produce an isotopic pattern for
        charge_carrier : float, optional
            The charge carrier mass, defaults to |PROTON|

        Returns
        -------
        :class:`~.IsotopicFitRecord`
            The fitted isotopic pattern
        """
        tid = self.averagine.isotopic_cluster(
            peak.mz, charge, charge_carrier=charge_carrier,
            truncate_after=truncate_after, ignore_below=ignore_below)
        eid = self.match_theoretical_isotopic_distribution(
            tid, error_tolerance=error_tolerance)
        record = self._evaluate_theoretical_distribution(eid, tid, peak, charge)
        return record

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON, truncate_after=0.8,
                              ignore_below=IGNORE_BELOW):
        """Given a set of candidate monoisotopic peaks and charge states, and a PPM error tolerance,
        fit each putative isotopic pattern.

        Calls :meth:`fit_theoretical_distribution` on each candidate.

        If a fit does not satisfy :attr:`scorer` `.reject`, it is discarded. If a fit has only one real peak
        and has a charge state greater than 1, it will also be discarded.

        Parameters
        ----------
        peak_charge_set : set
            The set of candidate (:class:`~.FittedPeak`, charge) tuples to try to fit
        error_tolerance : float
            Matching error tolerance
        charge_carrier : float, optional
            The charge carrier to use. Defaults to |PROTON|

        Returns
        -------
        set
            The set of :class:`~.IsotopicFitRecord` instances produced
        """
        results = []
        for peak, charge in peak_charge_set:
            if peak.mz < 1:
                continue
            fit = self.fit_theoretical_distribution(
                peak, error_tolerance, charge, charge_carrier, truncate_after,
                ignore_below)
            fit.missed_peaks = count_placeholders(fit.experimental)
            if not self._check_fit(fit):
                continue
            results.append(fit)
        return set(results)


try:
    from ms_deisotope._c.deconvoluter_base import DeconvoluterBase, AveragineDeconvoluterBase
except ImportError:
    pass

try:
    _quick_charge = quick_charge
    from ms_deisotope._c.deconvoluter_base import quick_charge
except ImportError:
    pass


def charge_range_(lo, hi, step=None):
    sign = -1 if lo < 0 else 1
    abs_lo, abs_hi = abs(lo), abs(hi)
    upper = max(abs_lo, abs_hi)
    lower = min(abs_lo, abs_hi)

    for c in range(upper, lower - 1, -1):
        yield c * sign


class ChargeIterator(object):
    def __init__(self, lo, hi):
        self.set_bounds(lo, hi)
        self.make_sequence()

    def set_bounds(self, lo, hi):
        self.sign = -1 if lo < 0 else 1
        abs_lo, abs_hi = abs(lo), abs(hi)
        if abs_lo < abs_hi:
            self.lower = abs_lo
            self.upper = abs_hi
        else:
            self.lower = abs_hi
            self.upper = abs_lo
        self.size = self.upper - self.lower + 1

    def make_sequence(self):
        self.index = 0
        self.size = self.upper - self.lower + 1
        self.values = [self.sign * (self.upper - i) for i in range(self.size)]

    def __len__(self):
        return self.size

    def reset(self):
        self.index = 0

    def sequence_from_quickcharge(self, peak_set, peak):
        charges = quick_charge(peak_set, peak.peak_count, abs(self.lower), abs(self.upper))
        n = charges.shape[0]
        self.index = 0
        if n == 0:
            self.size = 1
            self.values = [1 * self.sign]
        elif charges[0] != 1:
            self.size = n + 1
            self.values = [1 * self.sign] + list(self.sign * charges)
        else:
            self.size = n
            self.values = self.sign * charges

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == self.size:
            raise StopIteration()
        value = self.values[self.index]
        self.index += 1
        return value

    def next(self):
        return self.__next__()


class ExhaustivePeakSearchDeconvoluterBase(DeconvoluterBase):
    """Provides common methods for algorithms which attempt to find a deconvolution for every peak
    in a spectrum. This assumes no dependence between different peaks, instead it relies on subtraction,
    breadth of search, and order of encounter to avoid artefactual fits. This is usually not reasonable,
    so instead please use this class's extension, :class:`PeakDependenceGraphDeconvoluterBase` which can
    express dependence of fits on common resources.

    This class is not meant to be instantiated, but instead used as a mixin for classes that also
    inherit from :class:`DeconvoluterBase` and provide methods `fit_theoretical_distribution`
    and `_fit_peaks_at_charges`

    """
    def __init__(self, peaklist, *args, **kwargs):
        # Don't call superclass constructor from mixin?
        super(ExhaustivePeakSearchDeconvoluterBase, self).__init__(peaklist, *args, **kwargs)
        self.use_quick_charge = kwargs.get("use_quick_charge", False)

    def _get_all_peak_charge_pairs(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                   left_search_limit=3, right_search_limit=3,
                                   recalculate_starting_peak=True, use_quick_charge=False):
        """Construct the set of all unique candidate (monoisotopic peak, charge state) pairs using
        the provided search parameters.

        The search is performed using :func:`has_previous_peak_at_charge`, :func:`has_successor_peak_at_charge`,
        :meth:`_find_previous_putative_peak`, and :meth:`_find_next_putative_peak`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
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

        target_peaks = set()

        charges = ChargeIterator(*charge_range)

        if use_quick_charge:
            charges.sequence_from_quickcharge(self.peaklist, peak)

        if self.verbose:
            info("Considering charge range %r for %r" %
                 (list(charge_range_(*charge_range)), peak))
        for charge in charges:
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

    def _fit_all_charge_states(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), left_search_limit=3,
                               right_search_limit=3, recalculate_starting_peak=True, charge_carrier=PROTON,
                               truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Carry out the fitting process for `peak`.

        This method calls :meth:`_get_all_peak_charge_pairs` to collect all hypothetical solutions
        for `peak`, and invokes :meth:`_fit_peaks_at_charges` to evaluate them.

        The method :meth:`_fit_peaks_at_charges` is required by this interface, but is not defined by
        it, as it depends upon the underlying isotopic pattern fitting algorithm. See one of the
        Averagine-based algorithms for an implementation, such as :class:`AveragineDeconvoluterBase`,
        a complementary ancestor with this class.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        recalculate_starting_peak : bool, optional
            Whether or not to re-calculate the putative starting peak m/z based upon nearby
            peaks close to where isotopic peaks for `peak` should be. Defaults to True
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        set
            The set of :class:`~.IsotopicFitRecord` instances produced
        """
        target_peaks = self._get_all_peak_charge_pairs(
            peak, error_tolerance=error_tolerance,
            charge_range=charge_range,
            left_search_limit=left_search_limit,
            right_search_limit=right_search_limit,
            recalculate_starting_peak=recalculate_starting_peak,
            use_quick_charge=self.use_quick_charge)

        results = self._fit_peaks_at_charges(
            target_peaks, error_tolerance, charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)
        return (results)

    def charge_state_determination(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                   left_search_limit=3, right_search_limit=3,
                                   charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                                   ignore_below=IGNORE_BELOW):
        """Determine the optimal isotopic fit for `peak`, extracting it's charge state and monoisotopic peak.

        This method invokes :meth:`_fit_all_charge_states`, and then uses :attr:`scorer`'s `select` method to
        choose the optimal solution.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.IsotopicFitRecord`
            The best scoring isotopic fit
        """
        results = self._fit_all_charge_states(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)

        if self.verbose:
            info("Fits for %r" % peak)
            for rec in sorted(results)[-10:]:
                info(rec)
        try:
            result = self.scorer.select(results)
            return result
        except ValueError:
            return None

    def _make_deconvoluted_peak(self, fit, charge_carrier):
        score, charge, eid, tid = fit
        rep_eid = drop_placeholders(eid)
        total_abundance = sum(p.intensity for p in rep_eid)
        monoisotopic_mass = neutral_mass(
            tid.monoisotopic_mz, charge, charge_carrier)
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
            envelope=[(p.mz, p.intensity) for p in eid],
            mz=tid.monoisotopic_mz, fit=fit,
            area=sum(e.area for e in eid))
        return dpeak

    def deconvolute_peak(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                         left_search_limit=3, right_search_limit=3, charge_carrier=PROTON,
                         truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Perform a deconvolution for `peak`, generating a new :class:`ms_deisotope.peak_set.DeconvolutedPeak` instance
        corresponding to the optimal solution.

        This new peak has an m/z matching the monoisotopic peak of the pattern containing `peak`, and its intensity
        is the sum of all the matched peaks in its isotopic pattern. Its charge, isotopic fit, and other qualities
        are derived from the :class:`ms_deisotope.scoring.IsotopicFitRecord` instance corresponding to its best
        solution.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.DeconvolutedPeak`
        """
        fit = self.charge_state_determination(
            peak, error_tolerance=error_tolerance,
            charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)
        if fit is None:
            return
        tid = fit.theoretical
        dpeak = self._make_deconvoluted_peak(fit, charge_carrier)
        self._deconvoluted_peaks.append(dpeak)
        if self.use_subtraction:
            self.subtraction(tid, error_tolerance)
        return dpeak

    def targeted_deconvolution(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                               left_search_limit=3, right_search_limit=3,
                               charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                               ignore_below=IGNORE_BELOW):
        """Express the intent that this peak's deconvolution solution will be retrieved at a later point in the process
        and that it should be deconvoluted, and return a handle to retrieve the results with.

        This algorithm's implementation is simple enough that this is equivalent to just performing the deconvolution
        now and storing the result in a :class:`~.TrivialTargetedDeconvolutionResult` instance.

        Otherwise identical to :meth:`deconvolute_peak`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.TrivialTargetedDeconvolutionResult`
        """
        dpeak = self.deconvolute_peak(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)
        result = TrivialTargetedDeconvolutionResult(self, dpeak, peak)
        return result

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                    order_chooser=operator.attrgetter('index'),
                    left_search_limit=3, right_search_limit=3, charge_carrier=PROTON,
                    truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Completely deconvolute the spectrum.

        Visit each peak in the order chosen by `order_chooser`, and call :meth:`deconvolute_peak`
        on it with the provided arguments. This assumes all overlaps in isotopic pattern are captured
        by the search limits. This is usually not the case. For an alternative see
        :class:`PeakDependenceGraphDeconvoluterBase`

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        order_chooser : callable, optional:
            A callable used as a key function for sorting peaks into the order they will
            be visited during deconvolution. Defaults to :obj:`operator.attrgetter("index")`
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`
        """
        i = 0
        for peak in sorted(self.peaklist, key=order_chooser, reverse=True):
            if peak.mz < 2 or peak.intensity < self.minimum_intensity:
                continue
            self.deconvolute_peak(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit,
                right_search_limit=right_search_limit, charge_carrier=charge_carrier,
                truncate_after=truncate_after, ignore_below=ignore_below)
            i += 1

        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(
                self._deconvoluted_peaks)

        return DeconvolutedPeakSet(self._deconvoluted_peaks).reindex()


try:
    from ms_deisotope._c.deconvoluter_base import (
        _get_all_peak_charge_pairs as _c_get_all_peak_charge_pairs,
        _make_deconvoluted_peak as _c_make_deconvoluted_peak)
    ExhaustivePeakSearchDeconvoluterBase._get_all_peak_charge_pairs = _c_get_all_peak_charge_pairs
    ExhaustivePeakSearchDeconvoluterBase._make_deconvoluted_peak = _c_make_deconvoluted_peak
except ImportError as e:
    pass


class AveragineDeconvoluter(AveragineDeconvoluterBase, ExhaustivePeakSearchDeconvoluterBase):
    """A Deconvoluter which uses an :title-reference:`averagine` [1] model to generate theoretical
    isotopic patterns for each peak to consider. Combines :class:`AveragineDeconvoluterBase` and
    :class:`ExhaustivePeakSearchDeconvoluterBase` to create a working Deconvoluter type.


    Attributes
    ----------
    averagine : :class:`~.AveragineCache`
        The averagine model and associated theoretical isotopic pattern cache to use
        to build theoretical isotopic patterns.
    peaklist : :class:`~.PeakSet`
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
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
        self.peaklist = prepare_peaklist(peaklist)
        self.averagine = averagine
        self.scorer = scorer
        self._deconvoluted_peaks = []
        self.verbose = verbose

        super(AveragineDeconvoluter, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks=True, **kwargs)
        # ExhaustivePeakSearchDeconvoluterBase.__init__(self, peaklist, **kwargs)

    def config(self):
        return {
            "scale_method": self.scale_method,
            "use_subtraction": self.use_subtraction,
            "verbose": self.verbose,
            "scorer": self.scorer,
            "averagine": self.averagine
        }


class MultiAveragineDeconvoluterBase(DeconvoluterBase):

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON,
                              truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        results = []
        for peak, charge in peak_charge_set:
            for averagine in self.averagines:
                if peak.mz < 1:
                    continue
                fit = self.fit_theoretical_distribution(
                    peak, error_tolerance, charge, averagine, charge_carrier=charge_carrier,
                    truncate_after=truncate_after, ignore_below=ignore_below)
                fit.missed_peaks = count_placeholders(fit.experimental)
                fit.data = averagine
                if not self._check_fit(fit):
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
    averagine : list of :class:`~.ms_deisotope.averagine.AveragineCache`
        The averagine models and associated theoretical isotopic pattern caches to use
        to build theoretical isotopic patterns.
    peaklist : :class:`~.ms_peak_picker.PeakSet`
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : :class:`~.ms_deisotope.scoring.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    verbose : bool
        How much diagnostic information to provide

    References
    ----------
    [1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
        for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
        Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8
    """
    def __init__(self, peaklist, averagines=None, scorer=penalized_msdeconv,
                 use_subtraction=True, scale_method='sum',
                 merge_isobaric_peaks=True, minimum_intensity=5.,
                 verbose=False, *args, **kwargs):
        self.peaklist = prepare_peaklist(peaklist)
        self.scorer = scorer
        self.use_subtraction = use_subtraction
        self.scale_method = scale_method

        cache_backend = dict
        if averagines is None:
            averagines = [peptide, glycopeptide, glycan]
        averagines = [
            AveragineCache(avg, backend=cache_backend()) if not isinstance(
                avg, AveragineCache) else avg
            for avg in averagines]
        self.averagines = averagines
        self.verbose = verbose

        self._deconvoluted_peaks = []

        super(MultiAveragineDeconvoluter, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)
        # ExhaustivePeakSearchDeconvoluterBase.__init__(self, peaklist, **kwargs)


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
    peak_dependency_network : :class:`~.PeakDependenceGraph`
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
    """
    def __init__(self, peaklist, *args, **kwargs):
        max_missed_peaks = kwargs.get("max_missed_peaks", 1)
        self.subgraph_solver_type = kwargs.get("subgraph_solver", 'disjoint')
        super(PeakDependenceGraphDeconvoluterBase, self).__init__(peaklist, *args, **kwargs)
        self.peak_dependency_network = PeakDependenceGraph(
            self.peaklist, maximize=self.scorer.is_maximizing())
        self.max_missed_peaks = max_missed_peaks
        self.fit_postprocessor = kwargs.get("fit_postprocessor", None)
        self._priority_map = {}

    @property
    def max_missed_peaks(self):
        return self.peak_dependency_network.max_missed_peaks

    @max_missed_peaks.setter
    def max_missed_peaks(self, value):
        self.peak_dependency_network.max_missed_peaks = value

    def _explore_local(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                       ignore_below=IGNORE_BELOW):
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
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 1
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 0
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
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
            charge_carrier=charge_carrier, truncate_after=truncate_after, ignore_below=ignore_below)

        hold = set()
        for fit in results:
            if fit.charge > 1 and len(drop_placeholders(fit.experimental)) == 1:
                continue
            hold.add(fit)

        results = hold

        n = len(results)
        stop = max(min(n // 2, 100), 10)
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

    def populate_graph(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, charge_carrier=PROTON,
                       truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Visit each experimental peak and execute :meth:`_explore_local` on it with the provided
        parameters, populating the peak dependence graph with all viable candidates.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 1
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 0
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        """
        for peak in self.peaklist:
            if peak in self._priority_map or peak.intensity < self.minimum_intensity:
                continue
            self._explore_local(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                charge_carrier=charge_carrier,
                truncate_after=truncate_after, ignore_below=ignore_below)

    def postprocess_fits(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                         charge_carrier=PROTON, *args, **kwargs):
        if self.fit_postprocessor is None:
            return

    def select_best_disjoint_subgraphs(self, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
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
        if self.subgraph_solver_type == 'disjoint':
            solver = self._solve_subgraph_disjoint
        elif self.subgraph_solver_type == 'iterative':
            solver = self._solve_subgraph_iterative
        elif self.subgraph_solver_type == 'top':
            solver = self._solve_subgraph_top
        else:
            raise ValueError("Unknown solver type %r" % (self.subgraph_solver_type, ))

        for cluster in disjoint_envelopes:
            solutions = solver(cluster, error_tolerance, charge_carrier)
            for dpeak in solutions:
                self.peak_dependency_network.add_solution(dpeak.fit, dpeak)
                self._deconvoluted_peaks.append(dpeak)
                i += 1

    def _solve_subgraph_top(self, cluster, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Given a :class:`~.DependenceCluster`, return the single best fit from the collection of
        co-dependent fits.

        Parameters
        ----------
        cluster : :class:`~.DependenceCluster`
            The connected subgraph whose nodes will be searched
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result

        Returns
        -------
        list of :class:`~DeconvolutedPeak`
            The solved deconvolution solutions
        """
        fit = cluster[0]
        score, charge, eid, tid = fit
        rep_eid = drop_placeholders(eid)
        if len(rep_eid) == 0:
            return []
        dpeak = self._make_deconvoluted_peak(fit, charge_carrier)
        if self.use_subtraction:
            self.subtraction(tid, error_tolerance)
        return [dpeak]

    def _solve_subgraph_disjoint(self, cluster, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Given a :class:`~.DependenceCluster`, find a greedy disjoint set of isotopic fits.

        Parameters
        ----------
        cluster : :class:`~.DependenceCluster`
            The connected subgraph whose nodes will be searched
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result

        Returns
        -------
        list of :class:`~DeconvolutedPeak`
            The solved deconvolution solutions

        """
        disjoint_best_fits = cluster.disjoint_best_fits()
        i = 0
        solutions = []
        for fit in disjoint_best_fits:
            score, charge, eid, tid = fit
            rep_eid = drop_placeholders(eid)
            if len(rep_eid) == 0:
                continue
            dpeak = self._make_deconvoluted_peak(fit, charge_carrier)
            solutions.append(dpeak)
            i += 1
            if self.use_subtraction:
                self.subtraction(tid, error_tolerance)
        return solutions

    def _solve_subgraph_iterative(self, cluster, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Given a :class:`~.DependenceCluster`, build a :class:`~.ConnectedSubgraph` and incrementally
        subtract the best fitting solution and update its overlapping envelopes.

        Parameters
        ----------
        cluster : :class:`~.DependencyCluster`
            The connected subgraph whose nodes will be searched
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result

        Returns
        -------
        list of :class:`~DeconvolutedPeak`
            The solved deconvolution solutions
        """
        subgraph = cluster.build_graph()
        solutions = []
        mask = set()
        best_node = subgraph[0]
        peak = self._make_deconvoluted_peak(best_node.fit, charge_carrier)
        solutions.append(peak)
        if self.use_subtraction:
            self.subtraction(best_node.fit.theoretical, error_tolerance)
        mask.add(best_node.index)

        overlapped_nodes = list(best_node.overlap_edges)
        maximize = subgraph.maximize

        n = len(subgraph)
        while len(mask) != n:
            best_node = None
            best_score = 0 if maximize else float('inf')
            retained = []
            for node in overlapped_nodes:
                if node.index in mask:
                    continue
                missed_peaks = node.fit.missed_peaks = count_placeholders(node.fit.experimental)
                total_peaks = len(node.fit.experimental)
                invalid_peak_count = (missed_peaks >= total_peaks - 1 and abs(node.fit.charge) > 1
                                      ) or missed_peaks == total_peaks
                if invalid_peak_count or missed_peaks > self.max_missed_peaks:
                    mask.add(node.index)
                    continue
                fit = node.fit
                fit.theoretical.normalize().scale(fit.experimental, self.scale_method)
                fit.score = self.scorer.evaluate(self.peaklist, fit.experimental, fit.theoretical.peaklist)
                if self.scorer.reject(fit):
                    mask.add(node.index)
                    continue
                else:
                    retained.append(node)
                if maximize:
                    if fit.score > best_score:
                        best_node = node
                        best_score = fit.score
                else:
                    if fit.score < best_score:
                        best_node = node
                        best_score = fit.score

            if best_node is not None:
                peak = self._make_deconvoluted_peak(best_node.fit, charge_carrier)
                solutions.append(peak)
                if self.use_subtraction:
                    self.subtraction(best_node.fit.theoretical, error_tolerance)
                mask.add(best_node.index)
                overlapped_nodes = [node for node in best_node.overlap_edges if node.index not in mask]
            else:
                overlapped_nodes = []
            if not overlapped_nodes and len(mask) != n:
                overlapped_nodes = [node for node in subgraph if node.index not in mask]
        return solutions

    def targeted_deconvolution(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                               left_search_limit=3, right_search_limit=3, charge_carrier=PROTON,
                               truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Express the intent that this peak's deconvolution solution will be retrieved at a later point in the process
        and that it should be deconvoluted, and return a handle to retrieve the results with.

        As the operation does not immediately result in a deconvoluted peak but just adds the resulting fits to
        :attr:`peak_dependency_network`, this method constructs an instance of
        :class:`~.NetworkedTargetedDeconvolutionResult` which holds all
        the required information for recovering the best fit containing `peak`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|

        Returns
        -------
        :class:`~.NetworkedTargetedDeconvolutionResult`
        """
        self._explore_local(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit,
            right_search_limit=right_search_limit, charge_carrier=charge_carrier,
            truncate_after=truncate_after, ignore_below=ignore_below)
        result = NetworkedTargetedDeconvolutionResult(self, peak)
        self._priority_map[peak] = result
        return result

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                    left_search_limit=1, right_search_limit=0, iterations=MAX_ITERATION,
                    charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW,
                    convergence=CONVERGENCE):
        """Completely deconvolute the spectrum.

        For each iteration, clear :attr:`peak_depencency_network`, then invoke :meth:`populate_graph`
        followed by :meth:`select_best_disjoint_subgraphs` to populate the resulting
        :class:`~.DeconvolutedPeakSet`

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        order_chooser : callable, optional:
            A callable used as a key function for sorting peaks into the order they will
            be visited during deconvolution. Defaults to :obj:`operator.attrgetter("index")`
        left_search_limit : int, optional
            The number of steps to search to the left of :obj:`peak`. Defaults to 1
        right_search_limit : int, optional
            The number of steps to search to the right of :obj:`peak`. Defaults to 0
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        ignore_below : float, optional
            The minimum relative abundance to consider a peak in a theoretical isotopic
            pattern
        convergence : float, optional
            The threshold of the below which after the `(sum(intensity_before) - sum(
            intensity_after)) / sum(intensity_after)`

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`
        """
        if not self.use_subtraction:
            iterations = 1

        begin_signal = sum([p.intensity for p in self.peaklist])
        for i in range(iterations):
            self.peak_dependency_network.reset()
            self.populate_graph(
                error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                charge_carrier=charge_carrier,
                truncate_after=truncate_after, ignore_below=ignore_below)
            self.postprocess_fits(
                charge_range=charge_range,
                charge_carrier=charge_carrier,
                error_tolerance=error_tolerance)
            self.select_best_disjoint_subgraphs(error_tolerance, charge_carrier)
            self._slice_cache.clear()
            end_signal = sum([p.intensity for p in self.peaklist]) + 1

            if (begin_signal - end_signal) / end_signal < convergence:
                break
            begin_signal = end_signal

        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(
                self._deconvoluted_peaks)

        return DeconvolutedPeakSet(list(self._deconvoluted_peaks)).reindex()


class AveraginePeakDependenceGraphDeconvoluter(AveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    """A Deconvoluter which uses an :title-reference:`averagine` [1] model to generate theoretical
    isotopic patterns for each peak to consider, using a peak dependence graph to solve complex mass
    spectra.

    Extends :class:`AveragineDeconvoluter` to include features from
    :class:`PeakDependenceGraphDeconvoluterBase` making it suitable for deconvoluting complex spectra where
    peak overlaps are common.

    Attributes
    ----------
    peaklist : :class:`~.PeakSet`
        The centroided mass spectrum to deconvolute
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    averagine : :class:`~.AveragineCache`
        The averagine model and associated theoretical isotopic pattern cache to use
        to build theoretical isotopic patterns.
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : :class:`~.PeakDependenceGraph`
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
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

    References
    ----------
    [1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
        for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
        Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8
    """
    def __init__(self, peaklist, *args, **kwargs):
        super(AveraginePeakDependenceGraphDeconvoluter, self).__init__(peaklist, *args, **kwargs)
        # AveragineDeconvoluter.__init__(self, peaklist, *args, **kwargs)
        # PeakDependenceGraphDeconvoluterBase.__init__(self, peaklist, **kwargs)


class MultiAveraginePeakDependenceGraphDeconvoluter(MultiAveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    """Extends :class:`MultiAveragineDeconvoluter` to include features from
    :class:`PeakDependenceGraphDeconvoluterBase` making it suitable for deconvoluting complex spectra where
    peak overlaps are common.

    Attributes
    ----------
    peaklist : :class:`~.ms_peak_picker.PeakSet`
        The centroided mass spectrum to deconvolute
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    averagine : list of :class:`~.ms_deisotope.averagine.AveragineCache`
        The averagine model and associated theoretical isotopic pattern cache to use
        to build theoretical isotopic patterns.
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : :class:`~.PeakDependenceGraph`
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
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
    def __init__(self, peaklist, *args, **kwargs):
        super(MultiAveraginePeakDependenceGraphDeconvoluter, self).__init__(peaklist, *args, **kwargs)
        # MultiAveragineDeconvoluter.__init__(self, peaklist, *args, **kwargs)
        # PeakDependenceGraphDeconvoluterBase.__init__(self, peaklist, **kwargs)


class CompositionListDeconvoluterBase(DeconvoluterBase):
    """A mixin class to provide common features for deconvoluters which process spectra
    using a list of targeted compositions.

    Attributes
    ----------
    composition_list : list of :class:`~.Mapping`
        A series of objects which represent elemental compositions and support
        the :class:`~.Mapping` interface to access their individual elements.
    """
    def __init__(self, composition_list, *args, **kwargs):
        self.composition_list = list(composition_list)
        super(CompositionListDeconvoluterBase, self).__init__(*args, **kwargs)

    def generate_theoretical_isotopic_cluster(self, composition, charge, truncate_after=TRUNCATE_AFTER,
                                              mass_shift=None, charge_carrier=PROTON,
                                              ignore_below=IGNORE_BELOW):
        """Generate a theoretical isotopic pattern for ``composition``

        Parameters
        ----------
        composition : :class:`~.Mapping`
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
            each incremental change in charge state. Defaults to |PROTON|

        Returns
        -------
        :class:`~.TheoreticalIsotopicPattern`
            The theoretical isotopic pattern generated
        """
        tid = isotopic_variants(composition, charge=charge, charge_carrier=charge_carrier)
        tid = TheoreticalIsotopicPattern(tid, tid[0].mz)
        tid.truncate_after(truncate_after)
        tid.ignore_below(ignore_below)
        if mass_shift is not None:
            tid.shift(mass_shift / abs(charge))
        return tid

    def recalibrate_theoretical_mz(self, theoretical_distribution, experimental_mz):
        theoretical_distribution.shift(experimental_mz)
        return theoretical_distribution

    def fit_composition_at_charge(self, composition, charge, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON,
                                  truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW, mass_shift=None):
        """Produce an isotopic fit for `composition` at `charge` against the experimental peak set.

        This method requires that the instance also possess a method named `match_theoretical_isotopic_distribution`
        such as the one implemented in :class:`DeconvoluterBase`.

        Parameters
        ----------
        composition : :class:`~.Mapping`
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
            each incremental change in charge state. Defaults to |PROTON|

        Returns
        -------
        :class:`~.IsotopicFitRecord`
        """
        tid = self.generate_theoretical_isotopic_cluster(composition, charge=charge, truncate_after=truncate_after,
                                                         charge_carrier=charge_carrier, ignore_below=ignore_below,
                                                         mass_shift=mass_shift)
        monoisotopic_peak = self.peaklist.has_peak(tid[0].mz, error_tolerance)
        if monoisotopic_peak is not None:
            tid = self.recalibrate_theoretical_mz(tid, monoisotopic_peak.mz)
        eid = self.match_theoretical_isotopic_distribution(
            tid.peaklist, error_tolerance)

        missed_peaks = count_placeholders(eid)

        if missed_peaks > len(eid) / 2:
            return None

        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer.evaluate(self.peaklist, eid, tid.peaklist)
        fit = IsotopicFitRecord(None, score, charge, tid, eid)
        fit.missed_peaks = missed_peaks
        return fit

    def _make_deconvoluted_peak_solution(self, fit, composition, charge_carrier):
        eid = fit.experimental
        tid = fit.theoretical
        charge = fit.charge
        rep_eid = drop_placeholders(eid)
        total_abundance = sum(
            p.intensity for p in eid if p.intensity > 1)

        monoisotopic_mass = neutral_mass(
            tid.monoisotopic_mz, charge, charge_carrier)
        monoisotopic_mz = tid.monoisotopic_mz

        reference_peak = first_peak(eid)
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
        return peak

    def deconvolute_composition(self, composition, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW,
                                mass_shift=None):
        """For each charge state under consideration, fit the theoretical isotopic pattern for this composition,
        and if the fit is satisfactory, add it to the results set.

        Parameters
        ----------
        composition : :class:`~.Mapping`
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
            each incremental change in charge state. Defaults to |PROTON|
        """
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(composition, charge=charge, error_tolerance=error_tolerance,
                                                 truncate_after=truncate_after, charge_carrier=charge_carrier,
                                                 mass_shift=mass_shift, ignore_below=ignore_below)
            if fit is None:
                continue
            if not self.scorer.reject(fit):
                eid = fit.experimental
                tid = fit.theoretical
                rep_eid = drop_placeholders(eid)
                if (len(rep_eid) < 2) or (len(rep_eid) < (len(tid) / 2.)) or (len(rep_eid) == 1 and fit.charge > 1):
                    continue

                peak = self._make_deconvoluted_peak_solution(
                    fit, composition, charge_carrier)
                self._deconvoluted_peaks.append(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)


class CompositionListDeconvoluter(CompositionListDeconvoluterBase):
    '''
    Attributes
    ----------
    composition_list : list of :class:`~.Mapping`
        A series of objects which represent elemental compositions and support
        the :class:`~.Mapping` interface to access their individual elements.
    peaklist : :class:`~ms_peak_picker.PeakSet`
        The collection of :class:`~.ms_peak_picker.FittedPeak` instances and possible associated
        data to deconvolute.
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
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
    '''
    def __init__(self, peaklist, composition_list, scorer,
                 use_subtraction=False, scale_method='sum',
                 verbose=False):
        self.peaklist = prepare_peaklist(peaklist)
        self.scorer = scorer
        self.verbose = verbose
        self._deconvoluted_peaks = []
        super(CompositionListDeconvoluter, self).__init__(
            composition_list,
            use_subtraction=use_subtraction, scale_method=scale_method, merge_isobaric_peaks=True)

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), charge_carrier=PROTON,
                    truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW, mass_shift=None, **kwargs):
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance=error_tolerance,
                                         charge_range=charge_range, charge_carrier=charge_carrier,
                                         truncate_after=truncate_after, ignore_below=ignore_below,
                                         mass_shift=mass_shift)
        return DeconvolutedPeakSet(self._deconvoluted_peaks).reindex()


class CompositionListPeakDependenceGraphDeconvoluter(CompositionListDeconvoluter):
    '''
    Attributes
    ----------
    composition_list : list of :class:`~.Mapping`
        A series of objects which represent elemental compositions and support
        the :class:`~.Mapping` interface to access their individual elements.
    peaklist : :class:`~ms_peak_picker.PeakSet`
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : :class:`~PeakDependenceGraph`
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
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
    '''
    def __init__(self, peaklist, composition_list, scorer,
                 use_subtraction=False, scale_method='sum',
                 verbose=False, **kwargs):
        max_missed_peaks = kwargs.get("max_missed_peaks", 1)
        super(CompositionListPeakDependenceGraphDeconvoluter, self).__init__(
            peaklist, composition_list, scorer=scorer, use_subtraction=use_subtraction,
            scale_method=scale_method,
            verbose=verbose, **kwargs)

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

    def deconvolute_composition(self, composition, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW,
                                mass_shift=None):
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(
                composition, charge, error_tolerance, charge_carrier=charge_carrier,
                truncate_after=truncate_after, mass_shift=mass_shift, ignore_below=ignore_below)
            if fit is None:
                continue
            rep_eid = drop_placeholders(fit.experimental)
            if len(rep_eid) == 1 and fit.charge > 1:
                continue
            if not self.scorer.reject(fit):
                self.peak_dependency_network.add_fit_dependence(fit)

    def populate_graph(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), truncate_after=TRUNCATE_AFTER,
                       charge_carrier=PROTON, ignore_below=IGNORE_BELOW, mass_shift=None):
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance, charge_range,
                                         truncate_after=truncate_after, charge_carrier=charge_carrier,
                                         ignore_below=ignore_below, mass_shift=mass_shift)

    def select_best_disjoint_subgraphs(self, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        disjoint_envelopes = self.peak_dependency_network.find_non_overlapping_intervals()

        for cluster in disjoint_envelopes:
            for fit in cluster.disjoint_best_fits():
                eid = fit.experimental
                tid = fit.theoretical
                composition = fit.data
                rep_eid = drop_placeholders(eid)
                if len(rep_eid) < 2 or len(rep_eid) < len(tid) / 2.:
                    continue

                peak = self._make_deconvoluted_peak_solution(
                    fit, composition, charge_carrier)
                self._save_peak_solution(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), iterations=MAX_ITERATION,
                    truncate_after=TRUNCATE_AFTER, charge_carrier=PROTON, ignore_below=IGNORE_BELOW,
                    mass_shift=None, convergence=CONVERGENCE, **kwargs):
        if not self.use_subtraction:
            iterations = 1
        begin_signal = sum([p.intensity for p in self.peaklist])
        for i in range(iterations):
            self.populate_graph(error_tolerance, charge_range, charge_carrier=charge_carrier,
                                truncate_after=truncate_after, ignore_below=ignore_below,
                                mass_shift=mass_shift)
            self.select_best_disjoint_subgraphs(error_tolerance)
            self._slice_cache.clear()
            end_signal = sum([p.intensity for p in self.peaklist]) + 1
            if (begin_signal - end_signal) / end_signal < convergence:
                break
            begin_signal = end_signal
        return DeconvolutedPeakSet(self._deconvoluted_peaks).reindex()


_APDGD = AveraginePeakDependenceGraphDeconvoluter


class HybridAveragineCompositionListPeakDependenceGraphDeconvoluter(_APDGD, CompositionListDeconvoluterBase):
    def __init__(self, peaklist, composition_list, *args, **kwargs):
        super(HybridAveragineCompositionListPeakDependenceGraphDeconvoluter, self).__init__(
            peaklist, composition_list=composition_list, *args, **kwargs)

    def deconvolute_composition(self, composition, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                truncate_after=TRUNCATE_AFTER, charge_carrier=PROTON, ignore_below=IGNORE_BELOW,
                                mass_shift=None):
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(
                composition, charge, error_tolerance, charge_carrier=charge_carrier,
                truncate_after=truncate_after, mass_shift=mass_shift, ignore_below=ignore_below)
            if fit is None:
                continue
            rep_eid = drop_placeholders(fit.experimental)
            if len(rep_eid) == 1 and fit.charge > 1:
                continue
            if not self.scorer.reject(fit):
                fit.data = (composition, charge)
                self.peak_dependency_network.add_fit_dependence(fit)

    def populate_graph(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                       ignore_below=IGNORE_BELOW, mass_shift=None):
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
            charge_carrier=charge_carrier, truncate_after=truncate_after)


def deconvolute_peaks(peaklist, decon_config=None,
                      charge_range=(1, 8), error_tolerance=ERROR_TOLERANCE,
                      priority_list=None, left_search_limit=3, right_search_limit=3,
                      left_search_limit_for_priorities=None, right_search_limit_for_priorities=None,
                      verbose_priorities=False, verbose=False, charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                      deconvoluter_type=AveraginePeakDependenceGraphDeconvoluter, **kwargs):
    """Deconvolute a centroided mass spectrum

    This function constructs a deconvoluter object using the ``deconvoluter_type`` argument
    and deconvolutes the input ``peaklist`` by calling its :meth:`deconvolute` method.

    If ``priority_list`` is not :const:`None`, it is expected to be an iterable of either
    tuples of (:class:`~.FittedPeak`, ``(min charge, max charge)``) pairs, or instances of
    :class:`~.PriorityTarget`. These will be passed to :meth:`targeted_deconvolution` of
    the deconvoluter.

    Parameters
    ----------
    peaklist : :class:`~.PeakSet` or list of Peak-like objects
        The centroided mass spectrum to deconvolute.
    decon_config : dict, optional
        Parameters to use to initialize the deconvoluter instance produced by
        ``deconvoluter_type``
    charge_range : tuple of integers, optional
        The range of charge states to consider.
    error_tolerance : float, optional
        PPM error tolerance to use to match experimental to theoretical peaks
    priority_list : list, optional
        The set of peaks to target for deconvolution to be able to enforce external
        constraints on, such as selected precursors for fragmentation.
    left_search_limit : int, optional
        The maximum number of neutron shifts to search to the left  (decrease) from
        each query peak
    right_search_limit : int, optional
        The maximum number of neutron shifts to search to the right (increase) from
        each query peak
    left_search_limit_for_priorities : int, optional
        The maximum number of neutron shifts to search to the left (decrease) from
        each query peak for priority targets
    right_search_limit_for_priorities : int, optional
        The maximum number of neutron shifts to search to the right (increase) from
        each query peak for priority targets
    verbose_priorities : bool, optional
        Whether to turn on verbose mode for priority targets
    verbose : bool, optional
        Passed to the deconvoluter to enable verbose mode globally
    charge_carrier : float, optional
        The mass of the charge carrier. Defaults to |PROTON|
    truncate_after : float, optional
        The percentage of the isotopic pattern to include. Defaults to |TRUNCATE_AFTER|
    deconvoluter_type : type or callable, optional
        A callable returning a deconvoluter. Defaults to :class:`~.AveraginePeakDependenceGraphDeconvoluter`
    **kwargs
        Additional keywords included in ``decon_config``

    Returns
    -------
    :class:`~.DeconvolutionProcessResult`
    """
    if priority_list is None:
        priority_list = []
    if left_search_limit_for_priorities is None:
        left_search_limit_for_priorities = left_search_limit
    if right_search_limit_for_priorities is None:
        right_search_limit_for_priorities = right_search_limit

    decon_config = decon_config or {}
    decon_config.update(kwargs)
    decon_config.setdefault("use_subtraction", True)
    decon_config.setdefault("scale_method", SCALE_METHOD)
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
            p, error_tolerance=error_tolerance,
            charge_range=hinted_charge_range,
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
    errors = []
    for pr in priority_list_results:
        try:
            result = pr.get()
        except ValueError as e:
            result = None
            errors.append(e)
            logger.error("Could not extract a solution for %r", pr.query_peak, exc_info=True)
        acc.append(result)

    priority_list_results = acc

    return DeconvolutionProcessResult(
        decon, deconvoluted_peaks, priority_list_results, errors)
