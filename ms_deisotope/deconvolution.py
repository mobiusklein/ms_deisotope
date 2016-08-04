import operator

from .averagine import (
    AveragineCache, peptide, glycopeptide, glycan, neutral_mass, isotopic_variants,
    isotopic_shift, PROTON)
from .peak_set import DeconvolutedPeak, DeconvolutedPeakSolution, DeconvolutedPeakSet
from .scoring import IsotopicFitRecord, penalized_msdeconv, distinct_pattern_fitter
from .utils import range, Base
from .envelope_statistics import a_to_a2_ratio, average_mz, most_abundant_mz
from .peak_dependency_network import PeakDependenceGraph

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


def has_previous_peak_at_charge(peak_index, peak, charge=2, step=1):
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
    return peak_index.has_peak(prev)


def has_successor_peak_at_charge(peak_index, peak, charge=2, step=1):
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
    return peak_index.has_peak(nxt)


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
        experimental_distribution = [self.has_peak(
            p.mz, error_tolerance) for p in theoretical_distribution]
        return experimental_distribution

    def scale_theoretical_distribution(self, theoretical_distribution, experimental_distribution):
        if self.scale_method == 'sum':
            total_abundance = sum(p.intensity for p in experimental_distribution)
            for peak in theoretical_distribution:
                peak.intensity *= total_abundance
            return theoretical_distribution
        elif self.scale_method == 'max':
            i, peak = max(enumerate(theoretical_distribution), key=lambda x: x[1].intensity)
            scale_factor = experimental_distribution[i].intensity / peak.intensity
            for peak in theoretical_distribution:
                peak.intensity *= scale_factor
            return theoretical_distribution

    def subtraction(self, isotopic_cluster, error_tolerance=2e-5):
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
        shift = isotopic_shift(charge)
        prev_peak = mz - (shift)
        peaklist_slice = self.between(
            prev_peak - (prev_peak * tolerance),
            prev_peak + (prev_peak * tolerance))
        candidates = []
        for backward in peaklist_slice:
            prev_peak_mz = backward.mz
            if step == 1:
                candidates.extend(self._find_next_putative_peak(prev_peak_mz, charge, 1, tolerance))
            else:
                candidates.extend(
                    self._find_previous_putative_peak(prev_peak_mz, charge, step - 1, tolerance))
        return candidates

    def __repr__(self):
        type_name = self.__class__.__name__
        return "%s(peaklist=%s)" % (type_name, self.peaklist)


class AveragineDeconvoluterBase(DeconvoluterBase):
    def __init__(self, use_subtraction=False, scale_method="sum", merge_isobaric_peaks=True,
                 minimum_intensity=5., *args, **kwargs):
        super(AveragineDeconvoluterBase, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)

    def fit_theoretical_distribution(self, peak, error_tolerance, charge, charge_carrier=PROTON):
        tid = self.averagine.isotopic_cluster(peak.mz, charge, charge_carrier=charge_carrier)
        eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance=error_tolerance)
        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer(self.peaklist, eid, tid)
        return IsotopicFitRecord(peak, score, charge, tid, eid)

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON):
        results = []
        for peak, charge in peak_charge_set:
            if peak.mz < 1:
                continue
            fit = self.fit_theoretical_distribution(
                peak, error_tolerance, charge, charge_carrier)
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
    def _get_all_peak_charge_pairs(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=3,
                                   right_search_limit=3, use_charge_state_hint=False,
                                   recalculate_starting_peak=True):
        if use_charge_state_hint:
            upper_charge_limit = self.peaklist.predict_charge_state(peak) + 2
            if upper_charge_limit <= 1:
                upper_charge_limit = charge_range[-1]
            charge_range = (charge_range[0], min(upper_charge_limit, charge_range[-1]))

        target_peaks = set()
        if self.verbose:
            info("Considering charge range %r for %r" % (list(charge_range_(*charge_range)), peak))
        for charge in charge_range_(*charge_range):
            target_peaks.add((peak, charge))

            # Look Left
            for i in range(1, left_search_limit):
                prev_peak = has_previous_peak_at_charge(self.peaklist, peak, charge, i)
                if prev_peak is None:
                    continue
                target_peaks.add((prev_peak, charge))
                if recalculate_starting_peak:
                    target_peaks.update(self._find_previous_putative_peak(peak.mz, charge, i, 2 * error_tolerance))

            # Look Right
            for i in range(1, right_search_limit):
                nxt_peak = has_successor_peak_at_charge(self.peaklist, peak, charge, i)
                if nxt_peak is None:
                    continue
                target_peaks.add((nxt_peak, charge))

            if recalculate_starting_peak:
                for i in range(min(left_search_limit, 2)):
                    target_peaks.update(self._find_next_putative_peak(
                        peak.mz, charge, step=i, tolerance=2 * error_tolerance))

        return target_peaks

    def _fit_all_charge_states(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=3,
                               right_search_limit=3, use_charge_state_hint=False,
                               recalculate_starting_peak=True, charge_carrier=PROTON):
        target_peaks = self._get_all_peak_charge_pairs(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            use_charge_state_hint=use_charge_state_hint, recalculate_starting_peak=True)

        results = self._fit_peaks_at_charges(target_peaks, error_tolerance, charge_carrier=charge_carrier)
        return (results)

    def charge_state_determination(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=3,
                                   right_search_limit=3, use_charge_state_hint=False, charge_carrier=PROTON):
        results = self._fit_all_charge_states(
            peak, error_tolerance=error_tolerance, charge_range=charge_range, left_search_limit=left_search_limit,
            right_search_limit=right_search_limit, use_charge_state_hint=use_charge_state_hint,
            charge_carrier=charge_carrier)

        if self.verbose:
            info("\nFits for %r" % peak)
            for rec in sorted(results)[-10:]:
                info(rec)

        try:
            result = self.scorer.select(results)
            return result
        except ValueError:
            return None

    def deconvolute_peak(self, peak, error_tolerance=2e-5, charge_range=(1, 8), use_charge_state_hint=False,
                         left_search_limit=3, right_search_limit=3, charge_carrier=PROTON):
        charge_det = self.charge_state_determination(
            peak, error_tolerance=error_tolerance,
            charge_range=charge_range, use_charge_state_hint=use_charge_state_hint,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier)
        if charge_det is None:
            return
        score, charge, eid, tid = charge_det
        rep_eid = drop_placeholders(eid)
        total_abundance = sum(p.intensity for p in rep_eid)
        monoisotopic_mass = neutral_mass(eid[0].mz, charge, charge_carrier=charge_carrier)
        reference_peak = first_peak(eid)

        dpeak = DeconvolutedPeak(
            neutral_mass=monoisotopic_mass, intensity=total_abundance, charge=charge,
            signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
            index=reference_peak.index,
            full_width_at_half_max=mean(p.full_width_at_half_max for p in rep_eid),
            a_to_a2_ratio=a_to_a2_ratio(tid),
            most_abundant_mass=neutral_mass(most_abundant_mz(eid), charge),
            average_mass=neutral_mass(average_mz(eid), charge),
            score=score,
            envelope=[(p.mz, p.intensity) for p in rep_eid],
            mz=eid[0].mz)

        self._deconvoluted_peaks.append(dpeak)
        if self.use_subtraction:
            self.subtraction(tid, error_tolerance)
        return dpeak

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8),
                    order_chooser=operator.attrgetter('index'),
                    use_charge_state_hint=False, left_search_limit=3,
                    right_search_limit=3, charge_carrier=PROTON):
        i = 0
        for peak in sorted(self.peaklist, key=order_chooser, reverse=True):
            if peak.mz < 2 or peak.intensity < self.minimum_intensity:
                continue
            self.deconvolute_peak(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                use_charge_state_hint=use_charge_state_hint, left_search_limit=left_search_limit,
                right_search_limit=right_search_limit, charge_carrier=charge_carrier)
            i += 1

        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(self._deconvoluted_peaks)

        return DeconvolutedPeakSet(self._deconvoluted_peaks)._reindex()


class AveragineDeconvoluter(AveragineDeconvoluterBase, ExhaustivePeakSearchDeconvoluterBase):
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
    def fit_theoretical_distribution(self, peak, error_tolerance, charge, averagine, charge_carrier=PROTON):
        tid = averagine.isotopic_cluster(peak.mz, charge, charge_carrier=charge_carrier)
        eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance=error_tolerance)
        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer(self.peaklist, eid, tid)
        return IsotopicFitRecord(peak, score, charge, tid, eid)

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON):
        results = []
        for peak, charge in peak_charge_set:
            for averagine in self.averagine:
                if peak.mz < 1:
                    continue
                fit = self.fit_theoretical_distribution(
                    peak, error_tolerance, charge, averagine, charge_carrier=charge_carrier)
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

    def __init__(self, peaklist, *args, **kwargs):
        max_missed_peaks = kwargs.pop("max_missed_peaks", 1)
        super(PeakDependenceGraphDeconvoluterBase, self).__init__(peaklist=peaklist, *args, **kwargs)
        self.peak_dependency_network = PeakDependenceGraph(self.peaklist, maximize=self.scorer.is_maximizing())
        self.max_missed_peaks = max_missed_peaks

    @property
    def max_missed_peaks(self):
        return self.peak_dependency_network.max_missed_peaks

    @max_missed_peaks.setter
    def max_missed_peaks(self, value):
        self.peak_dependency_network.max_missed_peaks = value

    def _explore_local(self, peak, error_tolerance=2e-5, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, use_charge_state_hint=False, charge_carrier=PROTON):
        results = self._fit_all_charge_states(
            peak, error_tolerance=error_tolerance, charge_range=charge_range, left_search_limit=left_search_limit,
            right_search_limit=right_search_limit,
            use_charge_state_hint=use_charge_state_hint, charge_carrier=charge_carrier)

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
                       right_search_limit=0, use_charge_state_hint=False, charge_carrier=PROTON):
        for peak in self.peaklist:
            out = self._explore_local(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                use_charge_state_hint=use_charge_state_hint, charge_carrier=charge_carrier)

    def select_best_disjoint_subgraphs(self, error_tolerance=2e-5, charge_carrier=PROTON):
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
                monoisotopic_mass = neutral_mass(eid[0].mz, charge, charge_carrier)
                reference_peak = first_peak(eid)

                dpeak = DeconvolutedPeak(
                    neutral_mass=monoisotopic_mass, intensity=total_abundance, charge=charge,
                    signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
                    index=reference_peak.index,
                    full_width_at_half_max=mean(p.full_width_at_half_max for p in rep_eid),
                    a_to_a2_ratio=a_to_a2_ratio(tid),
                    most_abundant_mass=neutral_mass(most_abundant_mz(eid), charge),
                    average_mass=neutral_mass(average_mz(eid), charge),
                    score=score,
                    envelope=[(p.mz, p.intensity) for p in rep_eid],
                    mz=eid[0].mz)
                self._deconvoluted_peaks.append(dpeak)
                i += 1
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8),
                    use_charge_state_hint=False, left_search_limit=1,
                    right_search_limit=0, iterations=1, charge_carrier=PROTON):

        for i in range(iterations):
            self.populate_graph(
                error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                use_charge_state_hint=use_charge_state_hint, charge_carrier=charge_carrier)
            self.select_best_disjoint_subgraphs(error_tolerance, charge_carrier)

        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(self._deconvoluted_peaks)

        return DeconvolutedPeakSet(list(self._deconvoluted_peaks))._reindex()


class AveraginePeakDependenceGraphDeconvoluter(AveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    def __init__(self, peaklist, *args, **kwargs):
        AveragineDeconvoluter.__init__(self, peaklist, *args, **kwargs)
        PeakDependenceGraphDeconvoluterBase.__init__(self, peaklist, **kwargs)


class MultiAveraginePeakDependenceGraphDeconvoluter(MultiAveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    def __init__(self, peaklist, *args, **kwargs):
        MultiAveragineDeconvoluter.__init__(self, peaklist, *args, **kwargs)
        PeakDependenceGraphDeconvoluterBase.__init__(self, peaklist, **kwargs)


class CompositionListDeconvoluter(DeconvoluterBase):
    def __init__(self, peaklist, composition_list, scorer,
                 use_subtraction=False, scale_method='sum',
                 verbose=False):
        self.peaklist = peaklist.clone()
        self.composition_list = composition_list
        self.scorer = scorer
        self.verbose = verbose
        self._deconvoluted_peaks = []
        super(CompositionListDeconvoluter, self).__init__(
            use_subtraction=use_subtraction, scale_method=scale_method, merge_isobaric_peaks=True)

    def generate_theoretical_isotopic_cluster(self, composition, charge, truncate_after=0.95, charge_carrier=PROTON):
        cumsum = 0
        result = []
        for peak in isotopic_variants(composition, charge=charge, charge_carrier=charge_carrier):
            cumsum += peak.intensity
            result.append(peak)
            if cumsum >= truncate_after:
                break
        return result

    def fit_composition_at_charge(self, composition, charge, error_tolerance=2e-5, charge_carrier=PROTON,
                                  truncate_after=0.95):
        tid = self.generate_theoretical_isotopic_cluster(composition, charge=charge, truncate_after=truncate_after,
                                                         charge_carrier=charge_carrier)
        eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance)

        missed_peaks = count_placeholders(eid)

        if missed_peaks > len(eid) / 2:
            return None

        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer.evaluate(self.peaklist, eid, tid)
        fit = IsotopicFitRecord(None, score, charge, tid, eid)
        fit.missed_peaks = missed_peaks
        return fit

    def deconvolute_composition(self, composition, error_tolerance=2e-5, charge_range=(1, 8), charge_carrier=PROTON,
                                truncate_after=0.95):
        for charge in charge_range_(*charge_range):
            tid = self.generate_theoretical_isotopic_cluster(composition, charge=charge, truncate_after=truncate_after,
                                                             charge_carrier=charge_carrier)

            eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance)

            if count_placeholders(eid) > len(eid) / 2:
                continue

            self.scale_theoretical_distribution(tid, eid)
            score = self.scorer.evaluate(self.peaklist, eid, tid)
            fit = IsotopicFitRecord(None, score, charge, tid, eid)

            if not self.scorer.reject(fit):
                total_abundance = sum(p.intensity for p in eid if p.intensity > 1)

                monoisotopic_mass = neutral_mass(eid[0].mz, charge, charge_carrier)
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
                    full_width_at_half_max=mean(p.full_width_at_half_max for p in rep_eid),
                    a_to_a2_ratio=a_to_a2_ratio(tid),
                    most_abundant_mass=neutral_mass(most_abundant_mz(eid), charge),
                    average_mass=neutral_mass(average_mz(eid), charge),
                    score=score,
                    envelope=[(p.mz, p.intensity) for p in rep_eid],
                    mz=monoisotopic_mz)

                self._deconvoluted_peaks.append(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), charge_carrier=PROTON, truncate_after=0.95,
                    **kwargs):
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance=error_tolerance,
                                         charge_range=charge_range, charge_carrier=charge_carrier,
                                         truncate_after=truncate_after)
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

    def deconvolute_composition(self, composition, error_tolerance=2e-5, charge_range=(1, 8), truncate_after=0.95,
                                charge_carrier=PROTON):
        for charge in charge_range_(*charge_range):
            tid = self.generate_theoretical_isotopic_cluster(composition, charge=charge, charge_carrier=charge_carrier,
                                                             truncate_after=truncate_after)
            eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance)

            missed_peaks = count_placeholders(eid)

            if missed_peaks > len(eid) / 2 or len(eid) == 0:
                continue

            self.scale_theoretical_distribution(tid, eid)
            score = self.scorer.evaluate(self.peaklist, eid, tid)
            fit = IsotopicFitRecord(None, score, charge, tid, eid)
            fit.missed_peaks = missed_peaks
            fit.data = composition

            if not self.scorer.reject(fit):
                self.peak_dependency_network.add_fit_dependence(fit)

    def populate_graph(self, error_tolerance=2e-5, charge_range=(1, 8), truncate_after=0.95, charge_carrier=PROTON):
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance, charge_range,
                                         truncate_after=truncate_after, charge_carrier=charge_carrier)

    def select_best_disjoint_subgraphs(self, error_tolerance=2e-5, charge_carrier=PROTON):
        disjoint_envelopes = self.peak_dependency_network.find_non_overlapping_intervals()

        for cluster in disjoint_envelopes:
            for fit in cluster.disjoint_best_fits():
                eid = fit.experimental
                tid = fit.theoretical
                composition = fit.data
                charge = fit.charge
                score = fit.score

                total_abundance = sum(p.intensity for p in eid if p.intensity > 1)
                monoisotopic_mass = neutral_mass(eid[0].mz, charge, charge_carrier)
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
                    full_width_at_half_max=mean(p.full_width_at_half_max for p in rep_eid),
                    a_to_a2_ratio=a_to_a2_ratio(tid),
                    most_abundant_mass=neutral_mass(most_abundant_mz(eid), charge),
                    average_mass=neutral_mass(average_mz(eid), charge),
                    score=score,
                    envelope=[(p.mz, p.intensity) for p in rep_eid],
                    mz=monoisotopic_mz)

                self._save_peak_solution(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), iterations=1, truncate_after=0.95,
                    charge_carrier=PROTON, **kwargs):
        for i in range(iterations):
            self.populate_graph(error_tolerance, charge_range, charge_carrier=charge_carrier,
                                truncate_after=truncate_after)
            self.select_best_disjoint_subgraphs(error_tolerance)
        return DeconvolutedPeakSet(self._deconvoluted_peaks)._reindex()


def deconvolute_peaks(peaklist, deconvoluter_type=AveraginePeakDependenceGraphDeconvoluter,
                      decon_config=None, charge_range=(1, 8), error_tolerance=2e-5, priority_list=[],
                      use_charge_state_hint_for_priorities=True, left_search_limit=3, right_search_limit=3,
                      left_search_limit_for_priorities=None, right_search_limit_for_priorities=None,
                      verbose_priorities=False, verbose=False, charge_carrier=PROTON, **kwargs):
    if left_search_limit_for_priorities is None:
        left_search_limit_for_priorities = left_search_limit
    if right_search_limit_for_priorities is None:
        right_search_limit_for_priorities = right_search_limit

    decon_config = decon_config or {}
    decon_config.update(kwargs)
    decon = deconvoluter_type(peaklist=peaklist, **decon_config)

    if verbose_priorities or verbose:
        decon.verbose = True

    priority_list_results = []
    for p in priority_list:
        if not isinstance(p, FittedPeak):
            p = decon.peaklist.has_peak(p, error_tolerance)
        priority_result = decon.deconvolute_peak(
            p, error_tolerance=error_tolerance, charge_range=charge_range,
            use_charge_state_hint=use_charge_state_hint_for_priorities,
            left_search_limit=left_search_limit_for_priorities,
            right_search_limit=right_search_limit_for_priorities,
            charge_carrier=charge_carrier)
        priority_list_results.append(priority_result)

    if verbose_priorities and not verbose:
        decon.verbose = False

    deconvoluted_peaks = decon.deconvolute(
        error_tolerance=error_tolerance, charge_range=charge_range, left_search_limit=left_search_limit,
        right_search_limit=right_search_limit, charge_carrier=charge_carrier)
    return deconvoluted_peaks, priority_list_results
