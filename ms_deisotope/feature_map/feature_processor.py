from collections import defaultdict
from itertools import product

import numpy as np

from ms_peak_picker import FittedPeak

from .feature_map import (
    LCMSFeatureMap,
    DeconvolutedLCMSFeatureMap,
    smooth_overlaps_neutral)
from .lcms_feature import (
    LCMSFeature,
    EmptyFeature,
    FeatureSetIterator)
from .feature_fit import (
    LCMSFeatureSetFit,
    DeconvolutedLCMSFeatureTreeNode,
    DeconvolutedLCMSFeature)
from .dependence_network import FeatureDependenceGraph
from .profile_transform import binsearch, smooth_leveled
from ms_deisotope.peak_dependency_network.intervals import Interval, IntervalTreeNode
from ms_deisotope.averagine import (
    AveragineCache, PROTON, isotopic_shift,
    neutral_mass)
from ms_deisotope.peak_set import DeconvolutedPeak
from ms_deisotope.envelope_statistics import (
    most_abundant_mz, average_mz, a_to_a2_ratio)
from ms_deisotope.deconvolution import (
    charge_range_, drop_placeholders, first_peak,
    mean)
from ms_deisotope.task import LogUtilsMixin


def conform_envelopes(experimental, base_theoretical, minimum_theoretical_abundance=0.05):
    total = 0
    n_missing = 0
    i = 0
    cleaned_eid = []
    for peak in experimental:
        if peak is None:
            peak = FittedPeak(base_theoretical[i].mz, 1, 1, -1, -1, 0, 1, 0, 0)
            if base_theoretical[i].intensity > minimum_theoretical_abundance:
                n_missing += 1
        total += peak.intensity
        cleaned_eid.append(peak)
        i += 1

    tid = []
    for peak in base_theoretical:
        peak = peak.clone()
        peak.intensity *= total
        tid.append(peak)
    return cleaned_eid, tid, n_missing


class LCMSFeatureProcessorBase(object):
    def create_theoretical_distribution(self, mz, charge, charge_carrier=PROTON, truncate_after=0.8,
                                        ignore_below=0.05):
        base_tid = self.averagine.isotopic_cluster(
            mz, charge, truncate_after=truncate_after, charge_carrier=charge_carrier)
        total = 0
        kept_tid = []
        for i, p in enumerate(base_tid):
            if p.intensity < ignore_below and i > 1:
                continue
            else:
                total += p.intensity
                kept_tid.append(p)
        for p in kept_tid:
            p.intensity /= total
        return kept_tid

    def find_all_features(self, mz, error_tolerance=2e-5):
        return self.feature_map.find_all(mz, error_tolerance)

    def find_features(self, mz, error_tolerance=2e-5, interval=None):
        f = self.find_all_features(mz, error_tolerance)
        if not f:
            f = [None]
        elif interval is not None:
            f = [el for el in f if el.overlaps_in_time(interval)]
        return f

    def conform_envelopes(self, experimental, base_theoretical):
        return conform_envelopes(experimental, base_theoretical)

    def _find_thresholded_score(self, scores, percentage):
        scores = np.array(scores)
        maximum = scores.max()
        threshold = maximum * percentage
        return scores[scores > threshold].mean()

    def _fit_feature_set(self, mz, error_tolerance, charge,
                         charge_carrier=PROTON, truncate_after=0.8, max_missed_peaks=1,
                         threshold_scale=0.3, feature=None):
        base_tid = self.create_theoretical_distribution(mz, charge, charge_carrier, truncate_after)
        feature_groups = self.match_theoretical_isotopic_distribution(base_tid, error_tolerance, interval=feature)
        feature_fits = []
        for features in product(*feature_groups):
            if all(f is None for f in features):
                continue
            # If the monoisotopic feature wasn't actually observed, create a dummy feature
            # since the monoisotopic feature cannot be None
            if features[0] is None:
                features = list(features)
                features[0] = EmptyFeature(mz)
            feat_iter = FeatureSetIterator(features)
            scores = []
            times = []
            counter = 0
            for eid in feat_iter:
                cleaned_eid, tid, n_missing = self.conform_envelopes(eid, base_tid)
                if n_missing > max_missed_peaks:
                    continue
                score = self.scorer.evaluate(None, cleaned_eid, tid)
                if np.isnan(score):
                    continue
                scores.append(score)
                times.append(feat_iter.current_time)
                counter += 1
            final_score = self._find_thresholded_score(scores, threshold_scale)
            missing_features = 0
            for f in features:
                if f is None:
                    missing_features += 1
            fit = LCMSFeatureSetFit(
                features, base_tid, final_score, charge, missing_features,
                neutral_mass=neutral_mass(features[0].mz, charge, charge_carrier),
                missing_features=missing_features,
                scores=scores, times=times)
            if self.scorer.reject_score(fit.score):
                continue
            feature_fits.append(fit)
        return feature_fits

    def match_theoretical_isotopic_distribution(self, theoretical_distribution, error_tolerance=2e-5, interval=None):
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
        list of list of LCMSFeature
            The list of matched features for each theoretical m/z
        """
        experimental_distribution = [self.find_features(
            p.mz, error_tolerance, interval=interval) for p in theoretical_distribution]
        return experimental_distribution


try:
    has_c = True
    from ms_deisotope._c.feature_map.feature_processor import LCMSFeatureProcessorBase
except ImportError:
    has_c = False


def test_fit_features(self, features, charge, charge_carrier=PROTON, truncate_after=0.8, max_missed_peaks=1,
                      threshold_scale=0.3):
    base_tid = self.create_theoretical_distribution(features[0].mz, charge, charge_carrier, truncate_after)
    feat_iter = FeatureSetIterator(features)
    scores = []
    times = []
    counter = 0
    for eid in feat_iter:
        cleaned_eid, tid, n_missing = self.conform_envelopes(eid, base_tid)
        if n_missing > max_missed_peaks:
            continue
        score = self.scorer.evaluate(None, cleaned_eid, tid)
        if np.isnan(score):
            continue
        scores.append(score)
        times.append(feat_iter.current_time)
        counter += 1

    return np.array(scores), np.array(times)


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


def test_finalize(self, feature_fit, charge_carrier=PROTON, detection_threshold=0.1,
                  max_missed_peaks=1):
    start_time, end_time = find_bounds(feature_fit, detection_threshold)
    feat_iter = FeatureSetIterator(
        feature_fit.features, start_time, end_time)
    base_tid = feature_fit.theoretical
    charge = feature_fit.charge
    abs_charge = abs(charge)
    for eid in feat_iter:
        cleaned_eid, tid, n_missing = conform_envelopes(eid, base_tid)
        rep_eid = drop_placeholders(cleaned_eid)
        n_real_peaks = len(rep_eid)
        invalid = n_real_peaks == 0 or (n_real_peaks == 1 and abs_charge > 1) or n_missing > max_missed_peaks
        score = self.scorer.evaluate(None, cleaned_eid, tid)
        yield feat_iter.current_time, score, n_missing, invalid


class PrecursorMap(object):

    @classmethod
    def extract_ms1_scans_with_msn(cls, peak_loader):
        mapping = defaultdict(list)
        for msn_id, msn_info in peak_loader.extended_index.msn_ids.items():
            prec_id = msn_info['precursor_scan_id']
            mapping[prec_id].append(msn_info)
        return mapping

    @classmethod
    def map_precursor_information_to_peaks(cls, peak_loader, error_tolerance=2e-5):
        index = cls.extract_ms1_scans_with_msn(peak_loader)
        mapping = dict()
        missed = []
        for precursor_id, msn_info_set in index.items():
            precursor = peak_loader.get_scan_by_id(precursor_id)
            for msn_info in msn_info_set:
                peak = precursor.peak_set.has_peak(msn_info['mz'], error_tolerance)
                if peak is not None:
                    msn_info['intensity'] = peak.intensity
                    mapping[precursor.scan_time, peak.index] = msn_info
                else:
                    missed.append(msn_info)
        return mapping, missed

    @classmethod
    def build(cls, peak_loader, error_tolerance=2e-5):
        mapping, missed = cls.map_precursor_information_to_peaks(
            peak_loader, error_tolerance)
        return cls(mapping, missed)

    def __init__(self, mapping, missed=None):
        if missed is None:
            missed = []
        self.mapping = dict(mapping)
        self.missed = missed

    def assign_precursors(self, feature):
        hits = []
        for node in feature:
            for peak in node.members:
                key = (node.time, peak.index)
                if key in self.mapping:
                    pinfo = self.mapping.pop(key)
                    hits.append((node, pinfo))
        return hits

    def precursor_for_peak(self, time_peak_pair):
        if time_peak_pair in self.mapping:
            pinfo = self.mapping.pop(time_peak_pair)
        else:
            pinfo = None
        return pinfo


class LCMSFeatureProcessor(LCMSFeatureProcessorBase, LogUtilsMixin):
    def __init__(self, feature_map, averagine, scorer, precursor_map=None, minimum_size=3,
                 maximum_time_gap=0.25, prefer_multiply_charged=True):
        if precursor_map is None:
            precursor_map = PrecursorMap({})
        if isinstance(feature_map, LCMSFeatureMap):
            feature_map = feature_map.clone(deep=True)
        else:
            feature_map = LCMSFeatureMap(
                [f.clone(deep=True) for f in feature_map])
        self.feature_map = feature_map
        self.averagine = AveragineCache(averagine)
        self.prefer_multiply_charged = prefer_multiply_charged
        self.scorer = scorer
        self.precursor_map = precursor_map
        self.minimum_size = minimum_size
        self.maximum_time_gap = maximum_time_gap
        self.dependence_network = FeatureDependenceGraph(self.feature_map)
        self.orphaned_nodes = []

    def select_best_disjoint_subgraphs(self, disjoint_envelopes):
        solutions = []
        for cluster in disjoint_envelopes:
            if len(cluster) > 500:
                self.log("... Sub-graph with %d members: %r" % (len(cluster), cluster))
            disjoint_best_fits = cluster.disjoint_best_fits(max_size=5000)
            for fit in disjoint_best_fits:
                solutions.append(fit)
        return solutions

    def charge_state_determination(self, feature, error_tolerance, charge_range=(1, 8),
                                   left_search=1, right_search=1, charge_carrier=PROTON,
                                   truncate_after=0.8, max_missed_peaks=1, threshold_scale=0.3):
        fits = self.collect_all_fits(
            feature,
            error_tolerance=error_tolerance,
            charge_range=charge_range,
            left_search=left_search,
            right_search=right_search,
            charge_carrier=charge_carrier,
            truncate_after=truncate_after,
            max_missed_peaks=max_missed_peaks,
            threshold_scale=threshold_scale)
        return sorted(fits, key=lambda x: x.score, reverse=True)

    def collect_all_fits(self, feature, error_tolerance, charge_range=(1, 8),
                         left_search=1, right_search=1, charge_carrier=PROTON,
                         truncate_after=0.8, max_missed_peaks=1, threshold_scale=0.3):
        fits = []
        holdout = None
        if self.scorer.is_maximizing():
            best_fit_score = 0
        else:
            best_fit_score = float('inf')
        best_fit_charge = 0
        for charge in charge_range_(*charge_range):
            current_fits = self.fit_theoretical_distribution(
                feature,
                error_tolerance=error_tolerance,
                charge=charge,
                left_search=left_search,
                right_search=right_search,
                charge_carrier=charge_carrier,
                truncate_after=truncate_after,
                max_missed_peaks=max_missed_peaks,
                threshold_scale=threshold_scale)
            is_multiply_charged = abs(charge) > 1
            if self.scorer.is_maximizing():
                for fit in current_fits:
                    if fit.score > best_fit_score:
                        if is_multiply_charged and not fit.has_multiple_real_features():
                            continue
                        best_fit_score = fit.score
                        best_fit_charge = charge
            else:
                for fit in current_fits:
                    if fit.score < best_fit_score:
                        if is_multiply_charged and not fit.has_multiple_real_features():
                            continue
                        best_fit_score = fit.score
                        best_fit_charge = charge
            if abs(charge) == 1 and self.prefer_multiply_charged:
                holdout = current_fits
            else:
                for fit in current_fits:
                    if is_multiply_charged and not fit.has_multiple_real_features():
                        continue
                    self.dependence_network.add_fit_dependence(fit)
                    fits.append(fit)
        if holdout is not None and best_fit_charge == 1:
            for fit in holdout:
                self.dependence_network.add_fit_dependence(fit)
                fits.append(fit)
        return fits

    def fit_theoretical_distribution(self, feature, error_tolerance, charge, left_search=1, right_search=1,
                                     charge_carrier=PROTON, truncate_after=0.8, max_missed_peaks=1,
                                     threshold_scale=0.3):
        """Fits all theoretical distributions for `feature` with charge state `charge`.

        Parameters
        ----------
        feature : LCMSFeature
            The feature to use as the seed fit
        error_tolerance : float
            Permitted parts-per-million mass accuracy
        charge : int
            The charge state to fit at
        left_search : int, optional
            Description
        right_search : int, optional
            Description
        charge_carrier : float, optional
            Mass of the charge carrier
        truncate_after : float, optional
            Fraction of the total isotopic pattern to include
            in the set of theoretical peaks to search for.

        Returns
        -------
        list of LCMSFeatureSetFit
        """
        base_mz = feature.mz
        all_fits = []
        for offset in range(-left_search, right_search + 1):
            shift = isotopic_shift(charge) * offset
            mz = base_mz + shift
            feature_fits = self._fit_feature_set(
                mz,
                error_tolerance=error_tolerance,
                charge=charge,
                charge_carrier=charge_carrier,
                truncate_after=truncate_after,
                max_missed_peaks=max_missed_peaks,
                threshold_scale=threshold_scale,
                feature=feature)
            all_fits.extend(feature_fits)
        return all_fits

    def finalize_fit(self, feature_fit, charge_carrier=PROTON, subtract=True,
                     detection_threshold=0.1, max_missed_peaks=1):
        nodes = []
        start_time, end_time = find_bounds(feature_fit, detection_threshold)
        feat_iter = FeatureSetIterator(
            feature_fit.features, start_time, end_time)
        base_tid = feature_fit.theoretical
        charge = feature_fit.charge
        abs_charge = abs(charge)
        for eid in feat_iter:
            cleaned_eid, tid, n_missing = conform_envelopes(
                eid, base_tid.peaklist)
            rep_eid = drop_placeholders(cleaned_eid)
            n_real_peaks = len(rep_eid)
            if n_real_peaks == 0 or (n_real_peaks == 1 and abs_charge > 1) or \
               n_missing > max_missed_peaks:
                continue
            score = self.scorer.evaluate(None, cleaned_eid, tid)
            is_valid = True
            if np.isnan(score) or score < 0:
                is_valid = False
            envelope = [(e.mz, min(e.intensity, t.intensity)) for e, t in zip(cleaned_eid, tid)]
            if is_valid:
                total_abundance = sum(p[1] for p in envelope)
                monoisotopic_mass = neutral_mass(
                    base_tid.monoisotopic_mz, charge, charge_carrier=charge_carrier)
                reference_peak = first_peak(cleaned_eid)

                dpeak = DeconvolutedPeak(
                    neutral_mass=monoisotopic_mass, intensity=total_abundance,
                    charge=charge,
                    signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
                    index=reference_peak.index,
                    full_width_at_half_max=mean(p.full_width_at_half_max for p in rep_eid),
                    a_to_a2_ratio=a_to_a2_ratio(tid),
                    most_abundant_mass=neutral_mass(
                        most_abundant_mz(cleaned_eid), charge, charge_carrier=charge_carrier),
                    average_mass=neutral_mass(
                        average_mz(cleaned_eid), charge, charge_carrier=charge_carrier),
                    score=score,
                    envelope=envelope,
                    mz=base_tid.monoisotopic_mz,
                    area=sum(e.area for e in cleaned_eid))

                time = feat_iter.current_time
                precursor_info_set = []
                for peak in rep_eid:
                    pinfo = self.precursor_map.precursor_for_peak((time, peak.index))
                    if pinfo is not None:
                        precursor_info_set.append(pinfo)

                node = DeconvolutedLCMSFeatureTreeNode(time, [dpeak], precursor_info_set)
                nodes.append(node)
            if subtract:
                for fpeak, tpeak in zip(cleaned_eid, envelope):
                    # If a theoretical peak uses up more than 70%
                    # of the abundance of a single peak, this peak
                    # should not contribute meaninfully to any other
                    # fits from now on. Set it's abundance to 1.0 as
                    # if it were fully used up.
                    ruin = (fpeak.intensity * 0.7) < tpeak[1]
                    if ruin:
                        fpeak.intensity = 1.0
                    else:
                        fpeak.intensity -= tpeak[1]
                    if fpeak.intensity < 0:
                        fpeak.intensity = 1.0
        for feature in feature_fit.features:
            if feature is None or isinstance(feature, EmptyFeature):
                continue
            feature.invalidate()
        if len(nodes) < self.minimum_size:
            return None

        result_feature = DeconvolutedLCMSFeature(
            nodes, feature_fit.charge,
            score=feature_fit.score, n_features=len(feature_fit),
            supporters=feature_fit.supporters)

        return result_feature

    def _clean_solutions(self, solutions):
        out = []
        for solution in solutions:
            parts = solution.split_sparse(self.maximum_time_gap)
            for part in parts:
                if len(part) < self.minimum_size:
                    continue
                out.append(part)
        return out

    def _merge_isobaric(self, solutions):
        if len(solutions) == 0:
            return []
        solutions.sort(key=lambda x: (x.neutral_mass, x.charge))
        last = solutions[0]
        out = []
        for solution in solutions[1:]:
            if last.neutral_mass == solution.neutral_mass and last.charge == solution.charge:
                last.sum(solution)
            else:
                out.append(last)
                last = solution
        out.append(last)
        return out

    def remove_peaks_below_threshold(self, minimum_intensity):
        out = []
        for feature in self.feature_map:
            keep = []
            for node in feature:
                if node.total_intensity() > minimum_intensity:
                    keep.append(node)
                # else:
                #     self.orphaned_nodes.append(node)

            if keep:
                filtered_feature = LCMSFeature(keep, feature_id=feature.feature_id)
                out.extend(filter(
                    lambda f: len(f) >= self.minimum_size,
                    filtered_feature.split_sparse(
                        delta_rt=self.maximum_time_gap)))
        self.feature_map = LCMSFeatureMap(out)

    def store_solutions(self, fits, charge_carrier=PROTON, subtract=True, detection_threshold=0.1,
                        max_missed_peaks=1):
        solutions = []
        for fit in fits:
            extracted = extract_fitted_region(
                fit, detection_threshold=detection_threshold)
            if extracted is None:
                continue
            solution = self.finalize_fit(
                extracted, charge_carrier=charge_carrier, subtract=subtract,
                detection_threshold=detection_threshold, max_missed_peaks=max_missed_peaks)
            if solution is None:
                continue
            solutions.append(solution)
        return solutions

    def build_dependence_network(self):
        self.dependence_network = FeatureDependenceGraph(self.feature_map)

    def _map_precursors(self, error_tolerance):
        self.log("... Constructing Precursor Seeds")
        rt_map = RTMap(self.feature_map)
        cache = dict()
        seeds = set()

        for key, pinfo in (self.precursor_map.mapping.items()):
            time, ix = key
            try:
                candidates = cache[time]
            except KeyError:
                candidates = rt_map.rt_tree.contains_point(time)
                cache[time] = candidates
            mz = pinfo["mz"]
            hits = [c.members[0] for c in candidates if c.contains_mz(mz, error_tolerance)]
            seeds.update(hits)
        return seeds

    def _make_iterator_state(self, error_tolerance=2e-5, charge_range=(1, 8), left_search=1, right_search=0,
                             charge_carrier=PROTON, truncate_after=0.95, maxiter=10, minimum_intensity=100,
                             convergence=0.01, max_missed_peaks=1, threshold_scale=0.3, relfitter=None,
                             detection_threshold=0.1):
        state = FeatureDeconvolutionIterationState(
            self, error_tolerance, charge_range, left_search, right_search, charge_carrier,
            truncate_after, maxiter, minimum_intensity, convergence, max_missed_peaks, threshold_scale,
            relfitter, detection_threshold)
        return state

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), left_search=1, right_search=0,
                    charge_carrier=PROTON, truncate_after=0.95, maxiter=10, minimum_intensity=100,
                    convergence=0.01, max_missed_peaks=1, threshold_scale=0.3, relfitter=None,
                    detection_threshold=0.1):
        state = FeatureDeconvolutionIterationState(
            self, error_tolerance, charge_range, left_search, right_search, charge_carrier,
            truncate_after, maxiter, minimum_intensity, convergence, max_missed_peaks, threshold_scale,
            relfitter, detection_threshold)
        return state.run()


class FeatureDeconvolutionIterationState(LogUtilsMixin):
    def __init__(self, processor, error_tolerance=2e-5, charge_range=(1, 8), left_search=1, right_search=0,
                 charge_carrier=PROTON, truncate_after=0.95, maxiter=10, minimum_intensity=100,
                 convergence=0.01, max_missed_peaks=1, threshold_scale=0.3, relfitter=None,
                 detection_threshold=0.1):
        self.processor = processor
        self.last_signal_magnitude = sum(f.total_signal for f in self.processor.feature_map)
        self.next_signal_magnitude = 1.0
        self.total_signal_ratio = 1.0
        self.generation = None
        self.relations = None
        self.fits = None
        self.all_fits = None
        self.solutions = []
        self.disjoint_feature_clusters = None
        self.maxiter = maxiter
        self.iteration_count = 0

        self.error_tolerance = error_tolerance
        self.charge_range = charge_range
        self.left_search = left_search
        self.right_search = right_search
        self.charge_carrier = charge_carrier
        self.truncate_after = truncate_after
        self.minimum_intensity = minimum_intensity
        self.convergence = convergence
        self.max_missed_peaks = max_missed_peaks
        self.threshold_scale = threshold_scale
        self.relfitter = relfitter
        self.detection_threshold = detection_threshold
        self.debug = False

    def update_signal_ratio(self):
        self.next_signal_magnitude = sum(f.total_signal for f in self.processor.feature_map)
        self.total_signal_ratio = (self.last_signal_magnitude - self.next_signal_magnitude) / self.next_signal_magnitude
        self.log("Signal Ratio: %0.3e (%0.3e, %0.3e)" % (
            self.total_signal_ratio, self.last_signal_magnitude, self.next_signal_magnitude))
        self.last_signal_magnitude = self.next_signal_magnitude
        return self.total_signal_ratio

    def setup(self):
        self.log("Begin Iteration %d" % (self.iteration_count,))
        self.log("Total Signal: %0.3e" % (sum(f.total_signal for f in self.processor.feature_map),))
        self.processor.remove_peaks_below_threshold(self.minimum_intensity)
        self.processor.build_dependence_network()

    def map_fits(self, report_interval=5):
        i = 0
        n = len(self.processor.feature_map)
        interval = int(max(n // report_interval, 1000))
        for feature in sorted(self.processor.feature_map, key=lambda x: x.mz, reverse=True):
            fits = self.processor.charge_state_determination(
                feature,
                error_tolerance=self.error_tolerance,
                charge_range=self.charge_range,
                left_search=self.left_search,
                right_search=self.right_search,
                charge_carrier=self.charge_carrier,
                truncate_after=self.truncate_after,
                max_missed_peaks=self.max_missed_peaks,
                threshold_scale=self.threshold_scale)
            i += 1
            if i % interval == 0:
                self.log("... %0.1f%%" % ((100. * i) / n,))
        self.all_fits = list(self.processor.dependence_network.dependencies)
        self.log("... Fit %d Theoretical Patterns" % len(self.all_fits))
        self.disjoint_feature_clusters = self.processor.dependence_network.find_non_overlapping_intervals()
        self.log("... Found %d Disjoint Feature Subgraphs" % len(self.disjoint_feature_clusters))

        if self.relfitter is not None:
            self.relations = self.relfitter.fit(
                (d for cluster in self.disjoint_feature_clusters
                 for d in cluster), self.solutions)
            self.relfitter.predict((d for cluster in self.disjoint_feature_clusters
                                    for d in cluster))
        self.log("... Extracting Best Fits")
        self.fits = self.processor.select_best_disjoint_subgraphs(self.disjoint_feature_clusters)

    def postprocess(self, subtract=True, detection_threshold=0.1):
        self.generation = self.processor.store_solutions(
            self.fits, charge_carrier=self.charge_carrier,
            subtract=subtract, detection_threshold=detection_threshold,
            max_missed_peaks=self.max_missed_peaks)

        if self.relfitter is not None:
            self.relfitter.add_history((self.relations, self.generation))
        self.solutions.extend(self.generation)

    def step(self, report_interval=5, subtract=True):
        self.setup()

        self.map_fits(report_interval=report_interval)
        self.postprocess(
            subtract=subtract, detection_threshold=self.detection_threshold)

        self.iteration_count += 1
        self.update_signal_ratio()

        if self.total_signal_ratio < self.convergence:
            return True
        else:
            return False

    def run(self, callback=None):
        keep_going = True
        is_callable_callback = callable(callback)
        while keep_going:
            converged = self.step()
            if is_callable_callback:
                callback(self)
            if converged or self.iteration_count >= self.maxiter:
                keep_going = False
        self.solutions = self.processor._clean_solutions(self.solutions)
        return DeconvolutedLCMSFeatureMap(
            smooth_overlaps_neutral(
                self.solutions,
                self.error_tolerance))


def find_bounds(fit, detection_threshold=0.1, find_separation=True):
    start_time = 0
    end_time = float('inf')

    for f, p in zip(fit, fit.theoretical):
        if f is None:
            continue
        passed_threshold = p.intensity >= detection_threshold
        if f.start_time > start_time and passed_threshold:
            start_time = f.start_time
        if f.end_time < end_time and passed_threshold:
            end_time = f.end_time

    if fit.n_points > 0 and find_separation:
        last_score = float('inf')
        begin_i = 0
        end_i = len(fit.scores) - 1
        smoothed_scores = smooth_leveled(fit.times, fit.scores, 3)
        for i, score in enumerate(smoothed_scores):
            if score > 0 and last_score < 0:
                begin_i = i
            elif score < 0 and last_score > 0:
                end_i = i
            last_score = score
            pass
        if end_i < begin_i:
            end_i = len(fit.scores) - 1

        if start_time < fit.times[begin_i] and begin_i != 0:
            start_time = fit.times[begin_i]

        # if end_time > fit.times[end_i]:
        #     end_time = fit.times[end_i]

    return start_time, end_time


def extract_fitted_region(feature_fit, detection_threshold=0.1):
    fitted_features = []
    if feature_fit.n_points == 0:
        return None
    start_time, end_time = find_bounds(feature_fit, detection_threshold)
    for feature in feature_fit:
        if feature is None or isinstance(feature, EmptyFeature):
            fitted_features.append(feature)
            continue
        _, start_ix = feature.nodes.find_time(start_time)
        _, end_ix = feature.nodes.find_time(end_time)
        if feature.nodes[end_ix].time < end_time:
            end_ix += 1
        fitted = LCMSFeature(feature.nodes[start_ix:(end_ix + 1)], feature_id=feature.feature_id)
        if len(fitted) == 0:
            fitted_features.append(EmptyFeature(fitted.mz))
        else:
            fitted_features.append(fitted)
    start_ix = binsearch(feature_fit.times, start_time)
    end_ix = binsearch(feature_fit.times, end_time)
    extract = LCMSFeatureSetFit(
        fitted_features,
        feature_fit.theoretical,
        feature_fit.score,
        feature_fit.charge,
        feature_fit.missing_features,
        supporters=feature_fit.supporters,
        data=feature_fit.data,
        neutral_mass=feature_fit.neutral_mass,
        scores=feature_fit.scores[start_ix:end_ix + 1],
        times=feature_fit.times[start_ix:end_ix + 1])
    return extract


class RTFeatureNode(Interval):
    def __init__(self, feature):
        super(RTFeatureNode, self).__init__(
            feature.start_time, feature.end_time, [feature], mz=feature.mz)
        self.mz = feature.mz

    def contains_mz(self, mz, time, error_tolerance=2e-5):
        return abs((self.mz - mz) / mz) < error_tolerance


class RTMap(object):
    def __init__(self, features):
        self.rt_tree = IntervalTreeNode.build(map(RTFeatureNode, features))

    def locate_precursor_candidates(self, precursor_info):
        time = precursor_info.precursor.scan_time
        candidates = self.rt_tree.contains_point(time)
        return candidates

    def find_precursor(self, precursor_info, error_tolerance=2e-5):
        time = precursor_info.precursor.scan_time
        candidates = self.rt_tree.contains_point(time)
        matches = []
        for candidate in candidates:
            if candidate.contains_mz(precursor_info.mz, time, error_tolerance):
                matches.append(candidate)
        return matches

    def find_features_at_time(self, time):
        hits = self.rt_tree.contains_point(time)
        return [
            node.members[0] for node in hits
        ]


try:
    has_c = True
    _conform_envelopes = conform_envelopes
except ImportError:
    has_c = False
