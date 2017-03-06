from itertools import product

import numpy as np

from ms_peak_picker import FittedPeak

from .feature_map import (
    LCMSFeatureMap,
    DeconvolutedLCMSFeatureMap)
from .lcms_feature import (
    LCMSFeature,
    EmptyFeature,
    FeatureSetIterator)
from .feature_fit import (
    LCMSFeatureSetFit,
    DeconvolutedLCMSFeatureTreeNode,
    DeconvolutedLCMSFeature)
from .dependence_network import FeatureDependenceGraph
from ms_deisotope.averagine import AveragineCache, PROTON, isotopic_shift, mass_charge_ratio
from ms_deisotope.deconvolution import (
    charge_range_, DeconvolutedPeak, drop_placeholders,
    neutral_mass, average_mz, a_to_a2_ratio, first_peak,
    most_abundant_mz, mean)


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


def conform_envelopes(experimental, base_theoretical):
    total = 0
    n_missing = 0
    i = 0
    cleaned_eid = []
    for peak in experimental:
        if peak is None:
            peak = FittedPeak(base_theoretical[i].mz, 1, 1, -1, -1, 0, 1, 0, 0)
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


def merge_deconvoluted_features(last, solution):
    missing = []
    feat_iter = FeatureSetIterator([last, solution])
    for peaks in feat_iter:
        base = peaks[0]
        new = peaks[1]
        if base is None:
            missing.append((new, feat_iter.current_time))
        if new is not None:
            base.intensity += new.intensity
    if missing:
        for peak, time in missing:
            last.insert_node(DeconvolutedLCMSFeatureTreeNode(time, [peak]))
    return last


class LCMSFeatureProcessor(object):
    def __init__(self, feature_map, averagine, scorer):
        self.feature_map = LCMSFeatureMap([f.clone() for f in feature_map])
        self.averagine = AveragineCache(averagine)
        self.scorer = scorer
        self.dependence_network = FeatureDependenceGraph(self.feature_map)

    def find_all_features(self, mz, error_tolerance=2e-5):
        return self.feature_map.find_all(mz, error_tolerance)

    def conform_envelopes(self, experimental, base_theoretical):
        return conform_envelopes(experimental, base_theoretical)

    def select_best_disjoint_subgraphs(self, error_tolerance=2e-5, charge_carrier=PROTON):
        disjoint_envelopes = self.dependence_network.find_non_overlapping_intervals()

        solutions = []
        for cluster in disjoint_envelopes:
            disjoint_best_fits = cluster.disjoint_best_fits()
            for fit in disjoint_best_fits:
                solutions.append(fit)
        return solutions

    def create_theoretical_distribution(self, mz, charge, charge_carrier=PROTON, truncate_after=0.8):
        base_tid = self.averagine.isotopic_cluster(
            mz, charge, truncate_after=truncate_after, charge_carrier=charge_carrier)
        return base_tid

    def find_features(self, mz, error_tolerance=2e-5, interval=None):
        f = self.find_all_features(mz, error_tolerance)
        if not f:
            f = [None]
        elif interval is not None:
            f = [el for el in f if el.overlaps_in_time(interval)]
        return f

    def _fit_feature_set(self, mz, error_tolerance, charge, left_search=1, right_search=1,
                         charge_carrier=PROTON, truncate_after=0.8, feature=None):
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
            for eid in feat_iter:
                cleaned_eid, tid, n_missing = self.conform_envelopes(eid, base_tid)
                score = self.scorer.evaluate(None, cleaned_eid, tid)
                if np.isnan(score):
                    continue
                scores.append(score)
            final_score = sum(scores)
            missing_features = 0
            for f in features:
                if f is None:
                    missing_features += 1
            fit = LCMSFeatureSetFit(features, base_tid, final_score, charge, missing_features, )
            if self.scorer.reject_score(fit.score):
                continue
            feature_fits.append(fit)
        return feature_fits

    def fit_theoretical_distribution(self, feature, error_tolerance, charge, left_search=1, right_search=1,
                                     charge_carrier=PROTON, truncate_after=0.8):
        base_mz = feature.mz
        for offset in range(-left_search, right_search + 1):
            mz = base_mz + isotopic_shift(charge) * offset
            feature_fits = self._fit_feature_set(
                mz, error_tolerance, charge, left_search, right_search,
                charge_carrier, truncate_after, feature)
            for fit in feature_fits:
                self.dependence_network.add_fit_dependence(fit)
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

    def charge_state_determination(self, feature, error_tolerance, charge_range=(1, 8),
                                   left_search=1, right_search=1, charge_carrier=PROTON,
                                   truncate_after=0.8):
        fits = []
        for charge in charge_range_(*charge_range):
            fits.extend(self.fit_theoretical_distribution(
                feature, 2e-5, charge, left_search=left_search,
                right_search=right_search))
        return sorted(fits, key=lambda x: x.score, reverse=True)

    def finalize_fit(self, feature_fit, charge_carrier=PROTON, subtract=True):
        nodes = []
        feat_iter = FeatureSetIterator(feature_fit.features)
        base_tid = feature_fit.theoretical
        charge = feature_fit.charge
        for eid in feat_iter:
            cleaned_eid, tid, n_missing = conform_envelopes(eid, base_tid)
            rep_eid = drop_placeholders(cleaned_eid)
            if len(rep_eid) == 0:
                continue
            score = self.scorer.evaluate(None, cleaned_eid, tid)
            envelope = [(e.mz, min(e.intensity, t.intensity)) for e, t in zip(cleaned_eid, tid)]
            total_abundance = sum(p[1] for p in envelope)
            monoisotopic_mass = neutral_mass(
                cleaned_eid[0].mz, charge, charge_carrier=charge_carrier)
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
                mz=cleaned_eid[0].mz,
                area=sum(e.area for e in cleaned_eid))

            node = DeconvolutedLCMSFeatureTreeNode(feat_iter.current_time, [dpeak])
            nodes.append(node)
            if subtract:
                for fpeak, tpeak in zip(cleaned_eid, envelope):
                    fpeak.intensity -= tpeak[1]
                    if fpeak.intensity < 0:
                        fpeak.intensity = 1.0
        for feature in feature_fit.features:
            if feature is None or isinstance(feature, EmptyFeature):
                continue
            feature.invalidate()
        if len(nodes) == 0:
            return None
        return DeconvolutedLCMSFeature(
            nodes, feature_fit.charge,
            score=feature_fit.score, n_features=len(feature_fit))

    def _merge_isobaric(self, solutions):
        solutions.sort(key=lambda x: x.neutral_mass)
        last = solutions[0]
        out = []
        for solution in solutions[1:]:
            if last.neutral_mass == solution.neutral_mass and last.charge == solution.charge:
                merge_deconvoluted_features(last, solution)
            else:
                out.append(last)
                last = solution
        out.append(last)
        return out

    def remove_peaks_below_threshold(self, minimum_intensity):
        out = []
        for feature in self.feature_map:
            nodes = [node for node in feature if node.total_intensity() > minimum_intensity]
            if nodes:
                filtered_feature = LCMSFeature(nodes)
                out.append(filtered_feature)
        self.feature_map = LCMSFeatureMap(out)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), left_search=1, right_search=1,
                    charge_carrier=PROTON, truncate_after=0.8, iterations=1, minimum_intensity=500):
        solutions = []
        for j in range(iterations):
            print("Begin Iteration %d" % (j,))
            self.remove_peaks_below_threshold(minimum_intensity)
            self.dependence_network = FeatureDependenceGraph(self.feature_map)
            i = 0
            n = len(self.feature_map)
            for feature in sorted(self.feature_map, key=lambda x: x.mz, reverse=True):
                self.charge_state_determination(
                    feature, error_tolerance, charge_range, left_search, right_search,
                    charge_carrier, truncate_after)
                i += 1
                if i % 1000 == 0:
                    print("\t%0.2f%%" % ((100. * i) / n,))
            print("\tExtracting Fits")
            fits = self.select_best_disjoint_subgraphs(error_tolerance, charge_carrier)
            for fit in fits:
                extracted = extract_fitted_region(fit)
                solution = self.finalize_fit(extracted, charge_carrier=charge_carrier)
                if solution is None:
                    continue
                solutions.append(solution)
            print("\tRemoving Signal Below Threshold")
        # solutions = self._merge_isobaric(solutions)
        return DeconvolutedLCMSFeatureMap(solutions)


def extract_fitted_region(feature_fit):
    fitted_features = []
    start_time = max([f.start_time for f in feature_fit if f is not None])
    end_time = min([f.end_time for f in feature_fit if f is not None])
    for feature in feature_fit:
        if feature is None or isinstance(feature, EmptyFeature):
            fitted_features.append(feature)
            continue
        _, start_ix = feature.nodes.find_time(start_time)
        _, end_ix = feature.nodes.find_time(end_time)
        if feature.nodes[end_ix].retention_time < end_time:
            end_ix += 1
        fitted = LCMSFeature(feature.nodes[start_ix:(end_ix + 1)])
        if len(fitted) == 0:
            fitted_features.append(None)
        else:
            fitted_features.append(fitted)
    return LCMSFeatureSetFit(
        fitted_features, feature_fit.theoretical, feature_fit.score,
        feature_fit.charge, feature_fit.missing_features,
        feature_fit.data)


try:
    has_c = True
    _conform_envelopes = conform_envelopes
    from ms_deisotope._c.feature_map.feature_processor import conform_envelopes
except ImportError:
    has_c = False
