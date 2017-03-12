from collections import defaultdict
from itertools import product

import numpy as np

from ms_peak_picker import FittedPeak

from .feature_map import (
    LCMSFeatureMap,
    DeconvolutedLCMSFeatureMap)
from .lcms_feature import (
    LCMSFeature,
    EmptyFeature,
    FeatureSetIterator,
    NodeFeatureSetIterator)
from .feature_fit import (
    LCMSFeatureSetFit,
    DeconvolutedLCMSFeatureTreeNode,
    DeconvolutedLCMSFeature)
from .dependence_network import FeatureDependenceGraph
from ms_deisotope.peak_dependency_network.intervals import Interval, IntervalTreeNode
from ms_deisotope.averagine import AveragineCache, PROTON, isotopic_shift, mass_charge_ratio
from ms_deisotope.deconvolution import (
    charge_range_, DeconvolutedPeak, drop_placeholders,
    neutral_mass, average_mz, a_to_a2_ratio, first_peak,
    most_abundant_mz, mean)
from ms_deisotope.utils import printer


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


def sum_deconvoluted_features(last, solution):
    missing = []
    feat_iter = NodeFeatureSetIterator([last, solution])
    for nodes in feat_iter:
        base = nodes[0]
        new = nodes[1]
        if base is None:
            missing.append(new)
        elif new is not None:
            base.members[0].intensity += new.members[0].intensity
            base.precursor_information.extend(new.precursor_information)
    if missing:
        for node in missing:
            last.insert_node(DeconvolutedLCMSFeatureTreeNode(
                node.retention_time, list(node.members), list(node.precursor_information)))
    # last.supporters.extend(solution.supporters)
    return last


class LCMSFeatureProcessorBase(object):
    def create_theoretical_distribution(self, mz, charge, charge_carrier=PROTON, truncate_after=0.8):
        base_tid = self.averagine.isotopic_cluster(
            mz, charge, truncate_after=truncate_after, charge_carrier=charge_carrier)
        return base_tid

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
            final_score = sum(scores) / len(scores)
            missing_features = 0
            for f in features:
                if f is None:
                    missing_features += 1
            fit = LCMSFeatureSetFit(features, base_tid, final_score, charge, missing_features, data=scores)
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
    from ms_deisotope._c.feature_map.feature_processor import LCMSFeatureProcessorBase
except ImportError:
    pass


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
                key = (node.retention_time, peak.index)
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


class LCMSFeatureProcessor(LCMSFeatureProcessorBase):
    def __init__(self, feature_map, averagine, scorer, precursor_map=None):
        if precursor_map is None:
            precursor_map = PrecursorMap({})
        self.feature_map = LCMSFeatureMap([f.clone(deep=True) for f in feature_map])
        self.averagine = AveragineCache(averagine)
        self.scorer = scorer
        self.dependence_network = FeatureDependenceGraph(self.feature_map)
        self.precursor_map = precursor_map
        self.orphaned_nodes = []

    def select_best_disjoint_subgraphs(self, disjoint_envelopes):
        solutions = []
        for cluster in disjoint_envelopes:
            disjoint_best_fits = cluster.disjoint_best_fits()
            for fit in disjoint_best_fits:
                solutions.append(fit)
        return solutions

    def charge_state_determination(self, feature, error_tolerance, charge_range=(1, 8),
                                   left_search=1, right_search=1, charge_carrier=PROTON,
                                   truncate_after=0.8):
        fits = []
        for charge in charge_range_(*charge_range):
            fits.extend(self.fit_theoretical_distribution(
                feature, 2e-5, charge, left_search=left_search,
                right_search=right_search, truncate_after=truncate_after))
        return sorted(fits, key=lambda x: x.score, reverse=True)

    def fit_theoretical_distribution(self, feature, error_tolerance, charge, left_search=1, right_search=1,
                                     charge_carrier=PROTON, truncate_after=0.8):
        base_mz = feature.mz
        all_fits = []
        for offset in range(-left_search, right_search + 1):
            shift = isotopic_shift(charge) * offset
            mz = base_mz + shift
            feature_fits = self._fit_feature_set(
                mz, error_tolerance, charge, left_search,
                right_search, charge_carrier, truncate_after,
                feature)
            all_fits.extend(feature_fits)
            for fit in feature_fits:
                self.dependence_network.add_fit_dependence(fit)
        return all_fits

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
            score=feature_fit.score, n_features=len(feature_fit),
            supporters=feature_fit.supporters)

    def _merge_isobaric(self, solutions):
        solutions.sort(key=lambda x: (x.neutral_mass, x.charge))
        last = solutions[0]
        out = []
        for solution in solutions[1:]:
            if last.neutral_mass == solution.neutral_mass and last.charge == solution.charge:
                sum_deconvoluted_features(last, solution)
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
                else:
                    self.orphaned_nodes.append(node)

            if keep:
                filtered_feature = LCMSFeature(keep, feature_id=feature.feature_id)
                out.extend(filtered_feature.split_sparse())
        self.feature_map = LCMSFeatureMap(out)

    def store_solutions(self, fits, charge_carrier=PROTON):
        solutions = []
        for fit in fits:
            extracted = extract_fitted_region(fit)
            solution = self.finalize_fit(extracted, charge_carrier=charge_carrier)
            if solution is None:
                continue
            solutions.append(solution)
        return solutions

    def build_dependence_network(self):
        self.dependence_network = FeatureDependenceGraph(self.feature_map)

    def _map_precursors(self, error_tolerance):
        printer("\tConstructing Precursor Seeds")
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

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), left_search=1, right_search=1,
                    charge_carrier=PROTON, truncate_after=0.8, maxiter=10, minimum_intensity=500,
                    convergence=0.01, relfitter=None):
        solutions = []
        last_signal_magnitude = sum(f.total_signal for f in self.feature_map)
        next_signal_magnitude = 1.0
        total_signal_ratio = 1.0
        relations = None
        generation = None
        for j in range(maxiter):
            printer("Begin Iteration %d" % (j,))
            printer("Total Signal: %0.3e" % (sum(f.total_signal for f in self.feature_map),))
            self.remove_peaks_below_threshold(minimum_intensity)
            self.build_dependence_network()

            i = 0
            n = len(self.feature_map)
            interval = int(max(n // 5, 1000))
            for feature in sorted(self.feature_map, key=lambda x: x.mz, reverse=True):
                self.charge_state_determination(
                    feature, error_tolerance, charge_range, left_search, right_search,
                    charge_carrier, truncate_after)
                i += 1
                if i % interval == 0:
                    printer("\t%0.1f%%" % ((100. * i) / n,))
            disjoint_feature_clusters = self.dependence_network.find_non_overlapping_intervals()
            if relfitter is not None:
                relations = relfitter.fit((d for cluster in disjoint_feature_clusters for d in cluster), solutions)
            printer("\tExtracting Fits")
            fits = self.select_best_disjoint_subgraphs(disjoint_feature_clusters)
            generation = self.store_solutions(fits, charge_carrier=charge_carrier)
            if relfitter is not None:
                relfitter.add_history((relations, generation))
            solutions.extend(generation)

            next_signal_magnitude = sum(f.total_signal for f in self.feature_map)
            total_signal_ratio = (last_signal_magnitude - next_signal_magnitude) / next_signal_magnitude
            printer("Signal Ratio: %0.3e (%0.3e, %0.3e)" % (
                total_signal_ratio, last_signal_magnitude, next_signal_magnitude))
            if total_signal_ratio < convergence:
                break
            else:
                last_signal_magnitude = next_signal_magnitude
        else:
            printer("Failed to converge after %d iterations. (%0.4f)" % (j, total_signal_ratio))
        solutions = self._merge_isobaric(solutions)
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
        fitted = LCMSFeature(feature.nodes[start_ix:(end_ix + 1)], feature_id=feature.feature_id)
        if len(fitted) == 0:
            fitted_features.append(EmptyFeature(fitted.mz))
        else:
            fitted_features.append(fitted)
    return LCMSFeatureSetFit(
        fitted_features,
        feature_fit.theoretical,
        feature_fit.score,
        feature_fit.charge,
        feature_fit.missing_features,
        supporters=feature_fit.supporters,
        data=feature_fit.data,
        neutral_mass=feature_fit.neutral_mass)


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
    from ms_deisotope._c.feature_map.feature_processor import conform_envelopes
except ImportError:
    has_c = False
