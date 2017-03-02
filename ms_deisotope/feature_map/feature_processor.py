from itertools import product

from ms_peak_picker import FittedPeak

from .lcms_feature import LCMSFeature, FeatureSetIterator
from .feature_fit import LCMSFeatureSetFit
from ms_deisotope.scoring import PenalizedMSDeconVFitter
from ms_deisotope.averagine import AveragineCache, PROTON


class LCMSFeatureProcessor(object):
    def __init__(self, feature_map, averagine, scorer):
        self.feature_map = feature_map
        self.averagine = AveragineCache(averagine)
        self.scorer = scorer

    def has_feature(self, mz, error_tolerance=2e-5):
        return self.feature_map.search(mz, error_tolerance)

    def find_all_features(self, mz, error_tolerance=2e-5):
        return self.feature_map.find_all(mz, error_tolerance)

    def conform_envelopes(self, experimental, base_theoretical):
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

    def create_theoretical_distribution(self, mz, charge, charge_carrier=PROTON, truncate_after=0.8):
        base_tid = self.averagine.isotopic_cluster(
            mz, charge, truncate_after=truncate_after, charge_carrier=charge_carrier)
        return base_tid

    def fit_theoretical_distribution(self, feature, error_tolerance, charge, charge_carrier=PROTON, truncate_after=0.8):
        mz = feature.mz
        base_tid = self.create_theoretical_distribution(mz, charge, charge_carrier, truncate_after)
        feature_groups = []
        for tpeak in base_tid:
            f = self.find_all_features(tpeak.mz, 2e-5)
            if not f:
                f = [None]
            feature_groups.append(f)
        feature_fits = []
        for features in product(*feature_groups):
            feat_iter = FeatureSetIterator(features)
            scores = []
            for eid in feat_iter:
                cleaned_eid, tid, n_missing = self.conform_envelopes(eid, base_tid)
                assert len(cleaned_eid) == len(tid)
                score = self.scorer.evaluate(None, cleaned_eid, tid)
                scores.append(score)
            final_score = sum(scores)
            # TODO - properly account for missing feature
            fit = LCMSFeatureSetFit(features, base_tid, final_score, charge, 0)
            feature_fits.append(fit)
        return feature_fits

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
