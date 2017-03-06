from libc.stdlib cimport malloc, free

from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append
from ms_peak_picker._c.peak_set cimport FittedPeak
from brainpy._c.isotopic_distribution cimport TheoreticalPeak

from ms_deisotope._c.scoring cimport IsotopicFitterBase
from ms_deisotope._c.averagine cimport AveragineCache, PROTON

from ms_deisotope._c.feature_map.feature_map cimport LCMSFeatureMap
from ms_deisotope._c.feature_map.feature_fit cimport LCMSFeatureSetFit
from ms_deisotope._c.feature_map.lcms_feature cimport (
    LCMSFeature,
    FeatureSetIterator,
    EmptyFeature)

cimport numpy as cnp
import numpy as np


cdef class product_iterator(object):
    cdef:
        public list collections

    def __init__(self, collections):
        self.collections = list(map(list, collections))
        self.count = PyList_GET_SIZE(self.collections)

    cpdef list product(self):
        cdef:
            list result, next_result, new, layer, entry
            size_t layer_i, layer_n, entry_i, entry_n, value_i, value_n

        result = [[]]
        layer_n = self.count
        for layer_i in range(layer_n):
            layer = <list>PyList_GET_ITEM(self.collections, layer_i)
            next_result = []
            entry_n = PyList_GET_SIZE(result)
            for entry_i in range(entry_n):
                entry = <list>PyList_GET_ITEM(result, entry_i)
                value_n = PyList_GET_SIZE(layer)
                for value_i in range(value_n):
                    value = <object>PyList_GET_ITEM(layer, value_i)
                    new = list(entry)
                    new.append(value)
                    next_result.append(new)
            result = next_result
        return result

cdef list product(collections):
    cdef product_iterator piter
    piter = product_iterator(collections)
    return piter.product()


cdef tuple _conform_envelopes(list experimental, list base_theoretical):
    cdef:
        double total = 0
        size_t n_missing = 0
        size_t i = 0
        list cleaned_eid
        list tid
        FittedPeak fpeak
        TheoreticalPeak tpeak, peak

    cleaned_eid = []
    
    for i in range(PyList_GET_SIZE(experimental)):
        fpeak = <FittedPeak>PyList_GET_ITEM(experimental, i)
        if fpeak is None:
            peak = <TheoreticalPeak>PyList_GET_ITEM(base_theoretical, i)
            fpeak = FittedPeak(peak.mz, 1, 1, -1, -1, 0, 1, 0, 0)
            n_missing += 1
        total += fpeak.intensity
        cleaned_eid.append(fpeak)
        i += 1

    tid = []
    for i in range(PyList_GET_SIZE(base_theoretical)):
        tpeak = <TheoreticalPeak>PyList_GET_ITEM(base_theoretical, i)
        peak = tpeak.clone()
        peak.intensity *= total
        tid.append(peak)
    return cleaned_eid, tid, n_missing


cpdef tuple conform_envelopes(list experimental, list base_theoretical):
    return _conform_envelopes(experimental, base_theoretical)


cdef int validate_feature_set(list features):
    cdef:
        size_t i, n, cnt
        object feature

    n = PyList_GET_SIZE(features)
    cnt = 0
    for i in range(n):
        feature = <object>PyList_GET_ITEM(features, i)
        if feature is None:
            cnt += 1
    return cnt != n


cdef class FeatureProcessorBase(object):
    cdef:
        public LCMSFeatureMap feature_map
        public IsotopicFitterBase scorer
        public AveragineCache averagine

    cpdef list create_theoretical_distribution(self, double mz, int charge, double charge_carrier=PROTON, double truncate_after=0.8):
        base_tid = self.averagine.isotopic_cluster(
            mz, charge, truncate_after=truncate_after, charge_carrier=charge_carrier)
        return base_tid

    cpdef list find_all_features(self, double mz, double error_tolerance=2e-5):
        return self.feature_map._find_all(mz, error_tolerance)

    cpdef list find_features(self, double mz, double error_tolerance=2e-5, LCMSFeature interval=None):
        cdef:
            list f, overlapped
            LCMSFeature el
            size_t i, n
        f = self.find_all_features(mz, error_tolerance)
        if not f:
            f = [None]
        elif interval is not None:
            overlapped = []
            n = PyList_GET_SIZE(overlapped)
            for i in range(n):
                el = <LCMSFeature>PyList_GET_ITEM(f, i)
                if el.overlaps_in_time(interval):
                    overlapped.append(el)
            f = overlapped
        return f

    cpdef list match_theoretical_isotopic_distribution(self, list theoretical_distribution, double error_tolerance=2e-5, LCMSFeature interval=None):
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
        cdef:
            list experimental_distribution, found
            size_t i, n
            TheoreticalPeak p

        experimental_distribution = []
        n = PyList_GET_SIZE(theoretical_distribution)
        for i in range(n):
            p = <TheoreticalPeak>PyList_GET_ITEM(theoretical_distribution, i)
            found = self.find_features(p.mz, error_tolerance, interval=interval)
            experimental_distribution.append(found)
        return experimental_distribution

    cpdef list _fit_feature_set(self, mz, error_tolerance, charge, left_search=1, right_search=1,
                                charge_carrier=PROTON, truncate_after=0.8, feature=None):
        cdef:
            list base_tid, feature_groups, feature_fits, features, combn
            list scores, eid, cleaned_eid, tid, 
            size_t combn_size, combn_i, n_missing, missing_features
            double score, final_score
            FeatureSetIterator feat_iter


        base_tid = self.create_theoretical_distribution(mz, charge, charge_carrier, truncate_after)
        feature_groups = self.match_theoretical_isotopic_distribution(base_tid, error_tolerance, interval=feature)
        feature_fits = []

        combn = product(feature_groups)
        combn_size = PyList_GET_SIZE(combn)

        for combn_i in range(combn_size):
            features = <list>PyList_GET_ITEM(combn, combn_i)
            if not validate_feature_set(features):
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
