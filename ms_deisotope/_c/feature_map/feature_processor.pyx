# cython: embedsignature=True
# 

from libc.stdlib cimport malloc, free

cdef extern from "numpy/npy_math.h":
    bint npy_isnan(double x)

cimport cython

from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append, PyList_GetItem, PyList_SetItem, PyList_New
from ms_peak_picker._c.peak_set cimport FittedPeak
from brainpy._c.isotopic_distribution cimport TheoreticalPeak
from brainpy._c.double_vector cimport (
    DoubleVector as dvec,
    make_double_vector,
    free_double_vector,
    double_vector_append,
    double_vector_to_list)

from ms_deisotope._c.scoring cimport IsotopicFitterBase
from ms_deisotope._c.averagine cimport AveragineCache, PROTON, neutral_mass as calc_neutral_mass

from ms_deisotope._c.feature_map.feature_map cimport LCMSFeatureMap
from ms_deisotope._c.feature_map.feature_fit cimport LCMSFeatureSetFit
from ms_deisotope._c.feature_map.lcms_feature cimport (
    LCMSFeature,
    FeatureSetIterator,
    EmptyFeature)

cimport numpy as np
import numpy as np


cdef object zeros = np.zeros


cdef double INF = float('inf')


cdef object double_vector_to_ndarray(dvec* vec):
    cdef:
        np.ndarray[double, ndim=1, mode='c'] out
        np.npy_intp shape[1]
        size_t i, n

    n = shape[0] = vec.used
    # out = np.PyArray_SimpleNew(1, shape, np.NPY_DOUBLE)
    out = zeros(shape[0])
    for i in range(n):
        out[i] = vec.v[i]
    return out


cdef class CartesianProductIterator(object):
    cdef:
        public list collections
        public size_t size
        int* lengths
        int* indices
        public bint done
        public size_t total_combinations

    def __cinit__(self, collections):
        cdef:
            size_t i, n
            list sublist

        self.collections = list(map(list, collections))
        self.size = PyList_GET_SIZE(self.collections)
        self.lengths = <int*>malloc(sizeof(int) * self.size)
        self.indices = <int*>malloc(sizeof(int) * self.size)
        self.done = False
        self.total_combinations = 1

        for i in range(self.size):
            sublist = <list>PyList_GetItem(self.collections, i)
            self.lengths[i] = PyList_GET_SIZE(sublist)
            self.total_combinations *= self.lengths[i]
            if self.lengths[i] == 0:
                self.done = True
            self.indices[i] = 0

    def __dealloc__(self):
        free(self.lengths)
        free(self.indices)

    def get_lengths(self):
        return [self.lengths[i] for i in range(self.size)]

    def get_indices(self):
        return [self.indices[i] for i in range(self.size)]

    cpdef bint has_more(self):
        return not self.done

    cpdef list compose_next_value(self):
        cdef:
            int i, ix
            list result, sublist
            object value
        if self.done:
            return None
        result = []
        for i in range(self.size):
            sublist = <list>PyList_GetItem(self.collections, i)
            ix = self.indices[i]
            value = <object>PyList_GetItem(sublist, ix)
            result.append(value)
        return result

    cpdef advance(self):
        for i in range(self.size - 1, -1, -1):
            if self.indices[i] == self.lengths[i] - 1:
                self.indices[i] = 0
                if i == 0:
                    self.done = True
            else:
                self.indices[i] += 1
                break

    cpdef list get_next_value(self):
        cdef list value
        if self.done:
            return None
        else:
            value = self.compose_next_value()
            self.advance()
            return value

    def __iter__(self):
        return self

    def __next__(self):
        if self.done:
            raise StopIteration()
        value = self.compose_next_value()
        if value is None:
            raise StopIteration()
        else:
            self.advance()
        return value


cdef tuple _conform_envelopes(list experimental, list base_theoretical, size_t* n_missing_out):
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
            fpeak = FittedPeak._create(peak.mz, 1, 1, -1, -1, 0, 1, 0, 0)
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
    n_missing_out[0] = n_missing
    return cleaned_eid, tid


cdef class envelope_conformer:
    cdef:
        list experimental
        list theoretical
        size_t n_missing

    @staticmethod
    cdef envelope_conformer _create():
        return envelope_conformer.__new__(envelope_conformer)

    cdef void acquire(self, list experimental, list theoretical):
        self.experimental = experimental
        self.theoretical = theoretical
        self.n_missing = 0

    cdef void conform(self):
        cdef:
            double total = 0
            size_t n_missing = 0
            size_t i = 0
            list cleaned_eid
            list tid
            FittedPeak fpeak
            TheoreticalPeak tpeak, peak

        cleaned_eid = []
        
        for i in range(PyList_GET_SIZE(self.experimental)):
            fpeak = <FittedPeak>PyList_GET_ITEM(self.experimental, i)
            if fpeak is None:
                peak = <TheoreticalPeak>PyList_GET_ITEM(self.theoretical, i)
                fpeak = FittedPeak._create(peak.mz, 1, 1, -1, -1, 0, 1, 0, 0)
                n_missing += 1
            total += fpeak.intensity
            cleaned_eid.append(fpeak)
            i += 1

        tid = []
        for i in range(PyList_GET_SIZE(self.theoretical)):
            tpeak = <TheoreticalPeak>PyList_GET_ITEM(self.theoretical, i)
            peak = tpeak.clone()
            peak.intensity *= total
            tid.append(peak)
        self.experimental = cleaned_eid
        self.theoretical = tid
        self.n_missing = n_missing


cpdef tuple conform_envelopes(list experimental, list base_theoretical):
    cdef:
        size_t n_missing
        tuple out

    out = _conform_envelopes(experimental, base_theoretical, &n_missing)
    out += (n_missing,)
    return out


cdef int has_no_valid_features(list features):
    cdef:
        size_t i, n, cnt
        object feature

    n = PyList_GET_SIZE(features)
    cnt = 0
    for i in range(n):
        feature = <object>PyList_GET_ITEM(features, i)
        if feature is None:
            cnt += 1
    return cnt == n


cdef class LCMSFeatureProcessorBase(object):
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
            double search_width, hit_width
            LCMSFeature el
            size_t i, n
            bint has_miss
        has_miss = False
        f = self.find_all_features(mz, error_tolerance)
        if not f:
            f = [None]
            has_miss = True
        elif interval is not None:
            overlapped = []
            n = PyList_GET_SIZE(f)
            search_width = interval.get_end_time() - interval.get_start_time()
            if search_width == 0:
                return [None]
            for i in range(n):
                el = <LCMSFeature>PyList_GET_ITEM(f, i)
                if el.overlaps_in_time(interval):
                    hit_width = el.get_end_time() - el.get_start_time()
                    if (hit_width / search_width) < 0.05:
                        continue
                    overlapped.append(el)
            f = overlapped
            if PyList_GET_SIZE(f) == 0:
                f = [None]
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

    cpdef tuple conform_envelopes(self, list experimental, list base_theoretical):
        return conform_envelopes(experimental, base_theoretical)

    cdef tuple _conform_envelopes(self, list experimental, list base_theoretical, size_t* n_missing):
        return _conform_envelopes(experimental, base_theoretical, n_missing)

    cdef double _find_thresholded_score(self, dvec* scores, double percentage):
        cdef:
            size_t i, n, count
            double maximum, val, threshold, acc

        maximum = -INF
        count = 0
        acc = 0
        n = scores.used
        for i in range(n):
            val = scores.v[i]
            if val > maximum:
                maximum = val
        threshold = maximum * percentage
        for i in range(n):
            val = scores.v[i]
            if val > threshold:
                acc += val
                count += 1
        if count == 0:
            return 0
        return acc / count

    cpdef LCMSFeatureSetFit _fit_single_feature_set(self, list features, list base_tid, double error_tolerance, 
                                                    int charge, double charge_carrier=PROTON, int max_missed_peaks=1,
                                                    double threshold_scale=0.3):
        cdef:
            double score, final_score, score_acc, neutral_mass
            envelope_conformer conformer
            dvec* score_vec
            np.ndarray scores_array
            FeatureSetIterator feat_iter
            size_t feat_i, feat_n, n_missing, missing_features, counter

        conformer = envelope_conformer._create()
        # feat_iter = FeatureSetIterator._create_with_threshold(features, base_tid, 0.1)
        feat_iter = FeatureSetIterator._create(features)
        score_acc = 0
        counter = 0
        score_vec = make_double_vector()
        while feat_iter.has_more():
            eid = feat_iter.get_next_value()
            if eid is None:
                continue
            conformer.acquire(eid, base_tid)
            conformer.conform()
            cleaned_eid = conformer.experimental
            tid = conformer.theoretical
            n_missing = conformer.n_missing
            if n_missing > max_missed_peaks:
                continue
            score = self.scorer._evaluate(None, cleaned_eid, tid)
            if npy_isnan(score):
                continue
            score_acc += score
            counter += 1
            double_vector_append(score_vec, score)
        # Compute feature score from score_vec
        final_score = self._find_thresholded_score(score_vec, threshold_scale)
        missing_features = 0
        feat_n = PyList_GET_SIZE(features)
        for feat_i in range(feat_n):
            f = <LCMSFeature>PyList_GET_ITEM(features, feat_i)
            if f is None:
                missing_features += 1
        f = <LCMSFeature>PyList_GET_ITEM(features, 0)
        neutral_mass = calc_neutral_mass(
                f.get_mz(), charge, charge_carrier)
        scores_array = double_vector_to_ndarray(score_vec)
        fit = LCMSFeatureSetFit._create(
            features, base_tid, final_score, charge, missing_features,
            [], None, neutral_mass, counter, scores_array)
        free_double_vector(score_vec)
        return fit

    cpdef list _fit_feature_set(self, double mz, double error_tolerance, int charge, int left_search=1,
                                int right_search=1, double charge_carrier=PROTON, double truncate_after=0.8,
                                int max_missed_peaks=1, double threshold_scale=0.3, LCMSFeature feature=None):
        cdef:
            double score, final_score, score_acc, neutral_mass
            list base_tid, feature_groups, feature_fits, features
            list eid, cleaned_eid, tid
            np.ndarray scores_array
            tuple temp
            size_t combn_size, combn_i, n_missing, missing_features, counter
            size_t feat_i, feat_n
            FeatureSetIterator feat_iter
            LCMSFeature f
            CartesianProductIterator combn_iter
            envelope_conformer conformer
            dvec* score_vec

        base_tid = self.create_theoretical_distribution(
            mz, charge, charge_carrier, truncate_after)
        feature_groups = self.match_theoretical_isotopic_distribution(
            base_tid, error_tolerance, interval=feature)
        feature_fits = []

        conformer = envelope_conformer._create()

        combn_iter = CartesianProductIterator(feature_groups)
        combn_i = 0
        while combn_iter.has_more():
            features = combn_iter.get_next_value()
            if features is None:
                break
            combn_i += 1
            if combn_i % (100000) == 0:
                print("\t%d combinations have been considered for m/z %0.3f at charge %d in %s (%0.2f%%)" % (
                    combn_i, mz, charge, "-" if feature is None else (
                        "%0.2f-%0.2f" % (feature.start_time, feature.end_time)),
                    100. * combn_i / combn_iter.total_combinations))
            if has_no_valid_features(features):
                continue
            # If the monoisotopic feature wasn't actually observed, create a dummy feature
            # since the monoisotopic feature cannot be None
            if features[0] is None:
                features = list(features)
                features[0] = EmptyFeature._create(mz)
            # feat_iter = FeatureSetIterator._create_with_threshold(features, base_tid, 0.1)
            feat_iter = FeatureSetIterator._create(features)
            counter = 0
            score_acc = 0
            score_vec = make_double_vector()
            while feat_iter.has_more():
                eid = feat_iter.get_next_value()
                if eid is None:
                    continue
                counter += 1
                conformer.acquire(eid, base_tid)
                conformer.conform()
                cleaned_eid = conformer.experimental
                tid = conformer.theoretical
                n_missing = conformer.n_missing
                if n_missing > max_missed_peaks:
                    continue

                score = self.scorer._evaluate(None, cleaned_eid, tid)
                if npy_isnan(score):
                    continue
                score_acc += score
                double_vector_append(score_vec, score)
            # Compute feature score from score_vec
            final_score = self._find_thresholded_score(score_vec, threshold_scale)
            missing_features = 0
            feat_n = PyList_GET_SIZE(features)
            for feat_i in range(feat_n):
                f = <LCMSFeature>PyList_GET_ITEM(features, feat_i)
                if f is None:
                    missing_features += 1
            f = <LCMSFeature>PyList_GET_ITEM(features, 0)
            neutral_mass = calc_neutral_mass(
                    f.get_mz(), charge, charge_carrier)
            scores_array = double_vector_to_ndarray(score_vec)
            fit = LCMSFeatureSetFit._create(
                features, base_tid, final_score, charge, missing_features,
                [], None, neutral_mass, counter, scores_array)
            free_double_vector(score_vec)
            if self.scorer.reject_score(fit.score):
                continue
            feature_fits.append(fit)
        return feature_fits
