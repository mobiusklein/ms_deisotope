# cython: embedsignature=True

cimport cython
from libc.math cimport fabs, sqrt, log
from libc.stdlib cimport malloc, free
import operator

from cpython.list cimport PyList_New, PyList_GetItem, PyList_Size, PyList_GET_ITEM, PyList_SET_ITEM, PyList_GET_SIZE
from cpython.sequence cimport PySequence_List

from brainpy._c.isotopic_distribution cimport TheoreticalPeak
from ms_peak_picker._c.peak_set cimport PeakSet, FittedPeak


@cython.nonecheck(False)
@cython.cdivision(True)
cdef bint isclose(double x, double y, double rtol=1.e-5, double atol=1.e-8):
    return abs(x-y) <= (atol + rtol * abs(y))



@cython.freelist(500)
cdef class IsotopicFitRecord(object):

    def __cinit__(self, FittedPeak seed_peak, double score, int charge, list theoretical, list experimental):
        self.seed_peak = seed_peak
        self.score = score
        self.charge = charge
        self.experimental = experimental
        self.theoretical = theoretical
        self.monoisotopic_peak = experimental[0]

    def clone(self):
        return self.__class__(self.seed_peak, self.score, self.charge, self.theoretical, self.experimental)

    def __reduce__(self):
        return self.__class__, (self.seed_peak, self.score, self.charge, self.theoretical, self.experimental)

    cpdef bint _eq(self, IsotopicFitRecord other):
        return (self.score == other.score and
                self.charge == other.charge and
                self.experimental == other.experimental and
                self.theoretical == other.theoretical)

    cpdef bint _ne(self, IsotopicFitRecord other):
        return not (self == other)

    cpdef bint _lt(self, IsotopicFitRecord other):
        return self.score < other.score

    cpdef bint _gt(self, IsotopicFitRecord other):
        return self.score > other.score

    def __richcmp__(self, IsotopicFitRecord other, int code):
        if other is None:
            if code == 3:
                return True
            else:
                return False

        if code == 0:
            return self._lt(other)
        elif code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))
        elif code == 4:
            return self._gt(other)


    def __hash__(self):
        return hash((self.monoisotopic_peak.mz, self.charge))

    def __getitem__(self, index):
        if index == 0:
            return self.score
        else:
            raise KeyError(index)

    def __iter__(self):
        yield self.score
        yield self.charge
        yield self.experimental
        yield self.theoretical

    @property
    def npeaks(self):
        return len(self.experimental)

    def __repr__(self):
        return "IsotopicFitRecord(score=%0.5f, charge=%d, npeaks=%d, monoisotopic_mz=%0.5f)" % (
            self.score, self.charge, self.npeaks, self.monoisotopic_peak.mz)


cdef class FitSelectorBase(object):
    _default_minimum_score = 0

    def __init__(self, minimum_score=None):
        if minimum_score is None:
            minimum_score = self._default_minimum_score
        self.minimum_score = minimum_score

    cpdef IsotopicFitRecord best(self, object results):
        return NotImplemented

    def __call__(self, *args, **kwargs):
        return self.best(*args, **kwargs)

    cpdef bint reject(self, IsotopicFitRecord result):
        return NotImplemented


cdef class MinimizeFitSelector(FitSelectorBase):
    cpdef IsotopicFitRecord best(self, object results):
        return min(results, key=operator.attrgetter("score"))

    cpdef bint reject(self, IsotopicFitRecord fit):
        return fit.score > self.minimum_score


cdef class MaximizeFitSelector(FitSelectorBase):
    cpdef IsotopicFitRecord best(self, object results):
        return max(results, key=operator.attrgetter("score"))

    cpdef bint reject(self, IsotopicFitRecord fit):
        return fit.score < self.minimum_score


cdef class IsotopicFitterBase(object):
    _default_selector_type = MinimizeFitSelector

    def __init__(self, selector=None):
        if selector is None:
            selector = self._default_selector_type()
        self.select = selector

    def evaluate(self, PeakIndex peaklist, list observed, list expected, **kwargs):
        return self._evaluate(peaklist, observed, expected)

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected):
        return 0

    def __call__(self, *args, **kwargs):
        return self.evaluate(*args, **kwargs)

    cpdef bint reject(self, IsotopicFitRecord fit):
        return self.select.reject(fit)


cdef double sum_intensity_theoretical(list peaklist):
    cdef:
        double summed
        size_t i
        TheoreticalPeak peak
    summed = 0
    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <TheoreticalPeak>PyList_GET_ITEM(peaklist, i)
        summed += peak.intensity
    return summed    


cdef double sum_intensity_fitted(list peaklist):
    cdef:
        double summed
        size_t i
        FittedPeak peak
    summed = 0
    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <FittedPeak>PyList_GET_ITEM(peaklist, i)
        summed += peak.intensity
    return summed


@cython.cdivision
cdef double* normalize_intensity_theoretical(list peaklist, double* out, double total):
    cdef:
        size_t i
        TheoreticalPeak peak

    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <TheoreticalPeak>PyList_GET_ITEM(peaklist, i)
        out[i] = peak.intensity / total


@cython.cdivision
cdef double* normalize_intensity_fitted(list peaklist, double* out, double total):
    cdef:
        size_t i
        FittedPeak peak

    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <FittedPeak>PyList_GET_ITEM(peaklist, i)
        out[i] = peak.intensity / total


cdef class ScaledGTestFitter(IsotopicFitterBase):
    @cython.cdivision
    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected):
        cdef:
            double total_observed
            double total_expected
            double* normalized_observed
            double* normalized_expected
            double g_score, obs, theo, log_ratio
            size_t n
        total_observed = sum_intensity_fitted(observed)
        total_expected = sum_intensity_theoretical(expected)
        n = PyList_GET_SIZE(observed)

        normalized_observed = <double*>malloc(sizeof(double) * n)
        normalize_intensity_fitted(observed, normalized_observed, total_observed)

        normalized_expected = <double*>malloc(sizeof(double) * n)
        normalize_intensity_theoretical(expected, normalized_expected, total_expected)

        g_score = 0.
        for i in range(n):
            obs = normalized_observed[i]
            theo = normalized_expected[i]

            log_ratio = log(obs / theo)
            g_score += obs * log_ratio

        free(normalized_observed)
        free(normalized_expected)

        return g_score * 2.


cdef ScaledGTestFitter g_test_scaled

g_test_scaled = ScaledGTestFitter()


cdef double max_intensity_theoretical(list peaklist):
    cdef:
        double maximum
        size_t i
        TheoreticalPeak peak
    maximum = 0
    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <TheoreticalPeak>PyList_GET_ITEM(peaklist, i)
        if peak.intensity > maximum:
            maximum = peak.intensity
    return maximum


cdef double max_intensity_fitted(list peaklist):
    cdef:
        double maximum
        size_t i
        FittedPeak peak
    maximum = 0
    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <FittedPeak>PyList_GET_ITEM(peaklist, i)
        if peak.intensity > maximum:
            maximum = peak.intensity
    return maximum


cdef class LeastSquaresFitter(IsotopicFitterBase):
    
    @cython.cdivision
    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected):
        cdef:
            double exp_max, theo_max, sum_of_squared_errors, sum_of_squared_theoreticals
            double normed_theo, normed_expr
            size_t i
            TheoreticalPeak t
            FittedPeak e

        exp_max = max_intensity_fitted(observed)
        theo_max = max_intensity_theoretical(expected)

        sum_of_squared_errors = 0
        sum_of_squared_theoreticals = 0

        assert len(observed) == len(expected)

        # for e, t in zip(observed, expected):
        for i in range(PyList_GET_SIZE(observed)):
            e = <FittedPeak>PyList_GetItem(observed, i)
            t = <TheoreticalPeak>PyList_GetItem(expected, i)
            normed_expr = e.intensity / exp_max
            normed_theo = t.intensity / theo_max
            sum_of_squared_errors += (normed_theo - normed_expr) ** 2
            sum_of_squared_theoreticals += normed_theo ** 2
        return sum_of_squared_errors / sum_of_squared_theoreticals


cdef LeastSquaresFitter least_squares

least_squares = LeastSquaresFitter()


cdef class TopFixingFitterSelector(MaximizeFitSelector):
    cpdef IsotopicFitRecord best(self, object results):
        cdef:
            list ordered_results
            double lower_limit
            list filtered
            double best_score, new_score
            IsotopicFitRecord best_case
            IsotopicFitRecord case
            size_t i

        if len(results) == 0:
            raise ValueError("No options for selection")

        ordered_results = sorted(results, reverse=True)
        lower_limit = ordered_results[0].score * 0.75
        filtered = []

        for i in range(PyList_Size(ordered_results)):
            case = <IsotopicFitRecord>PyList_GetItem(ordered_results, i)
            if case.score > lower_limit:
                filtered.append(case)
            else:
                break

        best_score = self.minimum_score - 0.0000001
        best_case = None


        for i in range(PyList_Size(filtered)):
            case = <IsotopicFitRecord>PyList_GetItem(filtered, i)
            score = g_test_scaled._evaluate(None, case.experimental, case.theoretical)
            new_score = case.score * (1 - score)
            if new_score > best_score:
                best_score = new_score
                # case.score = new_score
                best_case = case

        return best_case


@cython.cdivision
cdef double score_peak(FittedPeak obs, TheoreticalPeak theo, double mass_error_tolerance=0.02, double minimum_signal_to_noise=1) nogil:
    cdef:
        double mass_error, abundance_diff
    if obs.signal_to_noise < minimum_signal_to_noise:
        return 0.

    mass_error = fabs(obs.mz - theo.mz)

    if mass_error <= mass_error_tolerance:
        mass_accuracy = 1 - mass_error / mass_error_tolerance
    else:
        mass_accuracy = 0

    if obs.intensity < theo.intensity and (((theo.intensity - obs.intensity) / obs.intensity) <= 1):
        abundance_diff = 1 - ((theo.intensity - obs.intensity) / obs.intensity)
    elif obs.intensity >= theo.intensity and (((obs.intensity - theo.intensity) / obs.intensity) <= 1):
        abundance_diff = sqrt(1 - ((obs.intensity - theo.intensity) / obs.intensity))
    else:
        abundance_diff = 0.
    score = sqrt(theo.intensity) * mass_accuracy * abundance_diff
    return score

cdef class MSDeconVFitter(IsotopicFitterBase):
    _default_selector_type = TopFixingFitterSelector

    def __init__(self, minimum_score=10):
        self.select = TopFixingFitterSelector()
        self.select.minimum_score = minimum_score

    @cython.cdivision
    cdef double score_peak(self, FittedPeak obs, TheoreticalPeak theo, double mass_error_tolerance=0.02, double minimum_signal_to_noise=1) nogil:
        return score_peak(obs, theo, mass_error_tolerance, minimum_signal_to_noise)

    @cython.cdivision
    cdef double _reweight(self, FittedPeak obs, TheoreticalPeak theo, double scale_obs, double scale_theo) nogil:
        cdef:
            double normed_obs, normed_theo

        normed_obs = obs.intensity / scale_obs
        normed_theo = theo.intensity / scale_theo
        return normed_obs * log(normed_obs / normed_theo)

    cpdef double reweight(self, FittedPeak obs, TheoreticalPeak theo, double scale_obs, double scale_theo):
        return self._reweight(obs, theo, scale_obs, scale_theo)

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected, double mass_error_tolerance=0.02):
        cdef:
            size_t i, n
            FittedPeak obs
            TheoreticalPeak theo
            double score

        n = PyList_GET_SIZE(observed)
        score = 0
        for i in range(n):
            obs = <FittedPeak>PyList_GET_ITEM(observed, i)
            theo = <TheoreticalPeak>PyList_GET_ITEM(expected, i)
            score += self.score_peak(obs, theo, mass_error_tolerance, 1)

        return score


cdef class PenalizedMSDeconVFitter(IsotopicFitterBase):
    def __init__(self, minimum_score=10):
        self.select = MaximizeFitSelector(minimum_score)
        self.msdeconv = MSDeconVFitter()
        self.penalizer = ScaledGTestFitter()

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected, double mass_error_tolerance=0.02):
        cdef:
            double score, penalty
        score = self.msdeconv._evaluate(peaklist, observed, expected, mass_error_tolerance)
        penalty = self.penalizer._evaluate(peaklist, observed, expected)
        return score * (1 - penalty)


cdef class FunctionScorer(IsotopicFitterBase):

    def __init__(self, function, minimum_score=10, selector_type=MaximizeFitSelector):
        self.function = function
        self.select = MaximizeFitSelector(minimum_score)

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected):
        return self.function(observed, expected)
