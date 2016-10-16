# cython: embedsignature=True

cimport cython
from libc.math cimport fabs, sqrt, log, ceil, floor
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
    """Describes a single isotopic pattern fit, comparing how well an
    experimentally observed sequence of peaks matches a theoretical isotopic
    pattern.
    
    IsotopicFitRecord instances are hashable and orderable (by score).

    Attributes
    ----------
    charge : int
        The charge state used to generate the theoretical pattern
    data : object
        An arbitrary Python object containing extra information
    experimental : list of FittedPeak
        The observed experimental peaks to be fitted
    missed_peaks : int
        The number of peaks in the theoretical pattern that do not
        have a matching experimental peak
    monoisotopic_peak : FittedPeak
        The fitted peak which corresponds to the monoisotopic peak
    score : float
        The score assigned to the fit by an IsotopicFitter object
    seed_peak : FittedPeak
        The peak that was used to initiate the fit. This may be unused
        if not using an Averagine method
    theoretical : list of TheoreticalPeak
        The theoretical isotopic pattern being fitted on the experimental data
    """
    def __cinit__(self, FittedPeak seed_peak, double score, int charge, list theoretical, list experimental,
                  object data=None, int missed_peaks=0):
        self.seed_peak = seed_peak
        self.score = score
        self.charge = charge
        self.experimental = experimental
        self.theoretical = theoretical
        self.monoisotopic_peak = experimental[0]
        self.data = data
        self.missed_peaks = missed_peaks

    def clone(self):
        return self.__class__(self.seed_peak, self.score, self.charge, self.theoretical,
                              self.experimental, self.data, self.missed_peaks)

    def __reduce__(self):
        return self.__class__, (self.seed_peak, self.score, self.charge, self.theoretical,
                                self.experimental, self.data, self.missed_peaks)

    cpdef bint _eq(self, IsotopicFitRecord other):
        cdef bint val
        val = (self.score == other.score and
               self.charge == other.charge and
               self.experimental == other.experimental and
               self.theoretical == other.theoretical)
        if self.data is not None or other.data is not None:
            val = val & (self.data == other.data)
        return val

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
    """An object that controls the filtering and
    selection of IsotopicFitRecord

    Attributes
    ----------
    minimum_score : int
        The minimum score needed to be a candidate for selection. If the
        FitSelector is `minimizing` it is the maximal score to be a candidate.
    """
    def __init__(self, minimum_score=0):
        self.minimum_score = minimum_score

    def __reduce__(self):
        return self.__class__, (self.minimum_score,)

    cpdef IsotopicFitRecord best(self, object results):
        return NotImplemented

    def __call__(self, *args, **kwargs):
        return self.best(*args, **kwargs)

    cpdef bint reject(self, IsotopicFitRecord result):
        return NotImplemented

    cpdef bint is_maximizing(self):
        return False


cdef class MinimizeFitSelector(FitSelectorBase):
    """A FitSelector which tries to minimize the score of the best fit.
    """
    cpdef IsotopicFitRecord best(self, object results):
        """Returns the IsotopicFitRecord with the smallest score
        
        Parameters
        ----------
        results : list of IsotopicFitRecord
            List of isotopic fits to select the most optimal case from
        
        Returns
        -------
        IsotopicFitRecord
            The most optimal fit
        """
        return min(results, key=operator.attrgetter("score"))

    cpdef bint reject(self, IsotopicFitRecord fit):
        """Decide whether the fit should be discarded for having too
        large a score. Compares against :attr:`minimum_score`
        
        Parameters
        ----------
        fit : IsotopicFitRecord
        
        Returns
        -------
        bool
        """
        return fit.score > self.minimum_score

    cpdef bint is_maximizing(self):
        """Returns that this is *not* a maximizing selector

        Returns
        -------
        False
        """
        return False


cdef class MaximizeFitSelector(FitSelectorBase):
    """A FitSelector which tries to maximize the score of the best fit.
    """
    cpdef IsotopicFitRecord best(self, object results):
        """Returns the IsotopicFitRecord with the largest score
        
        Parameters
        ----------
        results : list of IsotopicFitRecord
            List of isotopic fits to select the most optimal case from
        
        Returns
        -------
        IsotopicFitRecord
            The most optimal fit
        """
        return max(results, key=operator.attrgetter("score"))

    cpdef bint reject(self, IsotopicFitRecord fit):
        """Decide whether the fit should be discarded for having too
        small a score. Compares against :attr:`minimum_score`
        
        Parameters
        ----------
        fit : IsotopicFitRecord
        
        Returns
        -------
        bool
        """
        return fit.score < self.minimum_score

    cpdef bint is_maximizing(self):
        """Returns that this *is* a maximizing selector

        Returns
        -------
        False
        """
        return True


cdef class IsotopicFitterBase(object):

    def __init__(self, score_threshold=0.5):
        self.select = MinimizeFitSelector(score_threshold)

    def __reduce__(self):
        return self.__class__, (0,), self.__getstate__()

    def __getstate__(self):
        return (self.select,)

    def __setstate__(self, state):
        self.select, = state

    def evaluate(self, PeakIndex peaklist, list observed, list expected, **kwargs):
        return self._evaluate(peaklist, observed, expected)

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected):
        return 0

    def __call__(self, *args, **kwargs):
        return self.evaluate(*args, **kwargs)

    cpdef bint reject(self, IsotopicFitRecord fit):
        return self.select.reject(fit)

    cpdef bint is_maximizing(self):
        return self.select.is_maximizing()

    cpdef IsotopicFitterBase _configure(self, DeconvoluterBase deconvoluter, dict kwargs):
        return self

    def configure(self, DeconvoluterBase deconvoluter, **kwargs):
        return self._configure(deconvoluter, kwargs)


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

    def __init__(self, minimum_score=10):
        self.select = MaximizeFitSelector()
        self.select.minimum_score = minimum_score

    @cython.cdivision
    cdef double score_peak(self, FittedPeak obs, TheoreticalPeak theo, double mass_error_tolerance=0.02, double minimum_signal_to_noise=1) nogil:
        return score_peak(obs, theo, mass_error_tolerance, minimum_signal_to_noise)

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
    def __init__(self, minimum_score=10, penalty_factor=1):
        self.select = MaximizeFitSelector(minimum_score)
        self.msdeconv = MSDeconVFitter()
        self.penalizer = ScaledGTestFitter()
        self.penalty_factor = penalty_factor

    def __reduce__(self):
        return self.__class__, (0,), self.__getstate__()

    def __getstate__(self):
        return (self.select, self.msdeconv, self.penalizer, self.penalty_factor)

    def __setstate__(self, state):
        self.select, self.msdeconv, self.penalizer, self.penalty_factor = state

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected, double mass_error_tolerance=0.02):
        cdef:
            double score, penalty
        score = self.msdeconv._evaluate(peaklist, observed, expected, mass_error_tolerance)
        penalty = self.penalizer._evaluate(peaklist, observed, expected)
        return score * ((1 - penalty * self.penalty_factor))


cdef class FunctionScorer(IsotopicFitterBase):

    def __init__(self, function, minimum_score=10, selector_type=MaximizeFitSelector):
        self.function = function
        self.select = selector_type(minimum_score)

    def __getstate__(self):
        return self.select, self.function

    def __setstate__(self, state):
        self.select, self.function = state

    def __reduce__(self):
        return self.__class__, (None,), self.__getstate__()

    cpdef double _evaluate(self, PeakIndex peaklist, list observed, list expected):
        return self.function(observed, expected)


cdef class InterferenceDetection(object):
    def __init__(self, PeakIndex peaklist):
        self.peaklist = peaklist

    def __reduce__(self):
        return self.__class__, (None,), self.__getstate__()

    def __getstate__(self):
        return (self.peaklist,)

    def __setstate__(self, state):
        self.peaklist, = state

    cdef double detect_interference(self, list experimental_peaks):
        cdef:
            FittedPeak min_peak, max_peak
            PeakSet region
            double included_intensity, region_intensity, score
            size_t n

        n = PyList_GET_SIZE(experimental_peaks)
        if n == 0:
            return 1.0

        min_peak = <FittedPeak>PyList_GET_ITEM(experimental_peaks, 0)
        max_peak = <FittedPeak>PyList_GET_ITEM(experimental_peaks, PyList_GET_SIZE(experimental_peaks) - 1)

        region = self.peaklist._between(
            min_peak.mz - min_peak.full_width_at_half_max,
            max_peak.mz + max_peak.full_width_at_half_max)

        included_intensity = sum_intensity_fitted(experimental_peaks)
        region_intensity = sum_intensity_fitted(list(region))
        if region_intensity == 0:
            return 1.0

        score = 1 - (included_intensity / region_intensity)
        return score


cdef class DistinctPatternFitter(IsotopicFitterBase):

    def __init__(self, minimum_score=0.3, peak_count_scale=1.5, domain_scale=100.):
        self.select = MinimizeFitSelector(minimum_score)
        self.interference_detector = None
        self.g_test_scaled = ScaledGTestFitter()
        self.peak_count_scale = peak_count_scale
        self.domain_scale = domain_scale

    def __getstate__(self):
        return self.select, self.interference_detector, self.g_test_scaled, self.peak_count_scale, self.domain_scale

    def __setstate__(self, state):
        self.select, self.interference_detector, self.g_test_scaled, self.peak_count_scale, self.domain_scale = state

    cpdef IsotopicFitterBase _configure(self, DeconvoluterBase deconvoluter, dict kwargs):
        self.interference_detector = InterferenceDetection(deconvoluter.peaklist)
        return self    

    cpdef double _evaluate(self, PeakIndex peaklist, list experimental, list theoretical):
        cdef:
            double score
            double npeaks

        npeaks = PyList_GET_SIZE(experimental)

        if self.interference_detector is None:
            self.interference_detector = InterferenceDetection(peaklist)

        score = self.g_test_scaled._evaluate(peaklist, experimental, theoretical)
        score *= abs((self.interference_detector.detect_interference(experimental) + 0.00001) / (
            npeaks * self.peak_count_scale)) * self.domain_scale
        return score


cdef double percentile(double[:] N, double percent):
    cdef:
        double k, f, c, d0, d1
    k = (N.shape[0] - 1) * percent
    f = floor(k)
    c = ceil(k)
    if f == c:
        return N[<size_t>(k)]
    d0 = N[<size_t>(f)] * (c - k)
    d1 = N[<size_t>(c)] * (k - f)
    return d0 + d1


cdef class ScaledPenalizedMSDeconvFitter(IsotopicFitterBase):

    def __init__(self,  minimum_score=0.3, penalty_factor=1.):
        self.select = MaximizeFitSelector(minimum_score)
        self.scale_factor = 0.
        self.scorer = PenalizedMSDeconVFitter(penalty_factor=penalty_factor)

    def __getstate__(self):
        return self.select, self.scale_factor, self.scorer

    def __setstate__(self, state):
        self.select, self.scale_factor, self.scorer = state

    @property
    def penalty_factor(self):
        return self.scorer.penalty_factor

    @penalty_factor.setter
    def penalty_factor(self, value):
        self.scorer.penalty_factor = value

    cpdef double _calculate_scale_factor(self, PeakIndex peaklist):
        cdef:
            size_t i
            FittedPeak peak
            double maximum
            double intensity
            PeakSet peaks

        maximum = 0.
        peaks = peaklist.peaks
        for i in range(len(peaks)):
            peak = peaks.getitem(i)
            if peak.intensity > maximum:
                maximum = peak.intensity
        return maximum

    cdef void scale_fitted_peaks(self, list experimental, double factor):
        cdef:
            size_t i
            FittedPeak peak
        for i in range(len(experimental)):
            peak = <FittedPeak>PyList_GET_ITEM(experimental, i)
            peak.intensity *= factor        

    cdef void scale_theoretical_peaks(self, list theoretical, double factor):
        cdef:
            size_t i
            TheoreticalPeak peak
        for i in range(len(theoretical)):
            peak = <TheoreticalPeak>PyList_GET_ITEM(theoretical, i)
            peak.intensity *= factor

    cpdef double _evaluate(self, PeakIndex peaklist, list experimental, list theoretical, double mass_error_tolerance=0.02):
        cdef:
            double score
        if self.scale_factor < 1:
            self.scale_factor = self._calculate_scale_factor(peaklist)
        self.scale_fitted_peaks(experimental, 1. / self.scale_factor)
        self.scale_theoretical_peaks(theoretical, 1. / self.scale_factor)
        score = self.scorer._evaluate(peaklist, experimental, theoretical, mass_error_tolerance)
        self.scale_fitted_peaks(experimental, self.scale_factor)
        self.scale_theoretical_peaks(theoretical, self.scale_factor)
        return score
