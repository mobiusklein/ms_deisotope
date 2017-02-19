import math
import numpy as np
import operator

from .utils import Base

eps = 1e-4


class IsotopicFitRecord(object):
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
    __slots__ = ["seed_peak", "score", "charge", "experimental", "theoretical",
                 "monoisotopic_peak", "data", "missed_peaks"]

    def __init__(self, seed_peak, score, charge, theoretical, experimental, data=None,
                 missed_peaks=0, **kwargs):
        self.seed_peak = seed_peak
        self.score = score
        self.charge = charge
        self.experimental = experimental
        self.theoretical = theoretical
        self.monoisotopic_peak = experimental[0]
        self.data = data
        self.missed_peaks = missed_peaks

    def clone(self):
        return self.__class__(
            self.seed_peak, self.score, self.charge, self.theoretical, self.experimental, self.data, self.missed_peaks)

    def __reduce__(self):
        return self.__class__, (
            self.seed_peak, self.score, self.charge, self.theoretical, self.experimental, self.data, self.missed_peaks)

    def __eq__(self, other):
        val = (self.score == other.score and
               self.charge == other.charge and
               self.experimental == other.experimental and
               self.theoretical == other.theoretical)
        if self.data is not None or other.data is not None:
            val = val and (self.data == other.data)
        return val

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.score < other.score

    def __gt__(self, other):
        return self.score > other.score

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


class FitSelectorBase(Base):
    """An object that controls the filtering and
    selection of IsotopicFitRecord

    Attributes
    ----------
    minimum_score : int
        The minimum score needed to be a candidate for selection. If the
        FitSelector is `minimizing` it is the maximal score to be a candidate.
    """
    minimum_score = 0

    def __init__(self, minimum_score=0):
        self.minimum_score = minimum_score

    def best(self, results):
        return NotImplemented

    def __call__(self, *args, **kwargs):
        return self.best(*args, **kwargs)

    def reject(self, result):
        return NotImplemented

    def is_maximizing(self):
        return False


class MinimizeFitSelector(FitSelectorBase):
    """A FitSelector which tries to minimize the score of the best fit.
    """
    def best(self, results):
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

    def reject(self, fit):
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

    def is_maximizing(self):
        """Returns that this is *not* a maximizing selector

        Returns
        -------
        False
        """
        return False


class MaximizeFitSelector(FitSelectorBase):
    """A FitSelector which tries to maximize the score of the best fit.
    """
    def best(self, results):
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

    def reject(self, fit):
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

    def is_maximizing(self):
        """Returns that this *is* a maximizing selector

        Returns
        -------
        False
        """
        return True


class IsotopicFitterBase(Base):

    def __init__(self, score_threshold=0.5):
        self.select = MinimizeFitSelector(score_threshold)

    def evaluate(self, peaklist, observed, expected, **kwargs):
        return NotImplemented

    def _evaluate(self, peaklist, observed, expected, **kwargs):
        return self.evaluate(peaklist, observed, expected, **kwargs)

    def __call__(self, *args, **kwargs):
        return self.evaluate(*args, **kwargs)

    def reject(self, fit):
        return self.select.reject(fit)

    def is_maximizing(self):
        return self.select.is_maximizing()

    def configure(self, deconvoluter, **kwargs):
        return self


class GTestFitter(IsotopicFitterBase):

    def evaluate(self, peaklist, observed, expected, **kwargs):
        g_score = 2 * sum([obs.intensity * np.log(
            obs.intensity / theo.intensity) for obs, theo in zip(observed, expected)])
        return g_score


g_test = GTestFitter()


class ScaledGTestFitter(IsotopicFitterBase):

    def evaluate(self, peaklist, observed, expected, **kwargs):
        total_observed = sum(p.intensity for p in observed)
        total_expected = sum(p.intensity for p in expected)
        total_expected += eps
        normalized_observed = [obs.intensity /
                               total_observed for obs in observed]
        normalized_expected = [theo.intensity /
                               total_expected for theo in expected]
        g_score = 2 * sum([obs * np.log(obs / theo) for obs, theo in zip(
            normalized_observed, normalized_expected)])
        return g_score


g_test_scaled = ScaledGTestFitter()


class ChiSquareFitter(IsotopicFitterBase):

    def evaluate(self, peaklist, observed, expected, **kwargs):
        score = sum([(obs.intensity - theo.intensity)**2 / theo.intensity
                     for obs, theo in zip(observed, expected)])
        return score


chi_sqr_test = ChiSquareFitter()


class LeastSquaresFitter(IsotopicFitterBase):

    def evaluate(self, peaklist, observed, expected, **kwargs):
        exp_max = max(p.intensity for p in observed)
        theo_max = max(p.intensity for p in expected)

        sum_of_squared_errors = 0
        sum_of_squared_theoreticals = 0

        for e, t in zip(observed, expected):
            normed_expr = e.intensity / exp_max
            normed_theo = t.intensity / theo_max
            sum_of_squared_errors += (normed_theo - normed_expr) ** 2
            sum_of_squared_theoreticals += normed_theo ** 2
        return sum_of_squared_errors / sum_of_squared_theoreticals


least_squares = LeastSquaresFitter()


class MSDeconVFitter(IsotopicFitterBase):

    def __init__(self, minimum_score=10):
        self.select = MaximizeFitSelector()
        self.select.minimum_score = minimum_score

    def calculate_minimum_signal_to_noise(self, observed):
        snr = 0
        n = 0
        for obs in observed:
            if obs.signal_to_noise < 1:
                continue
            snr += obs.signal_to_noise
            n += 1
        return (snr / n) * 0.05

    def reweight(self, obs, theo, obs_total, theo_total):
        norm_obs = obs.intensity / obs_total
        norm_theo = theo.intensity / theo_total
        return norm_obs * np.log(norm_obs / norm_theo)

    def score_peak(self, obs, theo, mass_error_tolerance=0.02, minimum_signal_to_noise=1):
        if obs.signal_to_noise < minimum_signal_to_noise:
            return 0.

        mass_error = np.abs(obs.mz - theo.mz)

        if mass_error <= mass_error_tolerance:
            mass_accuracy = 1 - mass_error / mass_error_tolerance
        else:
            mass_accuracy = 0

        if obs.intensity < theo.intensity and (((theo.intensity - obs.intensity) / obs.intensity) <= 1):
            abundance_diff = 1 - \
                ((theo.intensity - obs.intensity) / obs.intensity)
        elif obs.intensity >= theo.intensity and (((obs.intensity - theo.intensity) / obs.intensity) <= 1):
            abundance_diff = np.sqrt(
                1 - ((obs.intensity - theo.intensity) / obs.intensity))
        else:
            abundance_diff = 0.
        score = np.sqrt(theo.intensity) * mass_accuracy * abundance_diff
        return score

    def evaluate(self, peaklist, observed, expected, mass_error_tolerance=0.02, **kwargs):
        score = 0
        for obs, theo in zip(observed, expected):
            inc = self.score_peak(obs, theo, mass_error_tolerance, 1)
            score += inc
        return score


class PenalizedMSDeconVFitter(IsotopicFitterBase):

    def __init__(self, minimum_score=10, penalty_factor=1.):
        self.select = MaximizeFitSelector(minimum_score)
        self.msdeconv = MSDeconVFitter()
        self.penalizer = ScaledGTestFitter()
        self.penalty_factor = penalty_factor

    def evaluate(self, peaklist, observed, expected, mass_error_tolerance=0.02, **kwargs):
        score = self.msdeconv.evaluate(
            observed, expected, mass_error_tolerance)
        penalty = abs(self.penalizer.evaluate(observed, expected))
        return score * (1 - penalty * self.penalty_factor)


def decon2ls_chisqr_test(peaklist, observed, expected, **kwargs):
    fit_total = 0
    sum_total = 0
    for obs, theo in zip(observed, expected):
        intensity_diff = obs.intensity - theo.intensity
        fit_total += (intensity_diff ** 2) / (theo.intensity + obs.intensity)
        sum_total += theo.intensity * obs.intensity
    return fit_total / (sum_total + 0.01)


class InterferenceDetection(object):

    def __init__(self, peaklist):
        self.peaklist = peaklist

    def detect_interference(self, experimental_peaks):
        min_peak = experimental_peaks[0]
        max_peak = experimental_peaks[-1]

        region = self.peaklist.between(
            min_peak.mz - min_peak.full_width_at_half_max,
            max_peak.mz + max_peak.full_width_at_half_max)

        included_intensity = sum(p.intensity for p in experimental_peaks)
        region_intensity = sum(p.intensity for p in region)

        score = 1 - (included_intensity / region_intensity)
        return score


class DistinctPatternFitter(IsotopicFitterBase):

    def __init__(self, minimum_score=0.3, peak_count_scale=1.5, domain_scale=100.):
        self.select = MinimizeFitSelector(minimum_score)
        self.interference_detector = None
        self.g_test_scaled = ScaledGTestFitter()
        self.peak_count_scale = peak_count_scale
        self.domain_scale = domain_scale

    def evaluate(self, peaklist, experimental, theoretical):
        npeaks = float(len(experimental))
        if self.interference_detector is None:
            self.interference_detector = InterferenceDetection(peaklist)

        score = self.g_test_scaled(peaklist, experimental, theoretical)
        score *= abs((self.interference_detector.detect_interference(experimental) + 0.00001) / (
            npeaks * self.peak_count_scale)) * self.domain_scale
        return score


def percentile(N, percent):
    if not N:
        return None
    k = (len(N) - 1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return N[int(k)]
    d0 = N[int(f)] * (c - k)
    d1 = N[int(c)] * (k - f)
    return d0 + d1


try:
    _c = True
    _IsotopicFitRecord = IsotopicFitRecord
    _LeastSquaresFitter = LeastSquaresFitter
    _MSDeconVFitter = MSDeconVFitter
    _ScaledGTestFitter = ScaledGTestFitter
    _PenalizedMSDeconVFitter = PenalizedMSDeconVFitter
    _DistinctPatternFitter = DistinctPatternFitter

    from ._c.scoring import (
        IsotopicFitRecord, LeastSquaresFitter, MSDeconVFitter,
        ScaledGTestFitter, PenalizedMSDeconVFitter, DistinctPatternFitter,
        ScaledPenalizedMSDeconvFitter)
except ImportError as e:
    print(e)
    _c = False

msdeconv = MSDeconVFitter()
least_squares = LeastSquaresFitter()
g_test_scaled = ScaledGTestFitter()
penalized_msdeconv = PenalizedMSDeconVFitter()
distinct_pattern_fitter = DistinctPatternFitter()


class MassShiftSupportPostProcessorBase(object):

    def __init__(self, shifts=None):
        if shifts is None:
            shifts = []
        self.shifts = shifts
        self._fits = None

    @property
    def fits(self):
        return self._fits

    @fits.setter
    def fits(self, value):
        self._fits = value

    def reweight(self, fit, shift):
        pass
