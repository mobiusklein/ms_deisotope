import numpy as np
import operator

from .utils import Base

eps = 1e-4


class IsotopicFitRecord(object):

    __slots__ = ["seed_peak", "score", "charge", "experimental", "theoretical", "monoisotopic_peak"]

    def __init__(self, seed_peak, score, charge, theoretical, experimental, **kwargs):
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

    def __eq__(self, other):
        return (self.score == other.score and
                self.charge == other.charge and
                self.experimental == other.experimental and
                self.theoretical == other.theoretical)

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.score < other.score

    def __gt__(self, other):
        return self.score > other.score

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


class FitSelectorBase(Base):
    minimum_score = 0

    def __init__(self, minimum_score=0):
        self.minimum_score = minimum_score

    def best(self, results):
        return NotImplemented

    def __call__(self, *args, **kwargs):
        return self.best(*args, **kwargs)

    def reject(self, result):
        return NotImplemented


class MinimizeFitSelector(FitSelectorBase):
    def best(self, results):
        return min(results, key=operator.attrgetter("score"))

    def reject(self, fit):
        return fit.score > self.minimum_score


class MaximizeFitSelector(FitSelectorBase):
    def best(self, results):
        return max(results, key=operator.attrgetter("score"))

    def reject(self, fit):
        return fit.score < self.minimum_score


class IsotopicFitterBase(Base):
    select = MinimizeFitSelector()

    def evaluate(self, peaklist, observed, expected, **kwargs):
        return NotImplemented

    def _evaluate(self, peaklist, observed, expected, **kwargs):
        return self.evaluate(peaklist, observed, expected, **kwargs)

    def __call__(self, *args, **kwargs):
        return self.evaluate(*args, **kwargs)

    def reject(self, fit):
        return self.select.reject(fit)


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
        normalized_observed = [obs.intensity/total_observed for obs in observed]
        normalized_expected = [theo.intensity/total_expected for theo in expected]
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


class TopFixingFitterSelector(MaximizeFitSelector):
    def best(self, results):
        if len(results) == 0:
            raise ValueError("No options for selection")
        results = sorted(results, reverse=True)
        lower_limit = results[0].score * 0.75
        filtered = [x for x in results if x.score > lower_limit]

        best_score = self.minimum_score - 0.0000001
        best_case = None

        for case in filtered:
            score = g_test_scaled(case.experimental, case.theoretical)
            new_score = case.score * (1 - score)
            if new_score > best_score:
                best_score = new_score
                # case.score = new_score
                best_case = case
        return best_case


class MSDeconVFitter(IsotopicFitterBase):
    select_type = TopFixingFitterSelector
    select = TopFixingFitterSelector()
    select.minimum_score = 10

    def __init__(self, minimum_score=10):
        self.select = self.select_type()
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
            abundance_diff = 1 - ((theo.intensity - obs.intensity) / obs.intensity)
        elif obs.intensity >= theo.intensity and (((obs.intensity - theo.intensity) / obs.intensity) <= 1):
            abundance_diff = np.sqrt(1 - ((obs.intensity - theo.intensity) / obs.intensity))
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
    def __init__(self, minimum_score=10):
        self.select = MaximizeFitSelector(minimum_score)
        self.msdeconv = MSDeconVFitter()
        self.penalizer = ScaledGTestFitter()

    def evaluate(self, peaklist, observed, expected, mass_error_tolerance=0.02, **kwargs):
        score = self.msdeconv.evaluate(observed, expected, mass_error_tolerance)
        penalty = self.penalizer.evaluate(observed, expected)
        return score * (1 - penalty)


def decon2ls_chisqr_test(peaklist, observed, expected, **kwargs):
    fit_total = 0
    sum_total = 0
    for obs, theo in zip(observed, expected):
        intensity_diff = obs.intensity - theo.intensity
        fit_total += (intensity_diff ** 2) / (theo.intensity + obs.intensity)
        sum_total += theo.intensity * obs.intensity
    return fit_total / (sum_total + 0.01)


try:
    _IsotopicFitRecord = IsotopicFitRecord
    _LeastSquaresFitter = LeastSquaresFitter
    _MSDeconVFitter = MSDeconVFitter
    _ScaledGTestFitter = ScaledGTestFitter
    _PenalizedMSDeconVFitter = PenalizedMSDeconVFitter
    from ._c.scoring import (
        IsotopicFitRecord, LeastSquaresFitter, TopFixingFitterSelector, MSDeconVFitter,
        ScaledGTestFitter, PenalizedMSDeconVFitter)
except ImportError:
    pass

msdeconv = MSDeconVFitter()
least_squares = LeastSquaresFitter()
g_test_scaled = ScaledGTestFitter()
penalized_msdeconv = PenalizedMSDeconVFitter()
