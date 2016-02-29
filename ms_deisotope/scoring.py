import numpy as np


def g_test(observed, expected, **kwargs):
    if len(observed) < 2:
        return float("inf")
    g_score = 2 * sum([obs.intensity * np.log(
        obs.intensity / theo.intensity) for obs, theo in zip(observed, expected)])
    return g_score


def g_test_scaled(observed, expected, **kwargs):
    if len(observed) < 2:
        return float("inf")
    total_observed = sum(p.intensity for p in observed)
    total_expected = sum(p.intensity for p in expected)
    normalized_observed = [obs.intensity/total_observed for obs in observed]
    normalized_expected = [theo.intensity/total_expected for theo in expected]
    g_score = 2 * sum([obs * np.log(obs / theo) for obs, theo in zip(
        normalized_observed, normalized_expected)])
    return g_score


def chisqr_test(observed, expected, **kwargs):
    if len(observed) < 2:
        return float("inf")
    score = sum([(obs.intensity - theo.intensity)**2 / theo.intensity
                 for obs, theo in zip(observed, expected)])
    return score


def decon2ls_chisqr_test(observed, expected, **kwargs):
    fit_total = 0
    sum_total = 0
    for obs, theo in zip(observed, expected):
        intensity_diff = obs.intensity - theo.intensity
        fit_total += (intensity_diff ** 2) / (theo.intensity + obs.intensity)
        sum_total += theo.intensity * obs.intensity
    return fit_total / (sum_total + 0.01)


def msdeconv(observed, expected, error_tolerance=1e-5, **kwargs):
    score = 0
    for obs, theo in zip(observed, expected):
        mass_accuracy = 1 - (obs.mz - theo.mz) / error_tolerance

        if obs.intensity > theo.intensity and ((theo.intensity - obs.intensity) / obs.intensity) <= 1:
            abundance_diff = 1 - ((theo.intensity - obs.intensity) / obs.intensity)
        elif theo.intensity >= obs.intensity and ((obs.intensity - theo.intensity) / obs.intensity) <= 1:
            abundance_diff = np.sqrt(1 - ((obs.intensity - theo.intensity) / obs.intensity))
        else:
            abundance_diff = 0
        score += np.sqrt(theo.intensity) * mass_accuracy * abundance_diff
    return score
