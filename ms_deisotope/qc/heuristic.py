# -*- coding: utf-8 -*-
"""Heuristic methods for estimating whether a peak list is likely
to produce a high quality interpretation.
"""
import numbers

import numpy as np


def xrea(peaks):
    """Calculate an approximation of the Xrea spectral quality metric.

    Xrea is essentially a measure of whether or not a spectrum is "flat", having no prominent
    peaks. Such spectra are more likely to be low quality.

    This implementation is based upon the implementation in SpectraST's denoising method [1].

    Parameters
    ----------
    peaks : :class:`Iterable` of :class:`PeakLike` or :class:`Iterable` of :class:`float`
        The peaks to calculate the quality measure of, or their intensities

    Returns
    -------
    :class:`float`

    References
    ----------
    [1] Shao, W., & Lam, H. (2013). Denoising Peptide Tandem Mass Spectra for Spectral Libraries: A Bayesian Approach.
        Journal of Proteome Research, 12(7), 3223–3232. https://doi.org/10.1021/pr400080b
    """
    n_peaks = len(peaks)
    if n_peaks == 0:
        return 0
    if isinstance(peaks, np.ndarray):
        intensities = np.array(peaks)
    else:
        try:
            intensities = np.array([p.intensity for p in peaks])
        except AttributeError:
            if isinstance(peaks[0], numbers.Number):
                intensities = np.array(peaks)
            else:
                raise

    intensities[::].sort()

    cumulative_intensity = np.cumsum(intensities)
    normalizer = np.cumsum(cumulative_intensity)[-1]
    cumulative_intensity = cumulative_intensity[-1]

    triangle = (cumulative_intensity * float(n_peaks + 1)) * 0.5
    xrea_ = (triangle - normalizer) / triangle
    return xrea_
