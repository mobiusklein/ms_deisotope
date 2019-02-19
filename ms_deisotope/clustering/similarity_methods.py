import math
from collections import defaultdict

import numpy as np


def peak_set_similarity(peak_set_a, peak_set_b, precision=0):
    """Computes the normalized dot product, also called cosine similarity between
    two peak sets, a similarity metric ranging between 0 (dissimilar) to 1.0 (similar).

    Parameters
    ----------
    peak_set_a : Iterable of Peak-like
    peak_set_b : Iterable of Peak-like
        The two peak collections to compare. It is usually only useful to
        compare the similarity of peaks of the same class, so the types
        of the elements of `peak_set_a` and `peak_set_b` should match.
    precision : int, optional
        The precision of rounding to use when binning spectra. Defaults to 0

    Returns
    -------
    float
        The similarity between peak_set_a and peak_set_b. Between 0.0 and 1.0
    """
    bin_a = defaultdict(float)
    bin_b = defaultdict(float)

    positions = set()

    for peak in peak_set_a:
        mz = round(peak.mz, precision)
        bin_a[mz] += peak.intensity
        positions.add(mz)

    for peak in peak_set_b:
        mz = round(peak.mz, precision)
        bin_b[mz] += peak.intensity
        positions.add(mz)

    return dot_product(positions, bin_a, bin_b)


def bin_peaks(peak_set_a, peak_set_b, precision):
    bin_a = defaultdict(float)
    bin_b = defaultdict(float)

    positions = set()

    for peak in peak_set_a:
        mz = round(peak.mz, precision)
        bin_a[mz] += peak.intensity
        positions.add(mz)

    for peak in peak_set_b:
        mz = round(peak.mz, precision)
        bin_b[mz] += peak.intensity
        positions.add(mz)
    return positions, bin_a, bin_b


def dot_product(positions, bin_a, bin_b, normalize=True):
    z = 0
    n_a = 0
    n_b = 0

    for mz in positions:
        a = bin_a[mz]
        b = bin_b[mz]
        z += a * b
        n_a += a ** 2
        n_b += b ** 2

    if not normalize:
        return z

    n_ab = math.sqrt(n_a) * math.sqrt(n_b)
    if n_ab == 0.0:
        return 0.0
    else:
        return z / n_ab


sparse_similarity = peak_set_similarity


try:
    _has_c = True
    from ms_deisotope._c import similarity_methods as csimilarity_methods

    def peak_set_similarity(peak_set_a, peak_set_b, precision=0):
        if precision > 2:
            return sparse_similarity(peak_set_a, peak_set_b, precision)
        else:
            return csimilarity_methods.peak_set_similarity(peak_set_a, peak_set_b, precision)

    peak_set_similarity.__doc__ = sparse_similarity.__doc__

except ImportError:
    _has_c = False
