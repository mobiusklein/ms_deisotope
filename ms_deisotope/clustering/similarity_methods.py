'''A collection of methods for comparing two peak sets for similarity.

Provides an adaptive cosine similarity implementation for both sparse and
dense spectra.
'''

import math
from collections import defaultdict


def top_n_filter(peak_set, n=40):
    """Keep only the top `n` most abundant peaks in the spectrum.

    Parameters
    ----------
    peak_set : :class:`Iterable` of :class:`~.PeakBase`
        The peaks to filter
    n : int, optional
        The maximum number of peaks to retain, by default 40

    Returns
    -------
    list
    """
    reduced_peaks = sorted(peak_set, key=lambda x: x.intensity, reverse=True)[:n]
    reduced_peaks.sort(key=lambda x: x.mz)
    return reduced_peaks


def sparse_peak_set_similarity(peak_set_a, peak_set_b, precision=0):
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

    return bin_dot_product(positions, bin_a, bin_b)


def bin_dot_product(positions, bin_a, bin_b, normalize=True):
    '''Compute a normalzied dot product between two aligned intensity maps
    '''
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


sparse_similarity = sparse_peak_set_similarity


def ppm_peak_set_similarity(peak_set_a, peak_set_b, error_tolerance=2e-5):
    ab = convolve_peak_sets(peak_set_a, peak_set_b, error_tolerance)
    aa = convolve_peak_sets(peak_set_a, peak_set_a, error_tolerance)
    bb = convolve_peak_sets(peak_set_b, peak_set_b, error_tolerance)
    return convolved_product(ab) / (math.sqrt(convolved_product(aa)) * math.sqrt(convolved_product(bb)))


def convolve_peak_sets(query_peaks, reference_peaks, error_tolerance=2e-5):
    peak_pairs = []
    for peak in query_peaks:
        pairs_for_peak = []
        low = peak.mz - peak.mz * error_tolerance
        high = peak.mz + peak.mz * error_tolerance
        for ref_peak in reference_peaks.between(low, high):
            pairs_for_peak.append((peak, ref_peak))
        peak_pairs.append(pairs_for_peak)
    return peak_pairs


def convolved_product(pairs):
    d = 0.0
    for pair_set in pairs:
        for p1, p2 in pair_set:
            d += p1.intensity * p2.intensity
    return d


try:
    _has_c = True
    from ms_deisotope._c import similarity_methods as csimilarity_methods

    def peak_set_similarity(peak_set_a, peak_set_b, precision=0):
        '''A thin dispatching wrapper for peak_set_similarity methods
        '''
        if peak_set_a is None or peak_set_b is None:
            raise TypeError("Peak sets cannot be None!")
        if precision > 2:
            return sparse_similarity(peak_set_a, peak_set_b, precision)
        else:
            return csimilarity_methods.peak_set_similarity(peak_set_a, peak_set_b, precision)

    ppm_peak_set_similarity = csimilarity_methods.ppm_peak_set_similarity
    align_peaks = csimilarity_methods.align_peaks
    SpectrumAlignment = csimilarity_methods.SpectrumAlignment

    peak_set_similarity.__doc__ = sparse_similarity.__doc__

except ImportError:
    _has_c = False
    peak_set_similarity = sparse_peak_set_similarity
