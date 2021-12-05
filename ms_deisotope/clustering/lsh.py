'''Implementation based upon [1].

.. warning::

    This implementation is *experimental* and uses a lot of RAM.

References
----------
[1] Kretizberg, P. A., Bern, M., Shu, Q., Yang, F., & Serang, O. (2019). The alphabet projection of spectra.
    Journal of Proteome Research, acs.jproteome.9b00216. https://doi.org/10.1021/acs.jproteome.9b00216
'''

# -*- coding: utf-8 -*-


import math

from collections import defaultdict

import numpy as np
from scipy import signal

from ms_deisotope.utils import uid



def _enumerate_powers_of_2(max_value):
    max_value = int(max_value)
    bins = np.zeros(max_value)
    i = 0
    j = 0
    j2 = 2 ** j
    while i < max_value:
        while i < j2 and i < max_value:
            bins[i] = j2
            i += 1
        j += 1
        j2 = 2 ** j
    return bins


precalc_powers_of_2 = _enumerate_powers_of_2(1000)


def _search_for_nearest_power_of_2(value):
    i = 0
    while value > 2 ** i:
        i += 1
    return 2 ** i


def power_of_2(value):
    if value < len(precalc_powers_of_2):
        return precalc_powers_of_2[value]
    else:
        return _search_for_nearest_power_of_2(value)


def random_unit_vectors(n_planes, n_dims, random_state=None):
    """Generate `n_planes` random vectors for locality sensitive hashing
    of binned spectra with `n_dims` bins.

    Parameters
    ----------
    n_planes : int
        The number of vectors to generate
    n_dims : int
        The size of each random vector

    Returns
    -------
    list
    """
    if random_state is None:
        random_state = np.random
    unit_vectors = []
    for _i in range(n_planes):
        plane = random_state.uniform(-1, 1, n_dims)
        magnitude = np.linalg.norm(plane)
        plane /= magnitude
        unit_vectors.append(plane)
    return unit_vectors



class LSH(object):
    def __init__(self, n_planes, span=2000.0, epsilon=0.02, random_state=None):
        if random_state is None:
            random_state = np.random
        self.span = span
        self.epsilon = epsilon
        self.n_dims = power_of_2(int(span / epsilon))
        self.n_planes = n_planes
        self.random_state = random_state
        self.unit_vectors = self._create_cutting_planes()
        self.hash_bins = defaultdict(list)

    def bin_peaks(self, peaks):
        return BinnedPeakSet(peaks, self.span, self.epsilon, key=peaks)

    def _create_cutting_planes(self):
        return random_unit_vectors(self.n_planes, self.n_dims, self.random_state)

    def add(self, peaks):
        if not isinstance(peaks, BinnedPeakSet):
            peaks = self.bin_peaks(peaks)
        hv = peaks.hash_unit_vectors(self.unit_vectors)
        self.hash_bins[hv.sum()].append(peaks)
        return hv


class BinnedPeakSet(object):
    def __init__(self, peak_set, span=None, epsilon=0.02, key=None):
        if span is None:
            span = 0
        if key is None:
            key = uid()
        self.key = key
        self.peak_set = peak_set
        self.span = max(span, self._guess_span())
        self.epsilon = epsilon

        self.indices = None
        self.binned_signal = None
        self.discretize()

    def _guess_span(self):
        return self.peak_set[-1].mz - self.peak_set[0].mz

    def __eq__(self, other):
        if other is None:
            return False
        try:
            return self.key == other.key
        except AttributeError:
            return self.peak_set == list(other)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.key)

    def bin_steps(self):
        min_mz = self.peak_set[0].mz
        max_mz = self.peak_set[-1].mz
        bins_per_mz = len(self.binned_signal) / (max_mz - min_mz)
        step = int(math.ceil((max_mz - min_mz) / 100.0) * 10)
        locs = []
        labels = []
        for i in range(10):
            locs.append(bins_per_mz * i * step)
            labels.append(min_mz + i * step)
        return locs, labels

    def discretize(self, span=None):
        if span is None:
            span = self.span
        else:
            span = self.span = max(span, self._guess_span())
        num_bins = power_of_2(int(math.ceil(span / self.epsilon)))
        bins = np.linspace(0, span + self.epsilon, num=num_bins)
        min_mz = self.peak_set[0].mz
        # This ignores the possibility that a peak delta falls outside of span
        self.indices = np.digitize([peak.mz - min_mz for peak in self.peak_set], bins)
        self.binned_signal = np.zeros(num_bins)
        for i, peak in enumerate(self.peak_set):
            # The original implementation just set the bin to exp(intensity) of the last
            # peak to fall in it, but exp() overflows with numbers higher than ~700.
            #
            # The absolute magnitude of the value here doesn't matter, just that it is
            # non-zero?
            self.binned_signal[self.indices[i]] += peak.intensity

    def similarity(self, other):
        self_clamped = self.binned_signal.clip(0, 1)
        other_clamped = other.binned_signal.clip(0, 1)
        convolved = signal.fftconvolve(self_clamped, other_clamped[::-1])
        i = np.argmax(convolved)
        score = convolved[i]
        return np.roll(other.binned_signal, i + 1), i + 1, score

    def hash_unit_vectors(self, unit_vectors):
        power_spectrum = np.abs(np.fft.fft(self.binned_signal))
        power_spectrum[0] = 0.0
        value = np.zeros_like(power_spectrum)
        for i, v in enumerate(unit_vectors):
            bit = v.dot(power_spectrum)
            value[i] = (bit > 0) * 2 ** i
        return value


def discretize(peak_set, span=None, epsilon=0.02):
    return BinnedPeakSet(peak_set, span, epsilon)
