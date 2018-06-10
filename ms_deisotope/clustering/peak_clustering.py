# pragma: no cover

import numpy as np
from ms_peak_picker import FittedPeak, PeakIndex, PeakSet
from ms_deisotope.peak_set import DeconvolutedPeak, DeconvolutedPeakSet

from collections import Counter, defaultdict


def unspool_peaks(peak_sets):
    return [p for peak_set in peak_sets for p in peak_set]


class PeakScanTree(object):
    def __init__(self, store=None):
        if store is None:
            store = dict()
        self.store = defaultdict(list, store)
        self._mz = None

    def add(self, peak, scan_id):
        self.store[scan_id].append(peak)
        if self._mz is None:
            self._mz = peak.mz

    def extend(self, tree):
        for key, peaks in tree.store.items():
            for peak in peaks:
                self.add(peak, key)

    def __iter__(self):
        for peaks in self.store.values():
            for peak in peaks:
                yield peak

    @property
    def mz(self):
        return self._mz

    def remove_scan(self, scan_id):
        self.store.pop(scan_id)
        self._reset()

    def _reset(self):
        vals = self.store.values()
        if vals:
            peak = vals[0][0]
            self._mz = peak.mz
        else:
            self._mz = None


class PeakCluster(object):
    def __init__(self, peaks=None):
        if peaks is None:
            peaks = []
        self.peaks = peaks
        self.mz = None
        self.probability = 1

    def add(self, peak):
        if isinstance(peak, PeakCluster):
            self.peaks.extend(peak)
        else:
            self.peaks.append(peak)
        if self.mz is None:
            self.mz = self.peaks[0].mz

    def __repr__(self):
        if len(self) == 0:
            return "PeakCluster()"
        return "PeakCluster(%f, %d)" % (self.mz, len(self.peaks))

    def __iter__(self):
        return iter(self.peaks)

    def __len__(self):
        return len(self.peaks)

    def weighted_mz(self):
        prod = 0
        total_intensity = 0
        for peak in self:
            prod += peak.mz * peak.intensity
            total_intensity += peak.intensity
        return (prod / total_intensity)

    def _scaling_factor(self):
        return 0.95 + 0.05 * (1 + self.probability) ** 5

    @property
    def intensity(self):
        return sum(peak.intensity for peak in self) * self._scaling_factor()


class PeakClusterBuilder(object):
    def __init__(self, error_tolerance=1e-5):
        self.error_tolerance = error_tolerance
        self.peak_clusters = []
        self.count = 0

    def __len__(self):
        return len(self.peak_clusters)

    def __iter__(self):
        return iter(self.peak_clusters)

    def __getitem__(self, i):
        if isinstance(i, (int, slice)):
            return self.peak_clusters[i]
        else:
            return [self.peak_clusters[j] for j in i]

    def has_peak(self, mz, error_tolerance=1e-5):
        ix, matched = binary_search_with_flag(self.peak_clusters, mz, error_tolerance)
        if matched:
            return self[ix]

    def find_insertion_point(self, peak):
        index, matched = binary_search_with_flag(
            self.peak_clusters, peak.mz, self.error_tolerance)
        return index, matched

    def find_minimizing_index(self, peak, indices):
        best_index = None
        best_error = float('inf')
        for index_case in indices:
            chroma = self[index_case]
            err = abs(chroma.mz - peak.mz) / peak.mz
            if err < best_error:
                best_index = index_case
                best_error = err
        return best_index

    def handle_peak(self, peak):
        if len(self) == 0:
            index = [0]
            matched = False
        else:
            index, matched = self.find_insertion_point(peak)
        if matched:
            pc = self.peak_clusters[self.find_minimizing_index(peak, index)]
            pc.add(peak)
        else:
            pc = PeakCluster()
            pc.add(peak)
            self.insert_peak_cluster(pc, index)
        self.count += 1

    def insert_peak_cluster(self, peak_cluster, index):
        if index[0] != 0:
            self.peak_clusters.insert(index[0] + 1, peak_cluster)
        else:
            if len(self) == 0:
                new_index = index[0]
            else:
                x = self.peak_clusters[index[0]]
                if x.mz < peak_cluster.mz:
                    new_index = index[0] + 1
                else:
                    new_index = index[0]
            self.peak_clusters.insert(new_index, peak_cluster)

    def cluster_peaks(self, peaks):
        for peak in sorted(peaks, key=lambda x: x.intensity):
            self.handle_peak(peak)

    def compute_probability(self, n_spectra):
        n_spectra = float(n_spectra)
        for pc in self:
            pc.probability = len(pc) / n_spectra


def cluster_peaks(peak_sets, error_tolerance=1e-5, n_spectra=None):
    if n_spectra is None:
        n_spectra = len(peak_sets)
    clusterer = PeakClusterBuilder(error_tolerance)
    clusterer.cluster_peaks(unspool_peaks(peak_sets))
    clusterer.compute_probability(n_spectra)
    return clusterer


def binary_search_with_flag(array, mz, error_tolerance=1e-5):
        lo = 0
        n = hi = len(array)
        while hi != lo:
            mid = (hi + lo) / 2
            x = array[mid]
            err = (x.mz - mz) / mz
            if abs(err) <= error_tolerance:
                i = mid - 1
                # Begin Sweep forward
                while i > 0:
                    x = array[i]
                    err = (x.mz - mz) / mz
                    if abs(err) <= error_tolerance:
                        i -= 1
                        continue
                    else:
                        break
                low_end = i
                i = mid + 1

                # Begin Sweep backward
                while i < n:
                    x = array[i]
                    err = (x.mz - mz) / mz
                    if abs(err) <= error_tolerance:
                        i += 1
                        continue
                    else:
                        break
                high_end = i
                return list(range(low_end, high_end)), True
            elif (hi - lo) == 1:
                return [mid], False
            elif err > 0:
                hi = mid
            elif err < 0:
                lo = mid
        return 0, False


def get_mz_window(peak_set, lo, hi, error_tolerance=1e-5):
    lo_ix, matched = binary_search_with_flag(peak_set, lo, error_tolerance)
    hi_ix, matched = binary_search_with_flag(peak_set, hi, error_tolerance)

    return peak_set[lo_ix[0]:hi_ix[0] + 1]


def top_n_peaks(peak_set, n=5):
    peak_set = sorted(peak_set, key=lambda x: x.intensity, reverse=True)
    try:
        for i in range(n):
            yield peak_set[i]
    except IndexError:
        pass


def filter_peak_set(peak_set, window_width=100, top_n=5):
    peaks = []
    upper_limit = max(peak_set, key=lambda x: x.mz).mz
    current_mz = 0.1
    while current_mz < upper_limit:
        peaks.extend(
            top_n_peaks(get_mz_window(
                peak_set, current_mz, current_mz + window_width),
                top_n))
        current_mz += window_width
    return sorted(peaks, key=lambda x: x.mz)


def average_fitted_peak(peak_cluster, divisor=False):
    if divisor:
        scale = divisor
    else:
        scale = 1
    weighted_mz = peak_cluster.weighted_mz()
    total_intensity = peak_cluster.intensity
    signal_to_noise = sum(p.signal_to_noise * p.intensity for p in peak_cluster) / total_intensity
    fwhm = sum(p.full_width_at_half_max * p.intensity for p in peak_cluster) / total_intensity
    area = sum(p.area * p.intensity for p in peak_cluster) / total_intensity
    return FittedPeak(weighted_mz, total_intensity / scale, signal_to_noise, -1, -1, fwhm, area)


def average_deconvoluted_peak(peak_cluster, divisor=False):
    if divisor:
        scale = divisor
    else:
        scale = 1
    weighted_mz = peak_cluster.weighted_mz()
    total_intensity = peak_cluster.intensity
    neutral_mass = sum(p.neutral_mass * p.intensity for p in peak_cluster) / total_intensity
    most_abundant_mass = sum(p.most_abundant_mass * p.intensity for p in peak_cluster) / total_intensity
    a_to_a2_ratio = sum(p.a_to_a2_ratio * p.intensity for p in peak_cluster) / total_intensity
    average_mass = sum(p.average_mass * p.intensity for p in peak_cluster) / total_intensity
    signal_to_noise = sum(p.signal_to_noise * p.intensity for p in peak_cluster) / total_intensity
    fwhm = sum(p.full_width_at_half_max * p.intensity for p in peak_cluster) / total_intensity
    area = sum(p.area * p.intensity for p in peak_cluster) / total_intensity
    score = sum(p.score * p.intensity for p in peak_cluster) / total_intensity
    charge_counts = Counter(p.charge for p in peak_cluster).items()
    charge = max(charge_counts, key=lambda x: x[1])[0]
    envelope = average_envelope([p.envelope for p in peak_cluster])
    return DeconvolutedPeak(
        neutral_mass=neutral_mass, intensity=total_intensity / scale, charge=charge, signal_to_noise=signal_to_noise,
        full_width_at_half_max=fwhm, index=-1, a_to_a2_ratio=a_to_a2_ratio, most_abundant_mass=most_abundant_mass,
        average_mass=average_mass, score=score, envelope=envelope, mz=weighted_mz, fit=None, chosen_for_msms=False,
        area=area)


def average_envelope(envelopes, error_tolerance=2e-5):
    for position in cluster_peaks(envelopes, error_tolerance):
        yield (position.mz, position.intensity)


def make_peak_index(fitted_peaks):
    ps = PeakSet(fitted_peaks)
    ps._index()
    return ps


def make_deconvoluted_peak_set(deconvoluted_peaks):
    ps = DeconvolutedPeakSet(deconvoluted_peaks)
    ps.reindex()
    return ps
