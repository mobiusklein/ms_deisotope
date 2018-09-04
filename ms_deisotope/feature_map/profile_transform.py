import numpy as np
from .lcms_feature import LCMSFeature


def peak_indices(x, min_height=0):
    """Find the index of local maxima.

    Parameters
    ----------
    x : np.ndarray
        Data to find local maxima in
    min_height : float, optional
        Minimum peak height

    Returns
    -------
    np.ndarray[int]
        Indices of maxima in x

    References
    ----------
    https://github.com/demotu/BMC/blob/master/functions/detect_peaks.py
    """
    if x.size < 3:
        return np.array([], dtype=int)
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    rising_edges = np.where((np.hstack((dx, 0)) <= 0) &
                            (np.hstack((0, dx)) > 0))[0]
    falling_edges = np.where((np.hstack((dx, 0)) < 0) &
                             (np.hstack((0, dx)) >= 0))[0]
    indices = np.unique(np.hstack((rising_edges, falling_edges)))
    if indices.size and min_height > 0:
        indices = indices[x[indices] >= min_height]
    return indices


class ValleyPoint(object):
    __slots__ = ["first_maximum", "minimum", "second_maximum", "first_maximum_index",
                 "minimum_index", "second_maximum_index", "total_distance"]

    def __init__(self, first_maximum, minimum, second_maximum, first_maximum_index,
                 minimum_index, second_maximum_index):
        self.first_maximum = first_maximum
        self.minimum = minimum
        self.second_maximum = second_maximum
        self.first_maximum_index = first_maximum_index
        self.minimum_index = minimum_index
        self.second_maximum_index = second_maximum_index
        self.total_distance = self.compute_distance()

    def compute_distance(self):
        return (self.first_maximum - self.minimum) + (self.second_maximum - self.minimum)

    def __repr__(self):
        return "ValleyPoint(%0.2f, %0.2f, %0.2f, %0.3e)" % (
            self.first_maximum_index, self.minimum_index, self.second_maximum_index, self.total_distance)

    def split(self, profile):
        index = profile.nodes.find_time(self.minimum_index)[1]
        before = LCMSFeature(profile[:index])
        after = LCMSFeature(profile[index:])
        return before, after


class PeakBoundary(object):
    def __init__(self, first_minimum, maximum, second_minimum,
                 first_minimum_index, maximum_index, second_minimum_index):
        self.first_minimum = first_minimum
        self.maximum = maximum
        self.second_minimum = second_minimum
        self.first_minimum_index = first_minimum_index
        self.maximum_index = maximum_index
        self.second_minimum_index = second_minimum_index
        self.total_distance = self.compute_distance()

    def compute_distance(self):
        return (self.maximum - self.first_minimum) + (self.maximum - self.second_minimum)

    def __repr__(self):
        return "PeakBoundary(%0.2f, %0.2f, %0.2f, %0.3e)" % (
            self.first_minimum_index, self.maximum_index, self.second_minimum_index,
            self.total_distance)

    def split(self, profile):
        first_minimum = profile.nodes.find_time(self.first_minimum_index)[1]
        second_minimum = profile.nodes.find_time(self.second_minimum_index)[1]
        if first_minimum == 0:
            first_minimum += 1
        prefix = profile[:first_minimum]
        peak = profile[(first_minimum - 1):second_minimum + 1]
        suffix = profile[second_minimum + 1:]
        cases = LCMSFeature(prefix), LCMSFeature(peak), LCMSFeature(suffix)
        return cases


def interpolate(xs, ys, n=200):
    new_xs = np.linspace(xs.min(), xs.max(), n)
    new_ys = np.interp(new_xs, xs, ys)
    return new_xs, new_ys


def sliding_mean(ys):
    return (np.concatenate((ys[1:], [0])) + ys + np.concatenate(([0], ys[:-1]))) / 3


def sliding_median(ys):
    arr = np.zeros_like(ys)
    n = len(ys)
    for i in range(n):
        if i == 0:
            value = np.median(np.concatenate(([0], ys[:2]),))
        elif i == n - 1:
            value = np.median(np.concatenate(([0], ys[n - 2:]),))
        else:
            value = np.median(ys[i - 1:i + 2])
        arr[i] = value
    return arr


def gauss(x, mean, variance, scale=2):
    return 1. / np.sqrt(2 * np.pi * variance) * np.exp(-(x - mean) ** 2 / (2 * variance))


def binsearch(array, value):
    lo = 0
    hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        point = array[mid]
        if value == point:
            return mid
        elif hi - lo == 1:
            return mid
        elif point > value:
            hi = mid
        else:
            lo = mid


def gaussian_smooth(x, y, width=0.05):
    smoothed = np.zeros_like(y)

    n = x.shape[0]
    for i in range(n):
        low_edge = binsearch(x, x[i] - width)
        high_edge = binsearch(x, x[i] + width)
        if high_edge - 1 == low_edge or high_edge == low_edge:
            smoothed[i] = y[i]
            continue
        x_slice = x[low_edge:high_edge]
        y_slice = y[low_edge:high_edge]
        center = np.mean(x_slice)
        spread = np.var(x_slice)
        weights = gauss(x_slice, center, spread)
        val = ((y_slice * weights).sum() / weights.sum())
        if np.isnan(val):
            raise ValueError("NaN")
        smoothed[i] = val
    return smoothed


def smooth_leveled(xs, ys, level=0):
    if level == 0:
        return ys
    elif level == 1:
        return sliding_mean(ys)
    elif level == 2:
        return sliding_mean(sliding_median(ys))
    elif level == 3:
        return gaussian_smooth(xs, ys, 0.05)
    else:
        return gaussian_smooth(xs, ys, level)


class ProfileSplitter(object):
    def __init__(self, profile):
        self.profile = profile
        self.xs, self.ys = profile.as_arrays()
        self.partition_sites = []

    def _extreme_indices(self, ys):
        maxima_indices = peak_indices(ys)
        minima_indices = peak_indices(-ys)
        return maxima_indices, minima_indices

    def interpolate(self, xs, ys, n=200):
        new_xs = np.linspace(xs.min(), xs.max(), n)
        new_ys = np.interp(new_xs, xs, ys)
        return new_xs, new_ys

    def locate_peak_boundaries(self, smooth=1, interpolate_past=200):
        xs = self.xs
        ys = self.ys
        if len(xs) > interpolate_past:
            xs, ys = self.interpolate(xs, ys, interpolate_past)
        if smooth:
            ys = smooth_leveled(xs, ys, smooth)

        maxima_indices, minima_indices = self._extreme_indices(ys)
        candidates = []
        for i in range(len(minima_indices)):
            min_i = minima_indices[i]
            for j in range(i + 1, len(minima_indices)):
                min_j = minima_indices[j]
                for k in range(len(maxima_indices)):
                    max_k = maxima_indices[k]
                    y_i = ys[min_i]
                    y_j = ys[min_j]
                    y_k = ys[max_k]
                    if min_i < max_k < min_j and (y_k - y_i) > (y_k * 0.01) and (
                            y_k - y_j) > (y_k * 0.01):
                        point = PeakBoundary(y_i, y_k, y_j, xs[min_i], xs[max_k], xs[min_j])
                        candidates.append(point)
        if candidates:
            candidates = sorted(candidates, key=lambda x: x.total_distance, reverse=True)
            best_point = candidates[0]
            self.partition_sites.append(best_point)

        return candidates

    def locate_valleys(self, scale=0.3, smooth=1, interpolate_past=200):
        xs = self.xs
        ys = self.ys

        if len(xs) > interpolate_past:
            xs, ys = self.interpolate(xs, ys, interpolate_past)

        if smooth:
            ys = smooth_leveled(xs, ys, smooth)

        maxima_indices, minima_indices = self._extreme_indices(ys)
        candidates = []

        for i in range(len(maxima_indices)):
            max_i = maxima_indices[i]
            for j in range(i + 1, len(maxima_indices)):
                max_j = maxima_indices[j]
                for k in range(len(minima_indices)):
                    min_k = minima_indices[k]
                    y_i = ys[max_i]
                    y_j = ys[max_j]
                    y_k = ys[min_k]
                    if max_i < min_k < max_j and (y_i - y_k) > (y_i * scale) and (
                            y_j - y_k) > (y_j * scale):
                        point = ValleyPoint(y_i, y_k, y_j, xs[max_i], xs[min_k], xs[max_j])
                        candidates.append(point)
        if candidates:
            candidates = sorted(candidates, key=lambda x: x.total_distance, reverse=True)
            best_point = candidates[0]
            self.partition_sites.append(best_point)

        return candidates

    def extract_peak(self, peak_boundary):
        cases = peak_boundary.split(self.profile)
        return cases

    def split_valley(self, valley):
        return valley.split(self.profile)


def split_valleys(profiles, scale=0.3, n_levels=2):
    for i in range(n_levels):
        out = []
        for case in profiles:
            ps = ProfileSplitter(case)
            vallies = ps.locate_valleys(scale)
            if vallies:
                splitter = vallies[0]
                out.extend(p for p in ps.split_valley(splitter) if len(p) > 0)
            else:
                out.append(case)
        profiles = out
    return profiles


def exhaustive_split(profiles, scale=0.5, maxiter=100, smooth=.15):
    hold = []
    for i in range(maxiter):
        n = len(profiles) + len(hold)
        out = []
        for case in profiles:
            ps = ProfileSplitter(case)
            vallies = ps.locate_valleys(scale, smooth=smooth)
            if vallies:
                splitter = vallies[0]
                out.extend(p for p in ps.split_valley(splitter) if len(p) > 0)
            else:
                hold.append(case)
        profiles = out
        if n == len(profiles) + len(hold):
            break
    else:
        print("Did not converge in %d iterations (%d != %d)" % (i, n, len(profiles) + len(hold)))
    return profiles + hold


def is_flat(feature, fraction=0.5, smooth=1):
    xs, ys = feature.as_arrays()
    ys = smooth_leveled(xs, ys, smooth)
    apex = ys.max()
    return apex * fraction < ys[0] and apex * fraction < ys[-1]


def clean_profiles(features, scale=0.5, maxiter=100, smooth=0.15):
    features = [c for f in features for c in f.split_sparse(0.25)]
    features = exhaustive_split(features, scale=scale, smooth=smooth)
    features = [fe for fe in features if len(fe) > 2 and not is_flat(
                fe, fraction=0.8, smooth=3)]
    return features


try:
    has_c = True
    _ProfileSplitter = ProfileSplitter
    _ValleyPoint = ValleyPoint
    _split_valleys = split_valleys
    _smooth_leveled = smooth_leveled
    _gaussian_smooth = gaussian_smooth

    from ms_deisotope._c.feature_map.profile_transform import (
        ProfileSplitter, ValleyPoint, split_valleys,
        gaussian_smooth, smooth_leveled)
except ImportError:
    has_c = False
