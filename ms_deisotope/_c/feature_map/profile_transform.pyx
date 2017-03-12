# cython: embedsignature=True
# 

cimport cython
cimport numpy as cnp
import numpy as np


from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature


cdef object locate_extrema(cnp.ndarray[double, ndim=1] x):
    cdef:
        size_t i, n
        list maxima, minima
        double val, prev, nxt

    maxima = []
    minima = []
    n = x.shape[0]
    for i in range(1, n - 1):
        val = x[i]
        prev = x[i - 1]
        nxt = x[i + 1]
        if val > prev and val > nxt:
            maxima.append(i)
        elif val < prev and val < nxt:
            minima.append(i)

    return np.array(maxima, dtype=np.int64), np.array(minima, dtype=np.int64)


def interpolate(xs, ys, n=200):
    new_xs = np.linspace(xs.min(), xs.max(), n)
    new_ys = np.interp(new_xs, xs, ys)
    return new_xs, new_ys


cdef object sliding_mean(cnp.ndarray[double, ndim=1, mode='c'] ys):
    return (np.concatenate((ys[1:], [0])) + ys + np.concatenate(([0], ys[:-1]))) / 3


cdef object sliding_median(cnp.ndarray[double, ndim=1, mode='c'] ys):
    cdef:
        cnp.ndarray[double, ndim=1, mode='c'] arr
        size_t i, n
        double value
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


cdef object smoother(cnp.ndarray[double, ndim=1, mode='c'] ys):
    return sliding_mean((ys))


cdef size_t binsearch(cnp.ndarray[double, ndim=1, mode='c'] array, double value):
    cdef:
        size_t lo, hi, mid
        double point
    lo = 0
    hi = array.shape[0]
    while hi != lo:
        mid = (hi + lo) / 2
        point = array[mid]
        if value == point:
            return mid
        elif hi - lo == 1:
            return mid
        elif point > value:
            hi = mid
        else:
            lo = mid


cdef object gauss(cnp.ndarray[double, ndim=1, mode='c'] x, double mean, double variance, double scale=2.0):
    return 1. / np.sqrt(2 * np.pi * variance) * np.exp(-(x - mean) ** 2 / (scale * variance))


cpdef object gaussian_smooth(cnp.ndarray[double, ndim=1, mode='c'] x,
                             cnp.ndarray[double, ndim=1, mode='c'] y, double width=0.05):
    cdef:
        cnp.ndarray[double, ndim=1, mode='c'] smoothed, x_slice, y_slice, weights
        size_t i, n, low_edge, high_edge
        double center, spread

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
        smoothed[i] = ((y_slice * weights).sum() / weights.sum())
    return smoothed


cdef class PeakBoundary(object):
    cdef:
        public double first_minimum
        public double maximum
        public double second_minimum
        public double first_minimum_index
        public double maximum_index
        public double second_minimum_index
        public double total_distance

    def __init__(self, first_minimum, maximum, second_minimum,
                 first_minimum_index, maximum_index, second_minimum_index):
        self.first_minimum = first_minimum
        self.maximum = maximum
        self.second_minimum = second_minimum
        self.first_minimum_index = first_minimum_index
        self.maximum_index = maximum_index
        self.second_minimum_index = second_minimum_index
        self.total_distance = self.compute_distance()

    cdef double compute_distance(self):
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


cdef class ValleyPoint(object):
    cdef:
        public double first_maximum
        public double minimum
        public double second_maximum
        public double first_maximum_index
        public double minimum_index
        public double second_maximum_index
        public double total_distance

    def __init__(self, first_maximum, minimum, second_maximum, first_maximum_index,
                 minimum_index, second_maximum_index):
        self.first_maximum = first_maximum
        self.minimum = minimum
        self.second_maximum = second_maximum
        self.first_maximum_index = first_maximum_index
        self.minimum_index = minimum_index
        self.second_maximum_index = second_maximum_index
        self.total_distance = self.compute_distance()

    cdef double compute_distance(self):
        return (self.first_maximum - self.minimum) + (self.second_maximum - self.minimum)

    def __repr__(self):
        return "ValleyPoint(%0.2f, %0.2f, %0.2f, %0.3e)" % (
            self.first_maximum_index, self.minimum_index, self.second_maximum_index, self.total_distance)

    def split(self, profile):
        index = profile.nodes.find_time(self.minimum_index)[1]
        before = LCMSFeature(profile[:index])
        after = LCMSFeature(profile[index:])
        return before, after


cdef class ProfileSplitter(object):
    cdef:
        public LCMSFeature profile
        public list partition_sites
        public cnp.ndarray xs
        public cnp.ndarray ys

    def __init__(self, profile):
        self.profile = profile
        self.xs, self.ys = profile.as_arrays()
        self.partition_sites = []

    def _extreme_indices(self, ys):
        maxima_indices, minima_indices = locate_extrema(ys)
        return maxima_indices.astype(np.int64), minima_indices.astype(np.int64)

    @cython.boundscheck(False)
    def locate_peak_boundaries(self):
        cdef:
            cnp.ndarray[double, ndim=1] xs, ys
            cnp.ndarray[cnp.int64_t, ndim=1] maxima_indices, minima_indices
            list candidates
            size_t i, j, k, max_i, max_j, min_k
            double y_i, y_j, y_k
            PeakBoundary point

        xs = self.xs
        ys = self.ys
        ys = sliding_mean(ys)
        if len(xs) > 200:
            xs, ys = interpolate(xs, ys)

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

    @cython.boundscheck(False)
    def locate_valleys(self, double scale=0.3, bint smooth=True):
        cdef:
            cnp.ndarray[double, ndim=1, mode='c'] xs
            cnp.ndarray[double, ndim=1, mode='c'] ys
            cnp.ndarray[cnp.int64_t, ndim=1, mode='c'] maxima_indices
            cnp.ndarray[cnp.int64_t, ndim=1, mode='c'] minima_indices
            list candidates
            size_t i, j, k, max_i, max_j, min_k, n_maxima, n_minima
            double y_i, y_j, y_k
            ValleyPoint point

        xs = self.xs
        ys = self.ys
        if len(xs) > 200:
            xs, ys = interpolate(xs, ys)

        if smooth:
            ys = smoother(ys)

        maxima_indices, minima_indices = self._extreme_indices(ys)
        candidates = []
        n_minima = len(minima_indices)
        n_maxima = len(maxima_indices)

        for i in range(n_maxima):
            max_i = maxima_indices[i]
            for j in range(i + 1, n_maxima):
                max_j = maxima_indices[j]
                for k in range(n_minima):
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
