# cython: embedsignature=True
#

from libc.math cimport sqrt, exp, pi

cimport cython
cimport numpy as np
import numpy as np


from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature
from ms_deisotope._c.feature_map.shape_fitter cimport (
    gaussian_shape,
    bigaussian_shape,
    skewed_gaussian_shape)


np.import_array()



cpdef object find_immediate_minima(np.ndarray[double, ndim=1] x, double scale):
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
        if val > prev * (1 + scale) and val > nxt * (1 + scale):
            maxima.append(i)
        elif val < prev * scale and val < nxt * scale:
            minima.append(i)

    return np.array(maxima, dtype=np.int64), np.array(minima, dtype=np.int64)


cpdef object locate_extrema(np.ndarray[double, ndim=1] x):
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


cpdef sliding_mean(np.ndarray[cython.floating, ndim=1, mode='c'] ys):
    cdef:
        size_t i, n
        np.npy_intp knd
        np.ndarray[double, ndim=1, mode='c'] result
        double vlo, v, vhi

    n = ys.shape[0]
    knd = n
    result = np.PyArray_ZEROS(1, &knd, np.NPY_FLOAT64, 0)
    for i in range(n):
        if i == 0:
            vlo = 0
        else:
            vlo = ys[i - 1]
        v = ys[i]
        if i == (n - 1):
            vhi = 0
        else:
            vhi = ys[i + 1]
        result[i] = (vlo + v + vhi) / 3
    return result



cdef object sliding_median(np.ndarray[double, ndim=1, mode='c'] ys):
    cdef:
        np.ndarray[double, ndim=1, mode='c'] arr
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


cpdef size_t binsearch(np.ndarray[double, ndim=1, mode='c'] array, double value):
    cdef:
        size_t lo, hi, mid
        double point
    lo = 0
    hi = array.shape[0]
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


cdef np.ndarray[double, ndim=1, mode='c'] gauss(np.ndarray[double, ndim=1, mode='c'] x,
                                                double mean, double variance, double scale=2.0):
    return 1. / np.sqrt(2 * np.pi * variance) * np.exp(-(x - mean) ** 2 / (scale * variance))


cpdef np.ndarray[double, ndim=1, mode='c'] gaussian_smooth(
        np.ndarray[double, ndim=1, mode='c'] x,
        np.ndarray[double, ndim=1, mode='c'] y, double width=0.05):
    cdef:
        np.ndarray[double, ndim=1, mode='c'] smoothed, x_slice, y_slice, weights
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


cpdef np.ndarray[double, ndim=1, mode='c'] smooth_leveled(
        np.ndarray[double, ndim=1, mode='c'] xs,
        np.ndarray[double, ndim=1, mode='c'] ys, double level=0):
    if level == 0:
        return ys
    elif level == 1:
        return sliding_mean[double](ys)
    elif level == 2:
        return sliding_mean[double](sliding_median(ys))
    elif level == 3:
        return gaussian_smooth(xs, ys, 0.05)
    else:
        return gaussian_smooth(xs, ys, level)


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
        public object profile
        public list partition_sites
        public np.ndarray xs
        public np.ndarray ys

    def __init__(self, profile):
        if isinstance(profile, LCMSFeature):
            self.profile = profile
            self.xs, self.ys = profile.as_arrays()
        else:
            self.xs, self.ys = profile
            self.profile = None
        self.partition_sites = []

    cpdef _extreme_indices(self, ys):
        maxima_indices, minima_indices = locate_extrema(ys)
        return maxima_indices.astype(np.int64), minima_indices.astype(np.int64)

    cpdef object interpolate(self, xs, ys, int n=200):
        new_xs = np.linspace(xs.min(), xs.max(), n)
        new_ys = np.interp(new_xs, xs, ys)
        return new_xs, new_ys

    @cython.boundscheck(False)
    def locate_peak_boundaries(self, double smooth=1, int interpolate_past=200):
        cdef:
            np.ndarray[double, ndim=1] xs, ys
            np.ndarray[np.int64_t, ndim=1] maxima_indices, minima_indices
            list candidates
            size_t i, j, k, max_i, max_j, min_k
            double y_i, y_j, y_k
            PeakBoundary point

        xs = self.xs
        ys = self.ys
        ys = sliding_mean[double](ys)
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

    @cython.boundscheck(False)
    cpdef locate_valleys(self, double scale=0.3, double smooth=1, int interpolate_past=200):
        cdef:
            np.ndarray[double, ndim=1, mode='c'] xs
            np.ndarray[double, ndim=1, mode='c'] ys
            np.ndarray[np.int64_t, ndim=1, mode='c'] maxima_indices
            np.ndarray[np.int64_t, ndim=1, mode='c'] minima_indices
            list candidates
            size_t i, j, k, max_i, max_j, min_k, n_maxima, n_minima
            double y_i, y_j, y_k
            ValleyPoint point

        xs = self.xs
        ys = self.ys

        if len(xs) > interpolate_past:
            xs, ys = self.interpolate(xs, ys, interpolate_past)

        if smooth:
            ys = smooth_leveled(xs, ys, smooth)

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
                        point = _create_ValleyPoint(y_i, y_k, y_j, xs[max_i], xs[min_k], xs[max_j])
                        candidates.append(point)
        if candidates:
            candidates = sorted(candidates, key=_get_total_distance, reverse=True)
            best_point = candidates[0]
            self.partition_sites.append(best_point)

        return candidates

    def extract_peak(self, peak_boundary):
        if self.profile is not None:
            cases = peak_boundary.split(self.profile)
        else:
            cases = None
        return cases

    def split_valley(self, valley):
        if self.profile is not None:
            return valley.split(self.profile)
        else:
            return None


cpdef double _get_total_distance(object obj):
    return (<ValleyPoint>obj).total_distance


cdef ValleyPoint _create_ValleyPoint(double first_maximum, double minimum, double second_maximum,
                                     double first_maximum_index, double minimum_index,
                                     double second_maximum_index):
    cdef:
        ValleyPoint self
    self = ValleyPoint.__new__(ValleyPoint)
    self.first_maximum = first_maximum
    self.minimum = minimum
    self.second_maximum = second_maximum
    self.first_maximum_index = first_maximum_index
    self.minimum_index = minimum_index
    self.second_maximum_index
    return self


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



@cython.cdivision(True)
@cython.boundscheck(False)
cdef size_t binsearch_view(cython.floating[::1] array, cython.floating x) nogil:
    cdef:
        size_t lo, hi, mid
        cython.floating y, err

    lo = 0
    hi = array.shape[0]

    while hi != lo:
        mid = (hi + lo) / 2
        y = array[mid]
        err = y - x
        if hi - lo == 1:
            return mid
        elif err > 0:
            hi = mid
        else:
            lo = mid
    return 0


ctypedef fused floating1:
    float
    double

ctypedef fused floating2:
    float
    double


@cython.cdivision(True)
@cython.boundscheck(False)
cpdef double interpolate_point(floating1[::1] time, floating2[::1] intensity, double x):
    cdef:
        size_t i, j, n
        double contrib, tmp
        double time_j, time_j1
        double inten_j, inten_j1

    j = binsearch_view(time, x)
    time_j = time[j]
    if time_j <= x and j + 1 < time.shape[0]:
        time_j1 = time[j + 1]
        inten_j = intensity[j]
        inten_j1 = intensity[j + 1]
    elif time_j > x and j > 0:
        time_j1 = time_j
        inten_j1 = intensity[j]
        time_j = time[j - 1]
        inten_j = time[j - 1]
    else:
        return 0.0
    tmp = time_j1 - time_j
    if tmp == 0:
        contrib = 0.0
    else:
        contrib = ((inten_j * (time_j1 - x)) + (inten_j1 * (x - time_j))) / (time_j1 - time_j)
    return contrib
