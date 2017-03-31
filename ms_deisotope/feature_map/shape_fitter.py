from collections import OrderedDict

import numpy as np
from scipy.optimize import leastsq
from numpy import pi, sqrt, exp
from scipy.special import erf

from ms_peak_picker import search

from .profile_transform import smooth_leveled


DEFAULT_SMOOTH = 3


class PeakShapeModelBase(object):
    def __repr__(self):
        return "{self.__class__.__name__}()".format(self=self)


class SkewedGaussianModel(PeakShapeModelBase):
    @staticmethod
    def fit(params, xs, ys):
        center, amplitude, sigma, gamma = params
        return ys - SkewedGaussianModel.shape(xs, center, amplitude, sigma, gamma) * (
            sigma / 2. if abs(sigma) > 2 else 1.)

    @staticmethod
    def guess(xs, ys):
        center = np.average(xs, weights=ys / ys.sum())
        height_at = np.abs(xs - center).argmin()
        apex = ys[height_at]
        sigma = np.abs(center - xs[[search.nearest_left(ys, apex / 2, height_at),
                                    search.nearest_right(ys, apex / 2, height_at + 1)]]).sum()
        gamma = 1
        return center, apex, sigma, gamma

    @staticmethod
    def params_to_dict(params):
        center, amplitude, sigma, gamma = params
        return OrderedDict((("center", center), ("amplitude", amplitude), ("sigma", sigma), ("gamma", gamma)))

    @staticmethod
    def shape(xs, center, amplitude, sigma, gamma):
        norm = (amplitude) / (sigma * sqrt(2 * pi)) * \
            exp(-((xs - center) ** 2) / (2 * sigma ** 2))
        skew = (1 + erf((gamma * (xs - center)) / (sigma * sqrt(2))))
        return norm * skew

    @staticmethod
    def center(params_dict):
        return params_dict['center']

    @staticmethod
    def spread(params_dict):
        return params_dict['sigma']


class PenalizedSkewedGaussianModel(SkewedGaussianModel):
    @staticmethod
    def fit(params, xs, ys):
        center, amplitude, sigma, gamma = params
        return ys - PenalizedSkewedGaussianModel.shape(xs, center, amplitude, sigma, gamma) * (
            sigma / 2. if abs(sigma) > 2 else 1.) * (gamma / 2. if abs(gamma) > 40 else 1.) * (
            center if center > xs[-1] or center < xs[0] else 1.)


class BiGaussianModel(PeakShapeModelBase):

    @staticmethod
    def center(params_dict):
        return params_dict['center']

    @staticmethod
    def spread(params_dict):
        return (params_dict['sigma_left'] + params_dict['sigma_right']) / 2.

    @staticmethod
    def shape(xs, center, amplitude, sigma_left, sigma_right):
        ys = np.zeros_like(xs, dtype=np.float32)
        left_mask = xs < center
        ys[left_mask] = amplitude * np.exp(-(xs[left_mask] - center) ** 2 / (2 * sigma_left ** 2)) * sqrt(2 * pi)
        right_mask = xs > center
        ys[right_mask] = amplitude * np.exp(-(xs[right_mask] - center) ** 2 / (2 * sigma_right ** 2)) * sqrt(2 * pi)
        return ys

    @staticmethod
    def fit(params, xs, ys):
        center, amplitude, sigma_left, sigma_right = params
        return ys - BiGaussianModel.shape(
            xs, center, amplitude, sigma_left, sigma_right) * (
            center if center > xs[-1] or center < xs[0] else 1.)

    @staticmethod
    def params_to_dict(params):
        center, amplitude, sigma_left, sigma_right = params
        return OrderedDict(
            (("center", center), ("amplitude", amplitude), ("sigma_left", sigma_left), ("sigma_right", sigma_right)))

    @staticmethod
    def guess(xs, ys):
        center = np.average(xs, weights=ys / ys.sum())
        height_at = np.abs(xs - center).argmin()
        apex = ys[height_at]
        sigma = np.abs(center - xs[[search.nearest_left(ys, apex / 2, height_at),
                                    search.nearest_right(ys, apex / 2, height_at + 1)]]).sum()
        return center, apex, sigma, sigma


class ChromatogramShapeFitterBase(object):
    def __init__(self, chromatogram, smooth=DEFAULT_SMOOTH, fitter=PenalizedSkewedGaussianModel()):
        self.chromatogram = chromatogram
        self.smooth = smooth
        self.xs = None
        self.ys = None
        self.line_test = None
        self.off_center = None
        self.shape_fitter = fitter

    def handle_invalid(self):
        self.line_test = 0.5

    def extract_arrays(self):
        self.xs, self.ys = self.chromatogram.as_arrays()
        if self.smooth:
            self.ys = smooth_leveled(self.xs, self.ys, self.smooth)
        if len(self.xs) > 2000:
            new_xs = np.linspace(self.xs.min(), self.xs.max(), 2000)
            new_ys = np.interp(new_xs, self.xs, self.ys)
            self.xs = new_xs
            self.ys = new_ys
            self.ys = smooth_leveled(self.xs, self.ys, self.smooth)

    def compute_residuals(self):
        return NotImplemented

    def perform_line_test(self):
        ys = self.ys
        residuals = self.compute_residuals()
        line_test = (residuals ** 2).sum() / (
            ((ys - ((ys.max() + ys.min()) / 2.)) ** 2).sum())
        self.line_test = line_test

    def plot(self, ax=None):
        if ax is None:
            from matplotlib import pyplot as plt
            fig, ax = plt.subplots(1)
        ax.plot(self.xs, self.ys, label='Observed')
        ax.scatter(self.xs, self.ys, label='Observed')
        ax.plot(self.xs, self.compute_fitted(), label='Fitted')
        ax.plot(self.xs, self.compute_residuals(), label='Residuals')
        return ax

    @property
    def fit_parameters(self):
        raise NotImplementedError()


class ChromatogramShapeFitter(ChromatogramShapeFitterBase):
    def __init__(self, chromatogram, smooth=DEFAULT_SMOOTH, fitter=PenalizedSkewedGaussianModel()):
        super(ChromatogramShapeFitter, self).__init__(chromatogram, smooth=smooth, fitter=fitter)

        self.params = None
        self.params_dict = None

        if len(chromatogram) < 5:
            self.handle_invalid()
        else:
            self.extract_arrays()
            self.peak_shape_fit()
            self.perform_line_test()
            self.off_center_factor()

    @property
    def fit_parameters(self):
        return self.params_dict

    def __repr__(self):
        return "ChromatogramShapeFitter(%s, %0.4f)" % (self.chromatogram, self.line_test)

    def off_center_factor(self):
        center = self.shape_fitter.center(self.params_dict)
        spread = self.shape_fitter.spread(self.params_dict)
        self.off_center = abs(1 - abs(1 - (2 * abs(
            self.xs[self.ys.argmax()] - center) / abs(spread))))
        if self.off_center > 1:
            self.off_center = 1. / self.off_center
        self.line_test /= self.off_center

    def compute_residuals(self):
        return self.shape_fitter.fit(self.params, self.xs, self.ys)

    def compute_fitted(self):
        return self.shape_fitter.shape(self.xs, **self.params_dict)

    def peak_shape_fit(self):
        xs, ys = self.xs, self.ys
        params = self.shape_fitter.guess(xs, ys)
        fit = leastsq(self.shape_fitter.fit,
                      params, (xs, ys))
        params = fit[0]
        self.params = params
        self.params_dict = self.shape_fitter.params_to_dict(params)

    def iterfits(self):
        yield self.compute_fitted()


def shape_fit_test(chromatogram, smooth=DEFAULT_SMOOTH):
    return ChromatogramShapeFitter(chromatogram, smooth).line_test


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


class MultimodalChromatogramShapeFitter(ChromatogramShapeFitterBase):
    def __init__(self, chromatogram, max_peaks=5, smooth=DEFAULT_SMOOTH, fitter=BiGaussianModel()):
        super(MultimodalChromatogramShapeFitter, self).__init__(chromatogram, smooth=smooth, fitter=fitter)
        self.max_peaks = max_peaks
        self.params_list = []
        self.params_dict_list = []

        if len(self.chromatogram) < 5:
            self.handle_invalid()
        else:
            self.extract_arrays()
            self.peak_shape_fit()
            self.perform_line_test()

    @property
    def fit_parameters(self):
        return self.params_dict_list

    def __repr__(self):
        return "MultimodalChromatogramShapeFitter(%s, %0.4f)" % (self.chromatogram, self.line_test)

    def peak_shape_fit(self):
        return self.set_up_peak_fit()

    def set_up_peak_fit(self, ys=None, min_height=0, peak_count=0):
        xs = self.xs
        if ys is None:
            ys = self.ys
        params = self.shape_fitter.guess(xs, ys)
        params_dict = self.shape_fitter.params_to_dict(params)

        indices = peak_indices(ys, min_height)
        center = xs[max(indices, key=lambda x: ys[x])]
        params_dict['center'] = center

        fit = leastsq(self.shape_fitter.fit,
                      params_dict.values(), (xs, ys))
        params = fit[0]
        params_dict = self.shape_fitter.params_to_dict(params)
        self.params_list.append(params)
        self.params_dict_list.append(params_dict)

        residuals = self.shape_fitter.fit(params, xs, ys)

        fitted_apex_index = search.get_nearest(xs, params_dict['center'], 0)
        fitted_apex = ys[fitted_apex_index]

        new_min_height = fitted_apex * 0.5

        if new_min_height < min_height:
            min_height *= 0.85
        else:
            min_height = new_min_height

        indices = peak_indices(residuals, min_height)

        peak_count += 1
        if indices.size and peak_count < self.max_peaks:
            residuals, params_dict = self.set_up_peak_fit(residuals, min_height, peak_count=peak_count)

        return residuals, params_dict

    def compute_fitted(self):
        xs = self.xs
        fitted = np.zeros_like(xs)
        for params_dict in self.params_dict_list:
            fitted += self.shape_fitter.shape(xs, **params_dict)
        return fitted

    def compute_residuals(self):
        return self.ys - self.compute_fitted()

    def iterfits(self):
        xs = self.xs
        for params_dict in self.params_dict_list:
            yield self.shape_fitter.shape(xs, **params_dict)


class AdaptiveMultimodalChromatogramShapeFitter(ChromatogramShapeFitterBase):
    def __init__(self, chromatogram, max_peaks=5, smooth=DEFAULT_SMOOTH, fitters=None):
        if fitters is None:
            fitters = (BiGaussianModel(), PenalizedSkewedGaussianModel(),)
        super(AdaptiveMultimodalChromatogramShapeFitter, self).__init__(
            chromatogram, smooth=smooth, fitter=fitters[0])
        self.max_peaks = max_peaks
        self.fitters = fitters
        self.params_list = []
        self.params_dict_list = []

        self.alternative_fits = []
        self.best_fit = None

        if len(self.chromatogram) < 5:
            self.handle_invalid()
        else:
            self.extract_arrays()
            self.peak_shape_fit()
            self.perform_line_test()

    @property
    def fit_parameters(self):
        return self.best_fit.fit_parameters

    def compute_fitted(self):
        return self.best_fit.compute_fitted()

    def compute_residuals(self):
        return self.best_fit.compute_residuals()

    def peak_shape_fit(self):
        for fitter in self.fitters:
            model_fit = ProfileSplittingMultimodalChromatogramShapeFitter(
                self.chromatogram, self.max_peaks, self.smooth, fitter=fitter)
            self.alternative_fits.append(model_fit)
            model_fit = MultimodalChromatogramShapeFitter(
                self.chromatogram, self.max_peaks, self.smooth, fitter=fitter)
            self.alternative_fits.append(model_fit)
        self.best_fit = min(self.alternative_fits, key=lambda x: x.line_test)
        self.params_list = self.best_fit.params_list
        self.params_dict_list = self.best_fit.params_dict_list
        self.shape_fitter = self.best_fit.shape_fitter

    def perform_line_test(self):
        self.line_test = self.best_fit.line_test

    def iterfits(self):
        xs = self.xs
        for params_dict in self.params_dict_list:
            yield self.shape_fitter.shape(xs, **params_dict)

    def __repr__(self):
        return "AdaptiveMultimodalChromatogramShapeFitter(%s, %0.4f)" % (self.chromatogram, self.line_test)


class SplittingPoint(object):
    __slots__ = ["first_maximum", "minimum", "second_maximum", "minimum_index", "total_distance"]

    def __init__(self, first_maximum, minimum, second_maximum, minimum_index):
        self.first_maximum = first_maximum
        self.minimum = minimum
        self.second_maximum = second_maximum
        self.minimum_index = minimum_index
        self.total_distance = self.compute_distance()

    def compute_distance(self):
        return (self.first_maximum - self.minimum) + (self.second_maximum - self.minimum)

    def __repr__(self):
        return "SplittingPoint(%0.4f, %0.4f, %0.4f, %0.2f, %0.3e)" % (
            self.first_maximum, self.minimum, self.second_maximum, self.minimum_index, self.total_distance)


class ProfileSplittingMultimodalChromatogramShapeFitter(ChromatogramShapeFitterBase):
    def __init__(self, chromatogram, max_splits=3, smooth=DEFAULT_SMOOTH, fitter=BiGaussianModel()):
        super(ProfileSplittingMultimodalChromatogramShapeFitter, self).__init__(
            chromatogram, smooth=smooth, fitter=fitter)
        self.max_splits = max_splits
        self.params_list = []
        self.params_dict_list = []
        self.partition_sites = []

        if len(self.chromatogram) < 5:
            self.handle_invalid()
        else:
            self.extract_arrays()
            self.peak_shape_fit()
            self.perform_line_test()

    def __repr__(self):
        return "ProfileSplittingMultimodalChromatogramShapeFitter(%s, %0.4f)" % (self.chromatogram, self.line_test)

    def _extreme_indices(self, ys):
        maxima_indices = peak_indices(ys)
        minima_indices = peak_indices(-ys)
        return maxima_indices, minima_indices

    def locate_extrema(self, xs=None, ys=None):
        if xs is None:
            xs = self.xs
        if ys is None:
            ys = self.ys

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
                    if max_i < min_k < max_j and (y_i - y_k) > (y_i * 0.01) and (
                            y_j - y_k) > (y_j * 0.01):
                        point = SplittingPoint(y_i, y_k, y_j, xs[min_k])
                        candidates.append(point)
        if candidates:
            best_point = max(candidates, key=lambda x: x.total_distance)
            self.partition_sites.append(best_point)

        return candidates

    def build_partitions(self):
        segments = []

        last_x = self.xs.min() - 1
        for point in self.partition_sites:
            mask = (self.xs <= point.minimum_index) & (self.xs > last_x)
            if any(mask):
                segments.append((self.xs[mask], self.ys[mask]))
            last_x = point.minimum_index
        mask = self.xs > last_x
        if any(mask):
            segments.append((self.xs[mask], self.ys[mask]))
        return segments

    def set_up_peak_fit(self, xs, ys):
        params = self.shape_fitter.guess(xs, ys)
        params_dict = self.shape_fitter.params_to_dict(params)
        if len(params) > len(xs):
            self.params_list.append(params)
            self.params_dict_list.append(params_dict)
            return ys, params_dict

        fit = leastsq(self.shape_fitter.fit,
                      params_dict.values(), (xs, ys))
        params = fit[0]
        params_dict = self.shape_fitter.params_to_dict(params)
        self.params_list.append(params)
        self.params_dict_list.append(params_dict)

        residuals = self.shape_fitter.fit(params, xs, ys)
        return residuals, params_dict

    def peak_shape_fit(self):
        self.locate_extrema()
        for segment in self.build_partitions():
            self.set_up_peak_fit(*segment)

    def compute_fitted(self):
        fitted = []
        for segment, params_dict in zip(self.build_partitions(), self.params_dict_list):
            fitted.append(self.shape_fitter.shape(segment[0], **params_dict))
        return np.concatenate(fitted)

    def compute_residuals(self):
        return self.ys - self.compute_fitted()

    def iterfits(self):
        for segment, params_dict in zip(self.build_partitions(), self.params_dict_list):
            yield self.shape_fitter.shape(segment[0], **params_dict)
