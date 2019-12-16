'''Peak shape fitting by non-linear least squares for chromatographic peaks.
'''
from collections import OrderedDict

import numpy as np
from numpy import pi, sqrt, exp

from scipy.optimize import leastsq, least_squares
from scipy.special import erf # pylint: disable=no-name-in-module

from ms_peak_picker import search

from .profile_transform import (smooth_leveled, ProfileSplitter)


DEFAULT_SMOOTH = 3


def linear_regression_residuals(x, y):
    r"""Calculate the squared residuals of ordinary least squares
    fit of :math:`y ~ \alpha + \beta x`

    Parameters
    ----------
    x : :class:`np.ndarray`
        The predictor, the time axis
    y : :class:`np.ndarray`
        The response, the intensity axis

    Returns
    -------
    :class:`np.ndarray`
        The squared residuals of the least squares fit of :math:`y` with :math:`\beta x`
    """
    X = np.vstack((np.ones(len(x)), np.array(x))).T
    Y = np.array(y)
    B = np.linalg.inv(X.T.dot(X)).dot(X.T.dot(Y))
    Yhat = X.dot(B)
    return (Y - Yhat) ** 2


class PeakShapeModelBase(object):
    """A simple abstract base class for chromatographic peak shape models
    """
    def __repr__(self):
        return "{self.__class__.__name__}()".format(self=self)

    def __init__(self, bounds=None):
        self._bounds = bounds

    @staticmethod
    def error(params, xs, ys):
        """A callback function for use with :func:`scipy.optimize.leastsq` to compute
        the residuals given a set of parameters, x, and y

        Parameters
        ----------
        params : list
            The model's parameters, a list of floats
        xs : :class:`np.ndarray`
            The time array to predict with respect to.
        ys : :class:`np.ndarray`
            The intensity array to predict against

        Returns
        -------
        :class:`np.ndarray`:
            The residuals of ``ys - self.shape(params, x)`` with optional penalty terms
        """
        raise NotImplementedError()

    @staticmethod
    def guess(xs, ys):
        """Get crude estimates of roughly where to start fitting parameters

        The results are by no means accurate, but will serve as a reasonable starting point
        for :func:`scipy.optimize.leastsq`.

        Parameters
        ----------
        xs : :class:`np.ndarray`
            The time array to predict with respect to.
        ys : :class:`np.ndarray`
            The intensity array to predict against

        Returns
        -------
        list
        """
        raise NotImplementedError()

    @staticmethod
    def shape(xs, *args):
        """Given some input array `X` and some parameters, compute the theoretical output `y`


        Parameters
        ----------
        xs : :class:`np.ndarray`
            The time array to predict over
        *args:
            A list of :class:`float` values which are individual parameters of the peak shape

        Returns
        -------
        :class:`np.ndarray`: ys
            The theoretical peak shape given the position in time and parameters
        """
        raise NotImplementedError()

    @staticmethod
    def center(params_dict):
        """Get the center position parameter for the dominant peak in a
        fitted peak shape

        Parameters
        ----------
        params_dict : dict
            The fitted peak shape parameters

        Returns
        -------
        float
        """
        return params_dict['center']

    @staticmethod
    def spread(params_dict):
        """Get an approximate symmetric variance term for the fitted peak shape.

        Parameters
        ----------
        params_dict : dict
            The fitted peak shape parameters

        Returns
        -------
        float
        """
        return params_dict['sigma']

    @classmethod
    def _default_bounds(cls, params):
        return ([-np.inf] * len(params), [np.inf] * len(params))

    def bounds(self, params):
        if self._bounds is None:
            return self._default_bounds(params)
        return self._bounds

    @classmethod
    def wrap_parameters(cls, params):
        return FittedPeakShape(cls.params_to_dict(params), cls)

    @classmethod
    def adapt_params_to_bounds(cls, params, bounds):
        new_params = params[:]
        for i, param in enumerate(params):
            if param < bounds[0][i]:
                new_params[i] = bounds[0][i]
            if param > bounds[1][i]:
                new_params[i] = bounds[1][i]
        return new_params

    def leastsq(self, xs, ys, params=None, method='leastsq'):
        if params is None:
            params = self.guess(xs, ys)
        if method == 'least_squares':
            bounds = self.bounds(params)
            params = self.adapt_params_to_bounds(params, bounds)
            result = least_squares(self.error, params, bounds=bounds, args=(xs, ys))
            return result['x'],
        else:
            result = leastsq(self.error, params, args=(xs, ys))
            return result[0],


class SkewedGaussianModel(PeakShapeModelBase):
    r"""Model for a Skewed Gaussian Peak Shape

    .. math::
        \frac{A}{\sigma\sqrt{2\pi}}\exp{-\frac{(x-\mu)^2}{2\sigma^2} +
        1 + \mathbf{erf}{\frac{\gamma * (x - \mu)}{\sigma\sqrt{2}}}`

    """
    @staticmethod
    def error(params, xs, ys):
        center, amplitude, sigma, gamma = params
        return ys - SkewedGaussianModel.shape(xs, center, amplitude, sigma, gamma) * (
            sigma / 2. if abs(sigma) > 2 else 1.)

    @classmethod
    def _default_bounds(cls, params):
        return ([-np.inf, 0, 0, -np.inf], [np.inf] * 4)

    @staticmethod
    def guess(xs, ys):
        weights = np.clip(ys, 0, np.infty)
        center = np.average(xs, weights=weights / weights.sum())
        if np.isnan(center):
            center = xs.mean()
        height_at = np.abs(xs - center).argmin()
        apex = ys[height_at]
        max_ix = ys.argmax()
        if apex < ys[max_ix] * 0.1:
            apex = ys[max_ix]
            height_at = max_ix
        sigma = np.abs(center - xs[[search.nearest_left(ys, apex / 2, height_at),
                                    search.nearest_right(ys, apex / 2, height_at + 1)]]).sum()
        gamma = 1
        return center, apex, sigma, gamma

    @staticmethod
    def params_to_dict(params):
        center, amplitude, sigma, gamma = params
        return OrderedDict((("center", center), ("amplitude", amplitude), ("sigma", sigma), ("gamma", gamma)))

    @staticmethod
    def shape(xs, center, amplitude, sigma, gamma):  # pylint: disable=arguments-differ
        norm = (amplitude) / (sigma * sqrt(2 * pi)) * \
            exp(-((xs - center) ** 2) / (2 * sigma ** 2))
        skew = (1 + erf((gamma * (xs - center)) / (sigma * sqrt(2))))
        return norm * skew


class PenalizedSkewedGaussianModel(SkewedGaussianModel):
    """A penalized version of :class:`SkewedGaussianModel` which applies an extra
    penalty on the fit when any of the shape or position parameters take on extreme
    values.
    """
    @staticmethod
    def error(params, xs, ys):
        center, amplitude, sigma, gamma = params
        return ys - PenalizedSkewedGaussianModel.shape(xs, center, amplitude, sigma, gamma) * (
            sigma / 2. if abs(sigma) > 2 else 1.) * (gamma / 2. if abs(gamma) > 40 else 1.) * (
                center if center > xs[-1] or center < xs[0] else 1.)


class BiGaussianModel(PeakShapeModelBase):
    """Fit an asymmetric Gaussian peak shape model with different variance on either
    side of the mean.
    """
    @staticmethod
    def center(params_dict):
        return params_dict['center']

    @staticmethod
    def spread(params_dict):
        return (params_dict['sigma_left'] + params_dict['sigma_right']) / 2.

    @staticmethod
    def shape(xs, center, amplitude, sigma_left, sigma_right):  # pylint: disable=arguments-differ
        ys = np.zeros_like(xs, dtype=np.float32)
        left_mask = xs < center
        ys[left_mask] = amplitude / \
            sqrt(2 * pi) * np.exp(-(xs[left_mask] -
                                    center) ** 2 / (2 * sigma_left ** 2))
        right_mask = xs >= center
        ys[right_mask] = amplitude / \
            sqrt(2 * pi) * np.exp(-(xs[right_mask] -
                                    center) ** 2 / (2 * sigma_right ** 2))
        return ys

    @staticmethod
    def error(params, xs, ys):
        center, amplitude, sigma_left, sigma_right = params
        return ys - BiGaussianModel.shape(
            xs, center, amplitude, sigma_left, sigma_right) * (
                center if center > xs[-1] or center < xs[0] else 1.)

    @classmethod
    def _default_bounds(cls, params):
        return ([-np.inf, 0, 0, 0], [np.inf, np.inf, 2, 2])

    @staticmethod
    def params_to_dict(params):
        center, amplitude, sigma_left, sigma_right = params
        return OrderedDict(
            (("center", center), ("amplitude", amplitude), ("sigma_left", sigma_left), ("sigma_right", sigma_right)))

    @staticmethod
    def guess(xs, ys):
        weights = np.clip(ys, 0, np.infty)
        center = np.average(xs, weights=weights / weights.sum())
        if np.isnan(center):
            center = xs.mean()
        height_at = np.abs(xs - center).argmin()
        apex = ys[height_at]
        max_ix = ys.argmax()
        if apex < ys[max_ix] * 0.1:
            apex = ys[max_ix]
            height_at = max_ix
        sigma = np.abs(center - xs[[search.nearest_left(ys, apex / 2, height_at),
                                    search.nearest_right(ys, apex / 2, height_at + 1)]]).sum()
        return center, apex, sigma, sigma


class GaussianModel(PeakShapeModelBase):
    """Fit an asymmetric Gaussian peak shape model with different variance on either
    side of the mean.
    """
    @staticmethod
    def center(params_dict):
        return params_dict['center']

    @staticmethod
    def spread(params_dict):
        return (params_dict['sigma'])

    @staticmethod
    def shape(xs, center, amplitude, sigma):  # pylint: disable=arguments-differ
        ys = amplitude / sqrt(2 * pi) * np.exp(-(xs - center) ** 2 / (2 * sigma ** 2))
        return ys

    @staticmethod
    def error(params, xs, ys):
        center, amplitude, sigma = params
        return ys - GaussianModel.shape(
            xs, center, amplitude, sigma) * (
                center if center > xs[-1] or center < xs[0] else 1.)

    @staticmethod
    def params_to_dict(params):
        center, amplitude, sigma = params
        return OrderedDict(
            (("center", center), ("amplitude", amplitude), ("sigma", sigma), ))

    @staticmethod
    def guess(xs, ys):
        weights = np.clip(ys, 0, np.infty)
        center = np.average(xs, weights=weights / weights.sum())
        if np.isnan(center):
            center = xs.mean()
        height_at = np.abs(xs - center).argmin()
        apex = ys[height_at]
        max_ix = ys.argmax()
        if apex < ys[max_ix] * 0.1:
            apex = ys[max_ix]
            height_at = max_ix
        sigma = np.abs(center - xs[[search.nearest_left(ys, apex / 2, height_at),
                                    search.nearest_right(ys, apex / 2, height_at + 1)]]).sum()
        return center, apex, sigma

    @classmethod
    def _default_bounds(cls, params):
        return ([-np.inf, 0, 0], [np.inf] * 3)


class SimpleGaussianModel(GaussianModel):

    @staticmethod
    def error(params, xs, ys):
        center, amplitude = params
        sigma = 1.0
        return ys - GaussianModel.shape(
            xs, center, amplitude, sigma) * (
                center if center > xs[-1] or center < xs[0] else 1.)

    @staticmethod
    def params_to_dict(params):
        center, amplitude = params
        return OrderedDict(
            (("center", center), ("amplitude", amplitude), ))

    @staticmethod
    def guess(xs, ys):
        weights = np.clip(ys, 0, np.infty)
        center = np.average(xs, weights=weights / weights.sum())
        if np.isnan(center):
            center = xs.mean()
        height_at = np.abs(xs - center).argmin()
        apex = ys[height_at]
        max_ix = ys.argmax()
        if apex < ys[max_ix] * 0.1:
            apex = ys[max_ix]
            height_at = max_ix
        return center, apex

    @classmethod
    def _default_bounds(cls, params):
        return ([-np.inf, 0, ], [np.inf] * 2)


try:
    from ms_deisotope._c.feature_map.shape_fitter import gaussian_shape, bigaussian_shape, skewed_gaussian_shape
    SkewedGaussianModel.shape = staticmethod(skewed_gaussian_shape)
    BiGaussianModel.shape = staticmethod(bigaussian_shape)
    GaussianModel.shape = staticmethod(gaussian_shape)
except ImportError as err:
    print(err)


class FittedPeakShape(object):
    def __init__(self, params, shape_model):
        self.params = params
        self.shape_model = shape_model

    def keys(self):
        return self.params.keys()

    def values(self):
        return self.params.values()

    def items(self):
        return self.params.items()

    def __iter__(self):
        return iter(self.params)

    def shape(self, xs):
        return self.shape_model.shape(xs, **self.params)

    def __getitem__(self, key):
        return self.params[key]

    def __setitem__(self, key, value):
        self.params[key] = value

    def __repr__(self):
        return "FittedPeakShape({params}, {self.shape_model})".format(
            self=self, params=", ".join("%s=%0.3f" % (k, v) for k, v in self.params.items()))

    @property
    def center(self):
        return self['center']

    @property
    def amplitude(self):
        return self['amplitude']

    def copy(self):
        """Create a deep copy of this object.

        Returns
        -------
        :class:`FittedPeakShapeModel`
        """
        return self.__class__(self.params.copy(), self.shape_model)


class ChromatogramShapeFitterBase(object):
    """Abstract base class for chromatographic peak shape fitting models.

    Attributes
    ----------
    chromatogram: :class:`~.ms_deisotope.feature_map.lcms_feature.LCMSFeature`-like
        The feature whose shape will be modeled. Expected to have an `as_arrays` method
        returning a time and intensity array.
    xs: :class:`np.ndarray`
        The time array along
    ys: :class:`np.ndarray`
        The intensity array to fit against. May be smoothed.
    line_test: float
        The final test score
    off_center: float
        A term describing the distance between the center of a modeled peak shape and the
        empirical maximum scaled by the width of the peak.
    shape_fitter: :class:`PeakShapeModelBase` class
        The peak shape model type to use.
    """
    def __init__(self, chromatogram, smooth=DEFAULT_SMOOTH, fitter=PenalizedSkewedGaussianModel()):
        self.chromatogram = chromatogram
        self.smooth = smooth
        self.xs = None
        self.ys = None
        self.line_test = None
        self.off_center = None
        self.shape_fitter = fitter

    def handle_invalid(self):
        """Set the score for the model when it would be invalid to use try to fit
        a peak shape.

        This defaults to setting :attr:`line_test` to 0.5

        """
        self.line_test = 0.5

    def extract_arrays(self):
        """Extract :attr:`xs` and :attr:`ys` from :attr:`chromatogram` and optionally
        smooth and downsample them.

        Downsampling is done if there are more than 2000 points to consider.
        """
        self.xs, self.ys = self.chromatogram.as_arrays()
        if self.smooth:
            self.ys = smooth_leveled(self.xs, self.ys, self.smooth)
        if len(self.xs) > 2000:
            new_xs = np.linspace(self.xs.min(), self.xs.max(), 2000)
            new_ys = np.interp(new_xs, self.xs, self.ys)
            self.xs = new_xs
            self.ys = new_ys
            self.ys = smooth_leveled(self.xs, self.ys, self.smooth)

    def compute_fitted(self):
        """Compute the aggregate peak shape model theoretical signal.

        Returns
        -------
        :class:`np.ndarray`
        """
        raise NotImplementedError()

    def compute_residuals(self):
        """Compute the difference between :attr:`ys` and the aggregate peak shape
        model theoretical signal.

        Returns
        -------
        :class:`np.ndarray`
        """
        raise NotImplementedError()

    def perform_line_test(self):
        """Compute a heuristic score describing how well the peak shape fits the
        signal versus a straight line fit.

        Sets :attr:`line_test` to the resulting value.

        """
        residuals = self.compute_residuals()
        null_residuals = linear_regression_residuals(self.xs, self.ys)
        line_test = (residuals ** 2).sum() / (null_residuals).sum()
        if line_test > 1.0:
            line_test = 1.0
        self.line_test = line_test

    def plot(self, ax=None):
        """Draws the peak shape fit.
        """
        if ax is None:
            from matplotlib import pyplot as plt
            _fig, ax = plt.subplots(1)
        ax.plot(self.xs, self.ys, label='Observed')
        ax.scatter(self.xs, self.ys, label='Observed')
        ax.plot(self.xs, self.compute_fitted(), label='Fitted')
        ax.plot(self.xs, self.compute_residuals(), label='Residuals')
        return ax

    @property
    def fit_parameters(self):
        """The fitted peak shapes with their parameters

        Returns
        -------
        list
        """
        raise NotImplementedError()


class ChromatogramShapeFitter(ChromatogramShapeFitterBase):
    """Fit a single peak in signal-over-time data with a peak shape model,
    with optional smoothing and provide a goodness-of-fit score.

    Used as a simplest baseline comparison.

    """
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
        return [self.params_dict]

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
        return self.shape_fitter.error(self.params, self.xs, self.ys)

    def compute_fitted(self):
        return self.shape_fitter.shape(self.xs, **self.params_dict)

    def peak_shape_fit(self):
        xs, ys = self.xs, self.ys
        params = self.shape_fitter.guess(xs, ys)
        fit = self.shape_fitter.leastsq(xs, ys, params)
        params = fit[0]
        self.params = params
        self.params_dict = self.shape_fitter.wrap_parameters(params)

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
    def __init__(self, chromatogram, max_peaks=5, smooth=DEFAULT_SMOOTH, fitter=BiGaussianModel(),
                 relative_peak_height_threshold=0.5, initial_fits=None):
        if initial_fits is None:
            initial_fits = []
        super(MultimodalChromatogramShapeFitter, self).__init__(chromatogram, smooth=smooth, fitter=fitter)
        self.max_peaks = max_peaks
        self.relative_peak_height_threshold = relative_peak_height_threshold
        self.params_dict_list = list(initial_fits)

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
        # If the guessed amplitude is less than 0, we are virtually guaranteed to never
        # improve.
        if params_dict['amplitude'] < 1:
            params_dict['amplitude'] = 1
        indices = peak_indices(ys, min_height)
        if len(indices) > 0:
            center = xs[max(indices, key=lambda x: ys[x])]
        else:
            center = xs[len(xs) // 2]
        params_dict['center'] = center

        fit = self.shape_fitter.leastsq(xs, ys, list(params_dict.values()))
        params = fit[0]
        params_dict = self.shape_fitter.wrap_parameters(params)
        self.params_dict_list.append(params_dict)

        residuals = self.shape_fitter.error(params, xs, ys)

        fitted_apex_index = search.get_nearest(xs, params_dict['center'], 0)
        fitted_apex = ys[fitted_apex_index]

        new_min_height = fitted_apex * self.relative_peak_height_threshold

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
    def __init__(self, chromatogram, max_peaks=5, smooth=DEFAULT_SMOOTH, fitters=None, relative_peak_height_threshold=0.5):
        if fitters is None:
            fitters = (BiGaussianModel(), PenalizedSkewedGaussianModel(),)
        super(AdaptiveMultimodalChromatogramShapeFitter, self).__init__(
            chromatogram, smooth=smooth, fitter=fitters[0])
        self.max_peaks = max_peaks
        self.relative_peak_height_threshold = relative_peak_height_threshold
        self.fitters = fitters
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
        return self.params_dict_list

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
                self.chromatogram, self.max_peaks, self.smooth, fitter=fitter,
                relative_peak_height_threshold=self.relative_peak_height_threshold)
            self.alternative_fits.append(model_fit)
        self.best_fit = min(self.alternative_fits, key=lambda x: x.line_test)
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


class ProfileSplittingMultimodalChromatogramShapeFitter(ChromatogramShapeFitterBase):
    def __init__(self, chromatogram, max_splits=3, smooth=DEFAULT_SMOOTH, fitter=BiGaussianModel(), initial_fits=None):
        super(ProfileSplittingMultimodalChromatogramShapeFitter, self).__init__(
            chromatogram, smooth=smooth, fitter=fitter)
        if initial_fits is None:
            initial_fits = []
        self.max_splits = max_splits
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

        splitter = ProfileSplitter((xs, ys))
        candidates = splitter.locate_valleys(scale=0.01, smooth=self.smooth, interpolate_past=2000)
        if candidates:
            best_point = candidates[0]
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
        params_dict = self.shape_fitter.wrap_parameters(params)
        if len(params) > len(xs):
            self.params_dict_list.append(params_dict)
            return ys, params_dict

        fit = self.shape_fitter.leastsq(xs, ys, list(params_dict.values()))
        params = fit[0]
        params_dict = self.shape_fitter.wrap_parameters(params)
        self.params_dict_list.append(params_dict)

        residuals = self.shape_fitter.error(params, xs, ys)
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


class CentroidFit(object):
    def __init__(self, center, weight, fits=None):
        if fits is None:
            fits = []
        self.center = center
        self.weight = weight
        self.fits = fits

    def __repr__(self):
        return "CentroidFit(%f, %0.3e, %d)" % (self.center, self.weight, len(self.fits))

    def add(self, fit_params):
        self.fits.append(fit_params)
        self.center, self.weight = self.reweight_center()

    def reweight_center(self):
        center_acc = 0
        weight_acc = 0
        for fit in self.fits:
            center_acc = fit['center'] * fit['amplitude']
            weight_acc = fit['amplitude']
        return center_acc / weight_acc, weight_acc

    @classmethod
    def fromparams(cls, params):
        center = params['center']
        weight = params['amplitude']
        return cls(center, weight, [params])


class ProfileSet(object):
    def __init__(self, features):
        self.features = list(features)
        self.fits = list(map(AdaptiveMultimodalChromatogramShapeFitter, self.features))
        self.baselines = self.compute_baselines()
        self.binned_fits = self.overlap_apexes()

    def compute_baselines(self):
        baselines = []
        for fit in self.fits:
            baselines.append(fit.ys[fit.ys < fit.ys.mean()].mean())
        return baselines

    def overlap_apexes(self, error_tolerance=2e-3):
        centroid_fits = []
        for fit in self.fits:
            centroid_fits.extend(fit.params_dict_list)

        centroid_fits.sort(key=lambda x: x['center'])

        binned_fits = []
        last_fit = CentroidFit.fromparams(centroid_fits[0])
        for fit in centroid_fits[1:]:
            if abs((last_fit.center - fit['center']) / fit['center']) < error_tolerance:
                last_fit.add(fit)
            else:
                binned_fits.append(last_fit)
                last_fit = CentroidFit.fromparams(fit)
        binned_fits.append(last_fit)
        binned_fits.sort(key=lambda x: x.weight, reverse=True)
        return binned_fits

    @staticmethod
    def find_right_intersect(vec, target_val, start_index=0):
        nearest_index = start_index
        next_index = start_index

        size = len(vec) - 1
        if next_index == size:
            return size

        next_val = vec[next_index]
        best_distance = np.abs(next_val - target_val)
        while (next_index < size):
            next_index += 1
            next_val = vec[next_index]
            dist = np.fabs(next_val - target_val)  # pylint: disable=assignment-from-no-return
            if dist < best_distance:
                best_distance = dist
                nearest_index = next_index
            if next_index == size or next_val < target_val:
                break
        return nearest_index

    @staticmethod
    def find_left_intersect(vec, target_val, start_index=0):
        nearest_index = start_index
        next_index = start_index

        size = len(vec) - 1
        if next_index == size:
            return size

        next_val = vec[next_index]
        best_distance = np.abs(next_val - target_val)
        while (next_index > 0):
            next_index -= 1
            next_val = vec[next_index]
            dist = np.fabs(  # pylint: disable=assignment-from-no-return
                next_val - target_val)
            if dist < best_distance:
                best_distance = dist
                nearest_index = next_index
            if next_index == size or next_val < target_val:
                break
        return nearest_index

    def find_intersects(self, fit_bin=0):
        starts = []
        ends = []
        bin_fit = self.binned_fits[fit_bin]
        for i, feat_fit in enumerate(self.fits):
            baseline = self.baselines[i]
            xs = feat_fit.xs
            ys = feat_fit.compute_fitted()
            center_ix = search.get_nearest(xs, bin_fit.center, 0)
            left_ix = self.find_left_intersect(ys, baseline, center_ix)
            right_ix = self.find_right_intersect(ys, baseline, center_ix)
            starts.append(xs[left_ix])
            ends.append(xs[right_ix])
        return starts, ends

    def find_bounds(self, fit_bin=0):
        starts, ends = self.find_intersects(fit_bin)
        start_acc = 0
        end_acc = 0
        weight = 0
        for i, f in enumerate(self.features):
            start_acc += starts[i] * f.intensity
            end_acc += ends[i] * f.intensity
            weight += f.intensity
        return start_acc / weight, end_acc / weight

    def split(self, fit_bin=0):
        start, end = self.find_bounds(fit_bin)
        before = []
        apex = []
        after = []
        for feat in self.features:
            before_part, rest = feat.split_at(start)
            if len(before_part) > 0:
                before.append(before_part)
            apex_part, after_part = rest.split_at(end)
            if len(apex_part) > 0:
                apex.append(apex_part)
            if len(after_part) > 0:
                after.append(after_part)
        return before, apex, after
