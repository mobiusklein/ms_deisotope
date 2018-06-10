from libc.math cimport sqrt, exp, pi

cimport cython
import numpy as np
cimport numpy as np

np.import_array()


cdef:
    double sqrt2pi = sqrt(2 * pi)
    double sqrt2 = sqrt(2)


cdef double erf(double x):
    '''https://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/'''
    # constants
    cdef double a1 =  0.254829592
    cdef double a2 = -0.284496736
    cdef double a3 =  1.421413741
    cdef double a4 = -1.453152027
    cdef double a5 =  1.061405429
    cdef double p  =  0.3275911

    # Save the sign of x
    cdef int sign = 1
    if x < 0:
        sign = -1
    x = abs(x)

    # A & S 7.1.26
    cdef double t = 1.0 / (1.0 + p * x)
    cdef double y = 1.0 - ((((( a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x)

    return sign*y


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef skewed_gaussian_shape(np.ndarray[np.float64_t, ndim=1] xs, double center,
                          double amplitude, double sigma, double gamma):
    cdef:
        np.ndarray[np.float64_t, ndim=1] temp
        size_t i, n
        double sigma_squared_2, xi, t, amp_over_sigma_sqrt2pi, sigma_sqrt2

    if sigma == 0:
        sigma = 1e-3

    sigma_squared_2 = (sigma * sqrt(2))
    amp_over_sigma_sqrt2pi = (amplitude) / (sigma * sqrt(2 * pi))
    sigma_sqrt2 = (sigma * sqrt(2))

    temp = (xs - center)

    n = temp.shape[0]
    for i in range(n):
        xi = xs[i]
        t = -(xi ** 2) / (sigma_squared_2)
        temp[i] *=  amp_over_sigma_sqrt2pi * exp(t)

    for i in range(n):
        xi = xs[i]
        t = 1 + erf((gamma * (xi - center)) / sigma_sqrt2)
        temp[i] *= t

    return temp

@cython.boundscheck(False)
@cython.cdivision(True)
cpdef bigaussian_shape(np.ndarray[np.float64_t, ndim=1] xs, double center,
                     double amplitude, double sigma_left, double sigma_right):
    cdef:
        np.ndarray[np.float64_t, ndim=1] temp
        size_t i, n, center_index
        double xi, t, two_sigma_left_squared, two_sigma_right_squared

    if sigma_left == 0:
        sigma_left = 1e-3
    if sigma_right == 0:
        sigma_right = 1e-3

    temp = xs * 0
    n = temp.shape[0]

    two_sigma_left_squared = (2 * sigma_left ** 2)
    two_sigma_right_squared = (2 * sigma_right ** 2)

    for i in range(n):
        xi = xs[i]
        # left side
        if xi < center:
            t = amplitude * exp(-(xi - center) ** 2 / two_sigma_left_squared) * sqrt2pi
        # right side
        else:
            t = amplitude * exp(-(xi - center) ** 2 / two_sigma_right_squared) * sqrt2pi
        temp[i] = t
    return temp


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef gaussian_shape(np.ndarray[np.float64_t, ndim=1] xs, double center,
                   double amplitude, double sigma):
    cdef:
        np.ndarray[np.float64_t, ndim=1] temp
        size_t i, n
        double xi, t, two_sigma_squared

    if sigma == 0:
        sigma = 1e-3

    two_sigma_squared = (2 * sigma ** 2)

    temp = xs * 0
    n = temp.shape[0]

    for i in range(n):
        xi = xs[i]
        t = amplitude * exp(-(xi - center) ** 2 / two_sigma_squared) * sqrt2pi
        temp[i] = t

    return temp
