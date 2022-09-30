cimport numpy as np


cdef double erf(double x)

cpdef gaussian_shape(np.ndarray[np.float64_t, ndim=1] xs, double center,
                     double amplitude, double sigma)

cpdef bigaussian_shape(np.ndarray[np.float64_t, ndim=1] xs, double center,
                       double amplitude, double sigma_left, double sigma_right)

cpdef skewed_gaussian_shape(np.ndarray[np.float64_t, ndim=1] xs, double center,
                            double amplitude, double sigma, double gamma)
