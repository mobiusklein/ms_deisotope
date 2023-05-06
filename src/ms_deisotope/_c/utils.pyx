# cython: embedsignature=True

cimport cython
from cpython cimport Py_INCREF
from cpython.list cimport PyList_Append, PyList_GET_ITEM, PyList_GET_SIZE, PyList_New, PyList_SET_ITEM

from libc.math cimport sqrt


import numpy as np
cimport numpy as np

from numpy.math cimport isnan, NAN

from ms_peak_picker._c.peak_index cimport PeakIndex
from ms_peak_picker._c.peak_set cimport PeakSet, FittedPeak, PeakSetIndexed

from ms_deisotope._c.peak_set cimport (
    Envelope, EnvelopePair, DeconvolutedPeak, DeconvolutedPeakSet, DeconvolutedPeakSetIndexed, PeakBase)

from ms_deisotope._c.averagine cimport neutral_mass


np.import_array()


@cython.cdivision(True)
cdef double ppm_error(double x, double y):
    return (x - y) / y


cdef double ppm2da(double mass, double error_tolerance):
    x = mass + mass * error_tolerance
    return x - mass


@cython.cdivision(True)
cdef double da2ppm(double mass, double error_tolerance):
    x = mass + error_tolerance
    return (x - mass) / mass


cpdef intensity_getter(peak):
    return (<PeakBase>peak).intensity


cpdef mz_getter(peak):
    return (<PeakBase>peak).mz


cpdef mass_getter(peak):
    return (<DeconvolutedPeak>peak).neutral_mass



@cython.boundscheck(False)
cpdef list decode_envelopes(np.ndarray[np.float32_t, ndim=1] array):
    cdef:
        list envelope_list, current_envelope
        tuple current_envelope_tuple
        size_t i, n, n_members
        np.float32_t a, b
    envelope_list = []
    current_envelope = []
    i = 0
    n = array.shape[0]
    n_members = 0
    while i < n:
        a = array[i]
        b = array[i + 1]
        i += 2
        if a == 0 and b == 0:
            if current_envelope is not None:
                if n_members > 0:
                    current_envelope_tuple = tuple(current_envelope)
                    PyList_Append(envelope_list, Envelope._create(current_envelope_tuple))
                elif i > 2:
                    PyList_Append(envelope_list, Envelope._create(
                            ()
                        )
                    )
                current_envelope = []
                n_members = 0
        else:
            PyList_Append(current_envelope, EnvelopePair._create(a, b))
            n_members += 1
    current_envelope_tuple = tuple(current_envelope)
    envelope_list.append(Envelope._create(current_envelope_tuple))
    return envelope_list


cpdef np.ndarray[np.float32_t] envelopes_to_array(list envelope_list):
    cdef:
        size_t n, m, i, j, point_count, k
        np.ndarray[np.float32_t] collection
        Envelope envelope
        EnvelopePair pair
    n = PyList_GET_SIZE(envelope_list)
    point_count = 0
    for i in range(n):
        envelope = <Envelope>PyList_GET_ITEM(envelope_list, i)
        point_count += (1 + envelope.get_size()) * 2
    collection = np.zeros(point_count, dtype=np.float32)
    k = 0
    for i in range(n):
        k += 2
        envelope = <Envelope>PyList_GET_ITEM(envelope_list, i)
        m = envelope.get_size()
        for j in  range(m):
            pair = envelope.getitem(j)
            collection[k] = pair.mz
            k += 1
            collection[k] = pair.intensity
            k += 1
    return collection


@cython.boundscheck(False)
cpdef DeconvolutedPeakSetIndexed deserialize_deconvoluted_peak_set(dict scan_dict, bint include_envelopes=True):
    cdef:
        list envelopes, peaks
        np.ndarray[np.float64_t] mz_array
        np.ndarray[np.float64_t] intensity_array
        np.ndarray[np.int8_t] charge_array
        np.ndarray[np.float64_t] score_array
        size_t n, i
        double mz, peak_neutral_mass
        int charge
        DeconvolutedPeak peak
        DeconvolutedPeakSetIndexed peak_set

    peaks = []
    if include_envelopes:
        envelopes = decode_envelopes(scan_dict["isotopic envelopes array"])
    else:
        envelopes = None
    mz_array = scan_dict['m/z array'].astype(np.float64, copy=False)
    intensity_array = scan_dict['intensity array'].astype(np.float64, copy=False)
    charge_array = scan_dict['charge array'].astype(np.int8, copy=False)
    score_array = scan_dict['deconvolution score array'].astype(np.float64, copy=False)
    n = mz_array.shape[0]
    i = 0
    for i in range(n):
        mz = mz_array[i]
        charge = charge_array[i]
        peak_neutral_mass = neutral_mass(mz, charge)

        peak = DeconvolutedPeak._create_simple(
            peak_neutral_mass,
            intensity_array[i],
            charge_array[i],
            score_array[i],
            mz,
            <Envelope>PyList_GET_ITEM(envelopes, i) if include_envelopes else None
        )
        peaks.append(peak)
    peak_set = DeconvolutedPeakSetIndexed(peaks)
    peak_set.reindex()
    return peak_set

cpdef PeakIndex deserialize_peak_set(dict scan_dict):
    cdef:
        np.ndarray[np.float64_t] mz_array
        np.ndarray[np.float64_t] intensity_array
        size_t n, i
        FittedPeak peak
        list peaks
    mz_array = scan_dict['m/z array'].astype(np.float64)
    intensity_array = scan_dict['intensity array'].astype(np.float64)
    n = len(mz_array)
    peaks = []
    for i in range(n):
        peak = FittedPeak._create(
            mz_array[i], intensity_array[i], intensity_array[i],
            0.01, 0.01, 0.1, 0, 0,
            intensity_array[i])
        peaks.append(peak)
    peak_set = PeakSetIndexed(peaks)
    peak_set.reindex()
    return PeakIndex(np.array([]), np.array([]), peak_set)


@cython.boundscheck(False)
cpdef DeconvolutedPeakSetIndexed build_deconvoluted_peak_set_from_arrays(np.ndarray[double, ndim=1] mz_array,
                                                                         np.ndarray[double, ndim=1] intensity_array,
                                                                         np.ndarray[long, ndim=1] charge_array):
    cdef:
        list peaks
        size_t n, i
        double mz, peak_neutral_mass
        int charge
        DeconvolutedPeak peak
        DeconvolutedPeakSetIndexed peak_set

    peaks = []
    n = mz_array.shape[0]
    i = 0
    for i in range(n):
        mz = mz_array[i]
        charge = charge_array[i]
        peak_neutral_mass = neutral_mass(mz, charge)

        peak = DeconvolutedPeak._create_simple(
            peak_neutral_mass,
            intensity_array[i],
            charge_array[i],
            intensity_array[i],
            mz,
            None)
        peaks.append(peak)
    peak_set = DeconvolutedPeakSetIndexed(peaks)
    peak_set.reindex()
    return peak_set


@cython.binding(True)
cpdef double _peak_sequence_tic(self, peak_collection peaks) except -1:
    cdef:
        size_t i, n
        PeakBase peak
        double tic
        list py_peaks

    tic = 0.0
    if peak_collection is PeakSet or peak_collection is PeakIndex or peak_collection is DeconvolutedPeakSet or peak_collection is DeconvolutedPeakSetIndexed:
        n = peaks.get_size()
    else:
        py_peaks = list(peaks)
        n = PyList_GET_SIZE(py_peaks)
        if n > 0:
            if not isinstance(<object>PyList_GET_ITEM(py_peaks, 0), PeakBase):
                raise TypeError("Cannot interpret %r as a PeakBase object" % (type(py_peaks[0])))
    for i in range(n):
        if peak_collection is PeakSet or peak_collection is DeconvolutedPeakSet or peak_collection is DeconvolutedPeakSetIndexed:
            peak = peaks.getitem(i)
        elif peak_collection is PeakIndex:
            peak = peaks.peaks.getitem(i)
        else:
            peak = <PeakBase>PyList_GET_ITEM(py_peaks, i)
        tic += peak.intensity
    return tic


@cython.binding(True)
cpdef PeakBase _peak_sequence_bp(self, peak_collection peaks):
    cdef:
        size_t i, n
        PeakBase peak, base_peak
        double max_intensity
        list py_peaks

    max_intensity = 0.0
    base_peak = None
    if peak_collection is PeakSet or peak_collection is PeakIndex or peak_collection is DeconvolutedPeakSet or peak_collection is DeconvolutedPeakSetIndexed:
        n = peaks.get_size()
    else:
        py_peaks = list(peaks)
        n = PyList_GET_SIZE(py_peaks)
        if n > 0:
            if not isinstance(<object>PyList_GET_ITEM(py_peaks, 0), PeakBase):
                raise TypeError("Cannot interpret %r as a PeakBase object" % (type(py_peaks[0])))
    for i in range(n):
        if peak_collection is PeakSet or peak_collection is DeconvolutedPeakSet or peak_collection is DeconvolutedPeakSetIndexed:
            peak = peaks.getitem(i)
        elif peak_collection is PeakIndex:
            peak = peaks.peaks.getitem(i)
        else:
            peak = <PeakBase>PyList_GET_ITEM(py_peaks, i)
        if max_intensity < peak.intensity:
            base_peak = peak
            max_intensity = peak.intensity
    return base_peak


@cython.boundscheck(False)
@cython.cdivision(True)
cpdef double correlation(cython.floating[:] x, cython.floating[:] y):
    cdef:
        size_t i, n
        double xsum, ysum, xmean, ymean
        double cov, varx, vary
    n = len(x)
    if n == 0:
        return NAN
    xsum = 0.
    ysum = 0.
    for i in range(n):
        xsum += x[i]
        ysum += y[i]
    xmean = xsum / n
    ymean = ysum / n
    cov = 0.
    varx = 0.
    vary = 0.
    for i in range(n):
        cov += (x[i] - xmean) * (y[i] - ymean)
        varx += (x[i] - xmean) ** 2
        vary += (y[i] - ymean) ** 2
    return cov / (sqrt(varx) * sqrt(vary))