# cython: embedsignature=True

cimport cython
from cpython.list cimport PyList_Append, PyList_GET_ITEM

from libc.math cimport sqrt, exp, pi

import numpy as np
cimport numpy as np

from ms_deisotope._c.peak_set cimport (
    Envelope, EnvelopePair, DeconvolutedPeak, DeconvolutedPeakSetIndexed)
from ms_deisotope._c.averagine cimport neutral_mass


np.import_array()


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
                current_envelope = []
                n_members = 0
        else:
            PyList_Append(current_envelope, EnvelopePair._create(a, b))
            n_members += 1
    current_envelope_tuple = tuple(current_envelope)
    envelope_list.append(Envelope._create(current_envelope_tuple))
    return envelope_list


@cython.boundscheck(False)
cpdef DeconvolutedPeakSetIndexed deserialize_deconvoluted_peak_set(dict scan_dict):
    cdef:
        list envelopes, peaks
        np.ndarray[np.float32_t] mz_array
        np.ndarray[np.float32_t] intensity_array
        np.ndarray[np.int8_t] charge_array
        np.ndarray[np.float32_t] score_array
        size_t n, i
        double mz, peak_neutral_mass
        int charge
        DeconvolutedPeak peak
        DeconvolutedPeakSetIndexed peak_set

    peaks = []
    envelopes = decode_envelopes(scan_dict["isotopic envelopes array"])
    mz_array = scan_dict['m/z array'].astype(np.float32)
    intensity_array = scan_dict['intensity array'].astype(np.float32)
    charge_array = scan_dict['charge array'].astype(np.int8)
    score_array = scan_dict['deconvolution score array'].astype(np.float32)
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
            <Envelope>PyList_GET_ITEM(envelopes, i))
        peaks.append(peak)
    peak_set = DeconvolutedPeakSetIndexed(peaks)
    peak_set.reindex()
    return peak_set
