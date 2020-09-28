# -*- coding: utf-8 -*-

import logging

import numpy as np

from ms_peak_picker import (
    FittedPeak, PeakSet, PeakIndex, simple_peak, is_peak)

from ms_deisotope.averagine import (isotopic_shift, neutral_mass)
from ms_deisotope.peak_set import DeconvolutedPeak, Envelope
from ms_deisotope.constants import ERROR_TOLERANCE


logger = logging.getLogger("deconvolution")
logger.addHandler(logging.NullHandler())
info = logger.info
debug = logger.debug
error = logger.error


def prepare_peaklist(peaks):
    '''Ensure ``peaks`` is a :class:`~.PeakSet` object,
    converting from other compatible types as needed. Additionally, make a deep
    copy of the peaks as signal subtraction methods will modify peaks in place.

    This function ensures that any of the following common input types are coerced
    to the appropriate type:

    1. :class:`ms_peak_picker.PeakSet` will be copied and indexed
    2. :class:`ms_peak_picker.PeakIndex` will have its peaks extracted and copied
    3. Any other *sequence* of :class:`PeakLike` objects (objects having an mz and
       intensity attribute) will be converted into a :class:`ms_peak_picker.PeakSet`
    4. Any *sequence* of :class:`tuple` or :class:`list` having at least two entries
       will be converted into a :class:`ms_peak_picker.PeakSet` with the m/z value
       of each peak being the the `p[0]` of each entry and the intensity `p[1]`. Any
       other entries will be ignored.

    Parameters
    ----------
    peaks: Sequence
        Any sequence of :class:`~.FittedPeak` objects, objects
        with ``mz`` and ``intensity`` attributes, or :class:`list` / :class:`tuple`
        objects containing paired values for ``mz`` and ``intensity``

    Returns
    -------
    :class:`~.PeakSet`
    '''
    if isinstance(peaks, PeakIndex):
        peaks = PeakSet(peaks.peaks).clone()
    else:
        peaks = tuple(peaks)
        if len(peaks) == 0:
            return PeakSet([])
        if not isinstance(peaks[0], FittedPeak):
            if is_peak(peaks[0]):
                peaks = [simple_peak(p.mz, p.intensity, 0.01) for p in peaks]
            elif isinstance(peaks[0], (list, tuple)):
                peaks = [simple_peak(p[0], p[1], 0.01) for p in peaks]
            else:
                raise TypeError("Cannot convert peaks into a PeakSet")

        peaks = PeakSet(peaks).clone()
    peaks.reindex()
    return peaks


def from_fitted_peak(peak, charge=1):
    """Convert a :class:`~.FittedPeak` into a :class:`~.DeconvolutedPeak`
    at the specified charge state.

    Parameters
    ----------
    peak : :class:`~.FittedPeak`
        The fitted peak to use as the template
    charge : int, optional
        The charge state to use, defaults to 1+

    Returns
    -------
    :class:`~.DeconvolutedPeak`
    """
    mass = neutral_mass(peak.mz, charge)
    dpeak = DeconvolutedPeak(
        mass, peak.intensity, charge,
        peak.signal_to_noise, -1, peak.full_width_at_half_max,
        0, mass, mass, 0, Envelope([(peak.mz, peak.intensity)]),
        peak.mz, area=peak.area)
    return dpeak


def mean(numbers):
    '''quick and dirty mean calculation
    without converting to a NumPy array

    Parameters
    ----------
    numbers: Iterable

    Returns
    -------
    float
    '''
    n = 0.
    total = 0
    for x in numbers:
        n += 1.
        total += x
    return total / n


def has_previous_peak_at_charge(peak_collection, peak, charge=2, step=1, error_tolerance=ERROR_TOLERANCE):
    """Get the `step`th *preceding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_collection : DeconvoluterBase
        Peak collection to look up peaks in. Calls :meth:`has_peak`
    peak : ms_peak_picker.FittedPeak
        The peak to use as a point of reference
    charge : int, optional
        The charge state to interpolate from. Defaults to 2.
    step : int, optional
        The number of peaks along the isotopic pattern to search.

    Returns
    -------
    FittedPeak
    """
    prev = peak.mz - isotopic_shift(charge) * step
    return peak_collection.has_peak(prev, error_tolerance)


def has_successor_peak_at_charge(peak_collection, peak, charge=2, step=1, error_tolerance=ERROR_TOLERANCE):
    """Get the `step`th *succeeding* peak from `peak` in a isotopic pattern at
    charge state `charge`, or return `None` if it is missing.

    Parameters
    ----------
    peak_collection : DeconvoluterBase
        Peak collection to look up peaks in. Calls :meth:`has_peak`
    peak : ms_peak_picker.FittedPeak
        The peak to use as a point of reference
    charge : int, optional
        The charge state to interpolate from. Defaults to 2.
    step : int, optional
        The number of peaks along the isotopic pattern to search.

    Returns
    -------
    FittedPeak
    """
    nxt = peak.mz + isotopic_shift(charge) * step
    return peak_collection.has_peak(nxt, error_tolerance)


def first_peak(peaks):
    """Get the first non-placeholder peak in a list of peaks

    Parameters
    ----------
    peaks : Iterable of FittedPeak

    Returns
    -------
    FittedPeak
    """
    for peak in peaks:
        if peak.intensity > 1 and peak.mz > 1:
            return peak


def drop_placeholders(peaks):
    """Removes all placeholder peaks from an iterable of peaks

    Parameters
    ----------
    peaks : Iterable of FittedPeak

    Returns
    -------
    list
    """
    return [peak for peak in peaks if peak.mz > 1 and peak.intensity > 1]


def count_placeholders(peaks):
    """Counts the number of placeholder peaks in an iterable
    of FittedPeaks

    Parameters
    ----------
    peaks : Iterable of FittedPeak

    Returns
    -------
    int
        Number of placeholder peaks
    """
    i = 0
    for peak in peaks:
        if peak.intensity <= 1:
            i += 1
    return i


def drop_placeholders_parallel(peaks, otherpeaks):
    """Given two parallel iterables of Peak objects, `peaks` and `otherpeaks`,
    for each position that is not a placeholder in `peaks`, include that Peak object
    and its counterpart in `otherpeaks` in a pair of output lists.

    Parameters
    ----------
    peaks : Iterable of FittedPeak
        Peak collection to filter against
    otherpeaks : Iterable
        Collection of objects (Peak-like) to include based upon
        contents of `peaks`

    Returns
    -------
    list
        Filtered form of `peaks`
    list
        Filtered form of `otherpeaks`
    """
    new_peaks = []
    new_otherpeaks = []
    for i, peak in enumerate(peaks):
        peak = peaks[i]
        if peak.intensity > 1:
            new_peaks.append(peak)
            new_otherpeaks.append(otherpeaks[i])
    return new_peaks, new_otherpeaks


def quick_charge(peak_set, index, min_charge, max_charge):
    """An implementation of Hoopman's :title-reference:`QuickCharge` [1] algorithm for quickly capping charge
    state queries

    Parameters
    ----------
    peak_set : :class:`ms_peak_picker.PeakSet
        The centroided peak set to search
    index : int
        The index of the peak to start the search from
    min_charge : int
        The minimum charge state to consider
    max_charge : int
        The maximum charge state to consider

    Returns
    -------
    np.ndarray
        The list of feasible charge states

    References
    ----------
    [1] Hoopmann, M. R., Finney, G. L., MacCoss, M. J., Michael R. Hoopmann, Gregory L. Finney,
        and, MacCoss*, M. J., … MacCoss, M. J. (2007). "High-speed data reduction, feature detection
        and MS/MS spectrum quality assessment of shotgun proteomics data sets using high-resolution
        Mass Spectrometry". Analytical Chemistry, 79(15), 5620–5632. https://doi.org/10.1021/ac0700833
    """
    size = len(peak_set)
    if size == 0:
        raise IndexError("peak_set cannot be empty!")
    if index > size:
        first_index = peak_set[0].peak_count
        if (index - first_index) < size and (index - first_index) >= 0:
            index -= first_index
        else:
            raise IndexError("%d is out of bounds for peak list of size %d in quick_charge" % (index, size))
    min_intensity = peak_set[index].intensity / 4.
    charges = np.zeros(max_charge, dtype=int)
    for j in range(index + 1, size):
        if peak_set[j].intensity < min_intensity:
            continue
        diff = peak_set[j].mz - peak_set[index].mz
        if diff > 1.1:
            break
        raw_charge = 1 / diff
        charge = int(raw_charge + 0.5)
        remain = raw_charge - int(raw_charge)
        if 0.2 < remain < 0.8:
            continue
        if charge < min_charge or charge > max_charge:
            continue
        charges[charge] = 1
    if not np.any(charges):
        return np.array([], dtype=int)
    for j in range(index - 1, -1, -1):
        diff = peak_set[index].mz - peak_set[j].mz
        if diff > 1.1:
            break
        raw_charge = 1 / diff
        charge = int(raw_charge + 0.5)
        remain = raw_charge - int(raw_charge)
        if 0.2 < remain < 0.8:
            continue
        if charge < min_charge or charge > max_charge:
            continue
        charges[charge] = 1
    return np.where(charges)[0]


try:
    _quick_charge = quick_charge
    from ms_deisotope._c.deconvoluter_base import quick_charge
except ImportError:
    pass


def charge_range_(lo, hi, step=None):
    """Generate a sequence of successive charge states, automatically deducing
    the polarity of the charge range to infer the step.

    Parameters
    ----------
    lo : int
        The smallest magnitude charge (positive or negative)
    hi : int
        The largest magnitude charge (positive or negative)
    step : int, optional
        The step size and sign to incrementally alter the charge by. Automatically
        inferred if omitted.

    Yields
    ------
    int
    """
    sign = -1 if lo < 0 else 1
    abs_lo, abs_hi = abs(lo), abs(hi)
    upper = max(abs_lo, abs_hi)
    lower = min(abs_lo, abs_hi)

    for c in range(upper, lower - 1, -1):
        yield c * sign


class ChargeIterator(object):
    def __init__(self, lo, hi):
        self.set_bounds(lo, hi)
        self.make_sequence()

    def set_bounds(self, lo, hi):
        self.sign = -1 if lo < 0 else 1
        abs_lo, abs_hi = abs(lo), abs(hi)
        if abs_lo < abs_hi:
            self.lower = abs_lo
            self.upper = abs_hi
        else:
            self.lower = abs_hi
            self.upper = abs_lo
        self.size = self.upper - self.lower + 1

    def make_sequence(self):
        self.index = 0
        self.size = self.upper - self.lower + 1
        self.values = [self.sign * (self.upper - i) for i in range(self.size)]

    def __len__(self):
        return self.size

    def reset(self):
        self.index = 0

    def sequence_from_quickcharge(self, peak_set, peak):
        charges = quick_charge(peak_set, peak.peak_count,
                               abs(self.lower), abs(self.upper))
        n = charges.shape[0]
        self.index = 0
        if n == 0:
            self.size = 1
            self.values = [1 * self.sign]
        elif charges[0] != 1:
            self.size = n + 1
            self.values = [1 * self.sign] + list(self.sign * charges)
        else:
            self.size = n
            self.values = self.sign * charges

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == self.size:
            raise StopIteration()
        value = self.values[self.index]
        self.index += 1
        return value

    def next(self):
        return self.__next__()
