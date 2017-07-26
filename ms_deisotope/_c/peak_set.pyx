# cython: embedsignature=True

import operator

cimport cython
from cpython.tuple cimport PyTuple_GET_ITEM, PyTuple_GetItem, PyTuple_GetSlice, PyTuple_GET_SIZE

from ms_deisotope._c.averagine cimport mass_charge_ratio
from ms_peak_picker._c.peak_set cimport PeakBase

@cython.cdivision
cdef double ppm_error(double x, double y):
    return (x - y) / y


cdef:
    object neutral_mass_getter
    object mz_getter


neutral_mass_getter = operator.attrgetter("neutral_mass")
mz_getter = operator.attrgetter('mz')


@cython.freelist(100000)
cdef class _Index(object):
    """
    Stores the ordered index of an object under sort by `mz` and
    sort by `neutral_mass`

    Attributes
    ----------
    mz : int
        Index of this object (or its parent), ordered by mz
    neutral_mass : int
        Index of this object (or its parent), ordered by neutral mass
    """
    def __init__(self, neutral_mass=None, mz=None):
        self.neutral_mass = neutral_mass
        self.mz = mz

    def __repr__(self):
        return "%d|%d" % (self.neutral_mass, self.mz)

    def clone(self):
        return self.__class__(self.neutral_mass, self.mz)

    def __reduce__(self):
        return _Index, (self.neutral_mass, self.mz)


@cython.freelist(1000000)
cdef class EnvelopePair:

    def __init__(self, mz, intensity):
        self.mz = mz
        self.intensity = intensity
    
    def __getitem__(self, i):
        if i == 0:
            return self.mz
        elif i == 1:
            return self.intensity
        else:
            raise IndexError(i)

    def __iter__(self):
        yield self.mz
        yield self.intensity

    def __reduce__(self):
        return EnvelopePair, (self.mz, self.intensity,)

    cpdef bint _eq(self, EnvelopePair other):
        return (abs(self.mz - other.mz) < 1e-5) and (
            abs(self.intensity - other.intensity) < 1e-5)

    def __richcmp__(self, other, int code):
        if code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))


    @staticmethod
    cdef EnvelopePair _create(double mz, double intensity):
        cdef:
            EnvelopePair inst
        inst = EnvelopePair.__new__(EnvelopePair)
        inst.mz = mz
        inst.intensity = intensity
        return inst


@cython.freelist(100000)
cdef class Envelope(object):
    """
    Represents a sequence of (mz, intensity) pairs which store peak positions
    that were matched. Since these peaks may later be mutated by subtraction,
    this structure stores copies of the numbers

    Supports the Sequence protocol, and is read-only.

    Attributes
    ----------
    pairs : tuple of tuple of (float, float)
        list of (mz, intensity) pairs
    """
    def __init__(self, pairs):
        self.pairs = tuple(map(lambda x: EnvelopePair(*x), pairs))

    def __getitem__(self, i):
        return self.pairs[i]

    def __iter__(self):
        return iter(self.pairs)

    def __repr__(self):
        return "[%s]" % (', '.join("(%0.4f, %0.2f)" % tuple(t) for t in self))

    def __len__(self):
        return len(self.pairs)

    cdef inline size_t get_size(self):
        return PyTuple_GET_SIZE(self.pairs)

    cdef inline EnvelopePair getitem(self, size_t i):
        return <EnvelopePair>PyTuple_GetItem(self.pairs, i)

    cpdef bint _eq(self, Envelope other):
        cdef:
            size_t i, n
            EnvelopePair a, b

        n = self.get_size()
        if n != other.get_size():
            return False

        i = 0
        for i in range(n):
            a = self.getitem(i)
            b = other.getitem(i)
            if a != b:
                return False
        return True

    def __richcmp__(self, other, int code):
        if code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))

    def clone(self):
        return self.__class__(self)

    def __reduce__(self):
        return Envelope, (self.pairs,)

    @staticmethod
    cdef Envelope _create(tuple pairs):
        cdef:
            Envelope inst
        inst = Envelope.__new__(Envelope)
        inst.pairs = pairs
        return inst


# @cython.freelist(100000)
cdef class DeconvolutedPeak(PeakBase):
    """
    Represent a single deconvoluted peak which represents an aggregated isotopic
    pattern collapsed to its monoisotopic peak, with a known charge state

    Attributes
    ----------
    a_to_a2_ratio : float
        Ratio of intensities of A peak to A+2 peak
    average_mass : float
        The averaged neutral mass of the composition, the weighted average
        of the envelope peaks.
    charge : int
        The signed charge state of the isotopic pattern
    envelope : Envelope
        The sequence of (mz, intensity) pairs which map to this peak
    full_width_at_half_max : float
        The averaged full width at half max of this isotopic pattern
    index : _Index
        The position of this peak in different orderings
    intensity : float
        The summed height of the peaks in the isotopic pattern this peak captures
    most_abundant_mass : float
        The neutral mass of the most abundant peak in the isotopic pattern this peak captures
    mz : float
        The mass-charge-ratio of the monoisotopic peak
    neutral_mass : float
        The neutral mass of the monoisotopic peak
    score : float
        An assigned value describing the quality of this peak's fit. The semantics of this score
        depends upon the scoring function
    signal_to_noise : float
        The average signal-to-noise ratio of the peaks in the isotopic pattern this peak captures
    """
    def __init__(self, neutral_mass, intensity, charge, signal_to_noise, index, full_width_at_half_max,
                 a_to_a2_ratio=0, most_abundant_mass=0, average_mass=0, score=0,
                 envelope=(), mz=0, fit=None, chosen_for_msms=False, area=0):
        if index is None:
            index = _Index()
        elif index == -1:
            index = _Index(0, 0)
        self.neutral_mass = neutral_mass
        self.intensity = intensity
        self.signal_to_noise = signal_to_noise
        self.index = index
        self.full_width_at_half_max = full_width_at_half_max
        self.charge = charge
        self.a_to_a2_ratio = a_to_a2_ratio
        self.most_abundant_mass = most_abundant_mass
        self.average_mass = average_mass
        self.score = score
        self.envelope = Envelope(envelope)
        self.mz = mz or mass_charge_ratio(self.neutral_mass, self.charge)
        self.fit = fit
        self.chosen_for_msms = chosen_for_msms
        self.area = area

    property index:
        def __get__(self):
            return self._index

        def __set__(self, value):
            if isinstance(value, _Index):
                self._index = value
            elif value >= 0:
                self._index = _Index(value, value)
            else:
                self._index = _Index(0, 0)

    def __hash__(self):
        return hash((self.mz, self.intensity, self.charge))

    cpdef PeakBase clone(self):
        return DeconvolutedPeak(self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
                                self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
                                self.most_abundant_mass, self.average_mass, self.score,
                                self.envelope, self.mz, self.fit, self.chosen_for_msms, self.area)

    def __reduce__(self):
        return DeconvolutedPeak, (self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
                                  self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
                                  self.most_abundant_mass, self.average_mass, self.score,
                                  self.envelope, self.mz, self.fit, self.chosen_for_msms, self.area)

    cpdef bint _eq(self, DeconvolutedPeak other):
        return (abs(self.neutral_mass - other.neutral_mass) < 1e-5) and (
            abs(self.intensity - other.intensity) < 1e-5)

    def __richcmp__(self, other, int code):
        if code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))

    def __repr__(self):
        return ("DeconvolutedPeak(a_to_a2_ratio={self.a_to_a2_ratio}, area={self.area}, "
            "average_mass={self.average_mass}, charge={self.charge}, chosen_for_msms={self.chosen_for_msms}, "
            "envelope={self.envelope}, full_width_at_half_max={self.full_width_at_half_max}, index={self.index}, "
            "intensity={self.intensity}, most_abundant_mass={self.most_abundant_mass}, mz={self.mz}, "
            "neutral_mass={self.neutral_mass}, score={self.score}, signal_to_noise={self.signal_to_noise})").format(self=self)

    
    @staticmethod
    cdef DeconvolutedPeak _create_simple(double neutral_mass, double intensity, int charge,
                                         double score, double mz, Envelope envelope):
        cdef:
            DeconvolutedPeak inst
        inst = DeconvolutedPeak.__new__(DeconvolutedPeak)
        inst.neutral_mass = neutral_mass
        inst.intensity = intensity
        inst.charge = charge
        inst.score = score
        inst.signal_to_noise = score
        inst.mz = mz
        inst.envelope = envelope
        inst.index = _Index(0, 0)
        inst.full_width_at_half_max = 0
        inst.a_to_a2_ratio = 0
        inst.most_abundant_mass = 0
        inst.average_mass = 0
        inst.area = 0

        return inst


cdef class DeconvolutedPeakSolution(DeconvolutedPeak):
    """
    Extends :class:`DeconvolutedPeak` to also include a reference to
    the :class:`IsotopicFitRecord` instance it is derived from, and optionally an
    object which the fit derives from.

    Attributes
    ----------
    fit : IsotpicFitRecord
        The isotopic pattern fit used to construct this peak
    solution : object
        Representation of the "source" of the isotopic fit
    """
    def __init__(self, solution, fit, *args, **kwargs):
        self.solution = solution
        self.fit = fit
        super(DeconvolutedPeakSolution, self).__init__(*args, **kwargs)

    cpdef PeakBase clone(self):
        return DeconvolutedPeakSolution(
            self.solution, self.fit, self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.chosen_for_msms, self.area)

    def __reduce__(self):
        return DeconvolutedPeakSolution, (
            self.solution, self.fit, self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.chosen_for_msms, self.area)

    def __iter__(self):
        yield self.solution
        yield self
        yield self.fit


cdef class DeconvolutedPeakSet:
    """
    Represents a collection of :class:`DeconvolutedPeak` instances under multiple orderings.

    Attributes
    ----------
    peaks : tuple of DeconvolutedPeak
        Collection of peaks ordered by `neutral_mass`
    _mz_ordered: tuple of DeconvolutedPeak
        Collection of peaks ordered by `mz`
    """
    def __init__(self, peaks):
        self.peaks = tuple(peaks)
        self._mz_ordered = None

    def reindex(self):
        """
        Updates the :attr:`index` of each peak in `self` and updates the
        sorted order.

        Returns
        -------
        self: DeconvolutedPeakSet
        """
        self._reindex()

    def _reindex(self):
        """
        Updates the :attr:`index` of each peak in `self` and updates the
        sorted order.

        Returns
        -------
        self: DeconvolutedPeakSet
        """
        cdef:
            size_t i, n
            DeconvolutedPeak peak
        self.peaks = tuple(sorted(self.peaks, key=neutral_mass_getter))
        self._mz_ordered = tuple(sorted(self.peaks, key=mz_getter))
        n = PyTuple_GET_SIZE(self.peaks)
        i = 0
        for i in range(n):
            peak = self.getitem(i)
            peak._index = _Index(0, 0)
            peak._index.neutral_mass = i
        i = 0
        for i in range(n):
            peak = <DeconvolutedPeak>PyTuple_GET_ITEM(self._mz_ordered, i)
            peak._index.mz = i
        return self

    def __iter__(self):
        return iter(self.peaks)

    def __len__(self):
        return PyTuple_GET_SIZE(self.peaks)

    def __repr__(self):
        return "<DeconvolutedPeakSet %d Peaks>" % (len(self))

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.__class__(tuple(p.clone() for p in self.peaks[item]))
        return self.peaks[item]

    def __iter__(self):
        return iter(self.peaks)

    def clone(self):
        return self.__class__(tuple(p.clone() for p in self))

    def __reduce__(self):
        return DeconvolutedPeakSet, (self.peaks, )

    cdef DeconvolutedPeak _has_peak(self, double neutral_mass, double error_tolerance=1e-5, bint use_mz=False):
        if use_mz:
            return binary_search_mz(self._mz_ordered, neutral_mass, error_tolerance)
        else:
            return binary_search_neutral_mass(self.peaks, neutral_mass, error_tolerance)

    cpdef DeconvolutedPeak has_peak(self, double neutral_mass, double error_tolerance=1e-5, bint use_mz=False):
        return self._has_peak(neutral_mass, error_tolerance, use_mz)

    cdef DeconvolutedPeak getitem(self, size_t i):
        return <DeconvolutedPeak>PyTuple_GET_ITEM(self.peaks, i)

    def all_peaks_for(self, double neutral_mass, double tolerance=1e-5):
        cdef:
            double lo, hi, lo_err, hi_err
            int lo_ix, hi_ix
            DeconvolutedPeak lo_peak
            DeconvolutedPeak hi_peak
        lo = neutral_mass - neutral_mass * tolerance
        hi = neutral_mass + neutral_mass * tolerance
        lo_peak = binary_search_nearest_neutral_mass(self.peaks, lo, &lo_err)
        hi_peak = binary_search_nearest_neutral_mass(self.peaks, hi, &hi_err)
        lo_ix = lo_peak._index.neutral_mass
        if lo_ix < PyTuple_GET_SIZE(self.peaks) and abs(
                ppm_error(lo_peak.neutral_mass, neutral_mass)) > tolerance:
            lo_ix += 1
        hi_ix = hi_peak._index.neutral_mass + 1
        if hi_ix != 0 and abs(
                ppm_error(hi_peak.neutral_mass, neutral_mass)) > tolerance:
            hi_ix -= 1
        return self[lo_ix:hi_ix]

    def get_nearest_peak(self, double neutral_mass):
        cdef:
            DeconvolutedPeak peak
            double errout
        peak = binary_search_nearest_neutral_mass(self.peaks, neutral_mass, &errout)
        return peak, errout

    def between(self, m1, m2, tolerance=1e-5, use_mz=False):
        acc = []
        collecting = False
        if not use_mz:
            getter = operator.attrgetter("neutral_mass")
            iterator = self.peaks
        else:
            getter = operator.attrgetter("mz")
            iterator = self._mz_ordered
        for peak in iterator:
            if not collecting and getter(peak) >= m1:
                collecting = True
            elif collecting and getter(peak) > m2:
                break

            if collecting:
                acc.append(peak.clone())

        return self.__class__(acc)._reindex()

cdef double INF
INF = float('inf')


cdef DeconvolutedPeak _sweep_solution_neutral_mass(tuple array, double value, size_t mid, double tolerance):
    cdef:
        size_t best_index, i, n
        double best_error, abs_error
        DeconvolutedPeak target

    best_index = mid
    best_error = INF
    n = PyTuple_GET_SIZE(array)

    i = 0
    while mid - i >= 0 and i <= mid:
        target = <DeconvolutedPeak>array[mid - i]
        abs_error = abs(ppm_error(value, target.neutral_mass))
        if abs_error < tolerance:
            if abs_error < best_error:
                best_index = mid - i
                best_error = abs_error
        else:
            break
        i += 1
    i = 1
    while (mid + i) < (n - 1):
        target = <DeconvolutedPeak>array[mid + i]
        abs_error = abs(ppm_error(value, target.neutral_mass))
        if abs_error < tolerance:
            if abs_error < best_error:
                best_index = mid + i
                best_error = abs_error
        else:
            break
        i += 1
    if best_error == INF:
        return None
    else:
        return array[best_index]


cdef DeconvolutedPeak binary_search_neutral_mass(tuple peak_set, double neutral_mass, double error_tolerance):
    cdef:
        size_t lo, hi, mid, i, j
        double err, found_mass
        DeconvolutedPeak found_peak

    lo = 0
    hi = PyTuple_GET_SIZE(peak_set)

    while hi != lo:
        mid = (hi + lo) / 2
        found_peak = <DeconvolutedPeak>peak_set[mid]
        found_mass = found_peak.neutral_mass

        if abs(ppm_error(found_mass, neutral_mass)) < error_tolerance:
            found_peak = _sweep_solution_neutral_mass(peak_set, neutral_mass, mid, error_tolerance)
            return found_peak
        elif hi - lo == 1:
            return None
        elif found_mass > neutral_mass:
            hi = mid
        else:
            lo = mid


cdef DeconvolutedPeak _sweep_solution_mz(tuple array, double value, size_t mid, double tolerance):
    cdef:
        size_t best_index, i, n
        double best_error, abs_error
        DeconvolutedPeak target

    best_index = mid
    best_error = INF
    n = PyTuple_GET_SIZE(array)

    i = 0
    while mid - i >= 0 and i <= mid:
        target = <DeconvolutedPeak>array[mid - i]
        abs_error = abs(ppm_error(value, target.mz))
        if abs_error < tolerance:
            if abs_error < best_error:
                best_index = mid - i
                best_error = abs_error
        else:
            break
        i += 1
    i = 1
    while (mid + i) < (n - 1):
        target = <DeconvolutedPeak>array[mid + i]
        abs_error = abs(ppm_error(value, target.mz))
        if abs_error < tolerance:
            if abs_error < best_error:
                best_index = mid + i
                best_error = abs_error
        else:
            break
        i += 1
    if best_error == INF:
        return None
    else:
        return array[best_index]


cdef DeconvolutedPeak binary_search_mz(tuple peak_set, double mz, double error_tolerance):
    cdef:
        size_t lo, hi, mid, i, j
        double err, found_mass
        DeconvolutedPeak found_peak

    lo = 0
    hi = len(peak_set)

    while hi != lo:
        mid = (hi + lo) / 2
        found_peak = <DeconvolutedPeak>peak_set[mid]
        found_mass = found_peak.mz

        if abs(ppm_error(found_mass, mz)) < error_tolerance:
            found_peak = _sweep_solution_mz(peak_set, mz, mid, error_tolerance)
            return found_peak
        elif hi - lo == 1:
            return None
        elif found_mass > mz:
            hi = mid
        else:
            lo = mid


cdef DeconvolutedPeak _sweep_nearest_match_neutral_mass(tuple array, double value, size_t lo, size_t hi, double* errout):
    cdef:
        size_t i
        size_t best_index
        double best_error, err
        double v

    best_error = float('inf')
    best_index = -1
    for i in range(hi - lo):
        i += lo
        v = array[i].neutral_mass
        err = abs(v - value)
        if err < best_error:
            best_error = err
            best_index = i
    errout[0] = best_error
    return array[best_index]


cdef DeconvolutedPeak binary_search_nearest_neutral_mass(tuple peak_set, double neutral_mass, double* errout):
    cdef:
        size_t lo, hi, mid, i, j
        double err, found_mass
        DeconvolutedPeak found_peak

    lo = 0
    hi = len(peak_set)
    while hi != lo:
        mid = (hi + lo) / 2
        found_peak = <DeconvolutedPeak>peak_set[mid]
        found_mass = found_peak.neutral_mass

        if abs(ppm_error(found_mass, neutral_mass)) < 1:
            found_peak = _sweep_nearest_match_neutral_mass(peak_set, neutral_mass, lo, hi, errout)
            return found_peak
        elif hi - lo == 1:
            errout[0] = found_mass - neutral_mass
            return found_peak
        elif found_mass > neutral_mass:
            hi = mid
        else:
            lo = mid


def convert(self):
    return DeconvolutedPeak(
        self.neutral_mass, self.intensity, self.charge,
        self.signal_to_noise, -1, self.full_width_at_half_max,
        self.a_to_a2_ratio, self.most_abundant_mass, self.average_mass,
        self.score, Envelope(self.envelope), self.mz, None, self.chosen_for_msms,
        self.area)


cimport numpy as np


def decode_envelopes(_array):
    cdef:
        np.ndarray array
        list envelope_list, current_envelope
        size_t i, n
        object a, b
    envelope_list = []
    current_envelope = []
    array = _array
    i = 0
    n = len(array)
    while i < n:
        a = array[i]
        b = array[i + 1]
        i += 2
        if a == 0 and b == 0:
            if current_envelope is not None:
                if current_envelope:
                    envelope_list.append(Envelope(current_envelope))
                current_envelope = []
        else:
            current_envelope.append(EnvelopePair(a, b))
    envelope_list.append(Envelope(current_envelope))
    return envelope_list
