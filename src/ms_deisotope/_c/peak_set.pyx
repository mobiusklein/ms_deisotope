# cython: embedsignature=True

import operator

cimport cython
from cpython.tuple cimport (PyTuple_GET_ITEM, PyTuple_GetItem, PyTuple_GetSlice,
PyTuple_GET_SIZE, PyTuple_New, PyTuple_SetItem)
from cpython.list cimport PyList_Append, PyList_AsTuple
from cpython.ref cimport Py_INCREF
from cpython.object cimport PyObject
from cpython.exc cimport PyErr_SetString
# from cpython.bytearray cimport PyByteArray_FromStringAndSize

# to support older versions of Cython which do not include this header
cdef extern from "Python.h":
    bytearray PyByteArray_FromStringAndSize(char *string, Py_ssize_t len)


from libc.string cimport memcpy
from libc.math cimport fabs
from libc.stdlib cimport malloc, free

from ms_deisotope._c.averagine cimport mass_charge_ratio
from ms_peak_picker._c.peak_set cimport PeakBase

cimport numpy as np

np.import_array()


@cython.cdivision
cdef double ppm_error(double x, double y):
    return (x - y) / y


cdef:
    object neutral_mass_getter
    object mz_getter

cdef double epsilon = 1e-3

neutral_mass_getter = operator.attrgetter("neutral_mass")
mz_getter = operator.attrgetter('mz')


cpdef _Index_reconstructor(object neutral_mass, object mz=None):
    if isinstance(neutral_mass, int):
        return _Index._create(neutral_mass, mz)
    elif isinstance(neutral_mass, bytearray):
        return _Index_from_bytes(neutral_mass)
    else:
        raise TypeError(type(neutral_mass))


cpdef _Index_to_bytes(_Index index):
    cdef:
        size_t j, m
        char* buff
        np.uint64_t n
        bytearray out

    m = sizeof(np.uint64_t) * 2
    buff = <char*>malloc(sizeof(char) * m)
    j = 0
    n = index.neutral_mass
    memcpy(&buff[j], &n, sizeof(np.uint64_t))
    j += sizeof(np.uint64_t)
    n = index.mz
    memcpy(&buff[j], &n, sizeof(np.uint64_t))
    j += sizeof(np.uint64_t)
    out = PyByteArray_FromStringAndSize(<char*>buff, m)
    free(buff)
    return out


cpdef _Index_from_bytes(bytearray view):
    cdef:
        size_t j
        np.uint64_t val
        char* buff
        _Index inst
    inst = _Index._create(0, 0)
    j = 0
    buff = view
    memcpy(&val, &buff[j], sizeof(np.uint64_t))
    j += sizeof(np.uint64_t)
    inst.neutral_mass = val
    memcpy(&val, &buff[j], sizeof(np.uint64_t))
    j += sizeof(np.uint64_t)
    inst.mz = val
    return inst

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
    def __init__(self, neutral_mass=0, mz=0):
        self.neutral_mass = neutral_mass
        self.mz = mz

    def __repr__(self):
        return "%d|%d" % (self.neutral_mass, self.mz)

    def clone(self):
        return _Index._create(self.neutral_mass, self.mz)

    def __reduce__(self):
        return _Index_reconstructor, (self.neutral_mass, self.mz)

    @staticmethod
    cdef _Index _create(size_t neutral_mass, size_t mz):
        cdef:
            _Index inst
        inst = _Index.__new__(_Index)
        inst.neutral_mass = neutral_mass
        inst.mz = mz
        return inst


@cython.freelist(1000000)
cdef class EnvelopePair(object):

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
        return (abs(self.mz - other.mz) < epsilon) and (
            abs(self.intensity - other.intensity) < epsilon)

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

    def __repr__(self):
        return "(%0.4f, %0.2f)" % (self.mz, self.intensity)


@cython.boundscheck(False)
cpdef bytearray Envelope_to_bytes(Envelope self):
    cdef:
        size_t i, j, m
        EnvelopePair pair
        char* buff
        np.uint32_t n
        np.float64_t val
        bytearray out

    n = self.get_size()
    m = sizeof(np.uint32_t) + sizeof(np.float64_t) * n * 2
    buff = <char*>malloc(sizeof(char) * m)
    j = 0
    memcpy(&buff[j], &n, sizeof(np.uint32_t))
    j += sizeof(np.uint32_t)
    for i in range(n):
        pair = self.getitem(i)
        val = pair.mz
        memcpy(&buff[j], &val, sizeof(np.float64_t))
        j += sizeof(np.float64_t)
        val = pair.intensity
        memcpy(&buff[j], &val, sizeof(np.float64_t))
        j += sizeof(np.float64_t)
    out = PyByteArray_FromStringAndSize(<char*>buff, m)
    free(buff)
    return out


cpdef Envelope_from_bytes(bytearray view):
    cdef:
        size_t i, j, m
        EnvelopePair pair
        np.uint32_t n
        np.float64_t val
        char* buff
        tuple result
    j = 0
    buff = view
    memcpy(&n, &buff[j], sizeof(np.uint32_t))
    j += sizeof(np.uint32_t)
    result = PyTuple_New(n)
    for i in range(n):
        pair = EnvelopePair._create(0, 0)
        memcpy(&pair.mz, &buff[j], sizeof(np.float64_t))
        j += sizeof(np.float64_t)
        memcpy(&pair.intensity, &buff[j], sizeof(np.float64_t))
        j += sizeof(np.float64_t)
        Py_INCREF(pair)
        PyTuple_SetItem(result, i, pair)
    return Envelope._create(result)

@cython.boundscheck(False)
cpdef Envelope_reconstructor(object value):
    if isinstance(value, tuple):
        return Envelope._create(value)
    elif isinstance(value, bytearray):
        return Envelope_from_bytes(value)
    else:
        raise TypeError(type(value))


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

    cpdef Envelope clone(self):
        cdef:
            tuple duplicate
            size_t i, n
            EnvelopePair pair, dup_pair
        n = self.get_size()
        duplicate = PyTuple_New(n)
        for i in range(n):
            pair = self.getitem(i)
            dup_pair = EnvelopePair._create(pair.mz, pair.intensity)
            Py_INCREF(dup_pair)
            PyTuple_SetItem(duplicate, i, dup_pair)
        return Envelope._create(duplicate)

    def __reduce__(self):
        return (
            Envelope_reconstructor,
            (Envelope_to_bytes(self),))
            # (self.pairs, ))

    @staticmethod
    cdef Envelope _create(tuple pairs):
        cdef:
            Envelope inst
        inst = Envelope.__new__(Envelope)
        inst.pairs = pairs
        return inst


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
            index = _Index._create(0, 0)
        elif index == -1:
            index = _Index._create(0, 0)
        if envelope is None:
            envelope = ()
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
        if not isinstance(envelope, Envelope):
            self.envelope = Envelope(envelope)
        else:
            self.envelope = envelope
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
                self._index = _Index._create(value, value)
            else:
                self._index = _Index._create(0, 0)

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
        return (abs(self.neutral_mass - other.neutral_mass) < epsilon) and (
            abs(self.intensity - other.intensity) < epsilon) and (
            self.charge == other.charge) and (
            abs(self.score - other.score) < epsilon)

    def __richcmp__(self, other, int code):
        if code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))

    def __repr__(self):
        return ("{self.__class__.__name__}(a_to_a2_ratio={self.a_to_a2_ratio}, area={self.area}, "
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
        inst.index = _Index._create(0, 0)
        inst.full_width_at_half_max = 0
        inst.a_to_a2_ratio = 0
        inst.most_abundant_mass = 0
        inst.average_mass = 0
        inst.area = 0

        return inst


cdef class IonMobilityDeconvolutedPeak(DeconvolutedPeak):
    def __init__(self, neutral_mass, intensity, charge, signal_to_noise, index, full_width_at_half_max,
                 a_to_a2_ratio=0, most_abundant_mass=0, average_mass=0, score=0,
                 envelope=(), mz=0, fit=None, chosen_for_msms=False, area=0, drift_time=0):
        super(IonMobilityDeconvolutedPeak, self).__init__(
            neutral_mass, intensity, charge, signal_to_noise, index, full_width_at_half_max,
            a_to_a2_ratio, most_abundant_mass, average_mass, score,
            envelope, mz, fit, chosen_for_msms, area)
        self.drift_time = drift_time

    cpdef PeakBase clone(self):
        return IonMobilityDeconvolutedPeak(
            self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.fit, self.chosen_for_msms, self.area,
            self.drift_time)

    def __reduce__(self):
        return IonMobilityDeconvolutedPeak, (
            self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.fit, self.chosen_for_msms, self.area,
            self.drift_time)

    cpdef bint _eq(self, DeconvolutedPeak other):
        return (abs(self.neutral_mass - other.neutral_mass) < epsilon) and (
            abs(self.intensity - other.intensity) < epsilon) and (
            self.charge == other.charge) and (
            abs(self.score - other.score) < epsilon) and (
            abs(self.drift_time - other.drift_time) < epsilon)

    def __repr__(self):
        return ("{self.__class__.__name__}(a_to_a2_ratio={self.a_to_a2_ratio}, area={self.area}, "
            "average_mass={self.average_mass}, charge={self.charge}, chosen_for_msms={self.chosen_for_msms}, "
            "envelope={self.envelope}, full_width_at_half_max={self.full_width_at_half_max}, index={self.index}, "
            "intensity={self.intensity}, most_abundant_mass={self.most_abundant_mass}, mz={self.mz}, "
            "neutral_mass={self.neutral_mass}, score={self.score}, signal_to_noise={self.signal_to_noise},"
            "drift_time={self.drift_time})").format(self=self)


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
        super(DeconvolutedPeakSolution, self).__init__(*args, **kwargs)
        self.fit = fit

    cpdef PeakBase clone(self):
        return self.__class__(
            self.solution, self.fit, self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.chosen_for_msms, self.area)

    def __reduce__(self):
        return self.__class__, (
            self.solution, self.fit, self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.chosen_for_msms, self.area)

    def __iter__(self):
        yield self.solution
        yield self
        yield self.fit

    @classmethod
    def from_peak(cls, peak, solution):
        inst = cls(
            solution, peak.fit, peak.neutral_mass, peak.intensity,
            peak.charge, peak.signal_to_noise, peak.index, peak.full_width_at_half_max,
            peak.a_to_a2_ratio, peak.most_abundant_mass, peak.average_mass,
            peak.score, peak.envelope, peak.mz, peak.chosen_for_msms, peak.area)
        return inst


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
        self.indexed = False

    cpdef reindex(self):
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
            peak._index = _Index._create(0, 0)
            peak._index.neutral_mass = i
        i = 0
        for i in range(n):
            peak = <DeconvolutedPeak>PyTuple_GET_ITEM(self._mz_ordered, i)
            peak._index.mz = i
        self.indexed = True
        return self

    def __iter__(self):
        return iter(self.peaks)

    def __len__(self):
        return self.get_size()

    cdef size_t get_size(self):
        return PyTuple_GET_SIZE(self.peaks)

    def __repr__(self):
        return "<DeconvolutedPeakSet %d Peaks>" % (len(self))

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.__class__(tuple(p.clone() for p in self.peaks[item]))
        return self.peaks[item]

    def __eq__(self, other):
        try:
            return self.peaks == other.peaks
        except AttributeError:
            return tuple(self) == tuple(other)

    def __ne__(self, other):
        try:
            return self.peaks != other.peaks
        except AttributeError:
            return tuple(self) != tuple(other)

    def __iter__(self):
        return iter(self.peaks)

    cpdef DeconvolutedPeakSet clone(self):
        return self.copy()

    cpdef DeconvolutedPeakSet copy(self):
        cdef:
            tuple duplicate_peaks
            DeconvolutedPeak peak
            size_t i, n
        n = self.get_size()
        duplicate_peaks = PyTuple_New(n)
        for i in range(n):
            peak = self.getitem(i)
            peak = peak.clone()
            Py_INCREF(peak)
            PyTuple_SetItem(duplicate_peaks, i, peak)
        return self.__class__(duplicate_peaks)

    def __reduce__(self):
        return self.__class__, (self.peaks, )

    cdef DeconvolutedPeak _has_peak(self, double neutral_mass, double error_tolerance=1e-5, bint use_mz=False):
        """Find the peak that best matches ``neutral_mass`` within ``error_tolerance`` mass accuracy ppm.

        If ``use_mz`` is True, instead of matching neutral masses, match peaks using m/z instead.

        Parameters
        ----------
        neutral_mass: double
            The mass to search for
        error_tolerance: double
            The PPM error tolerance to apply
        use_mz: bool
            Whether to search using m/z instead of neutral mass

        Returns
        -------
        DeconvolutedPeak
            The found peak, or None if no peak is found
        """
        if use_mz:
            return binary_search_mz(self._mz_ordered, neutral_mass, error_tolerance)
        else:
            return binary_search_neutral_mass(self.peaks, neutral_mass, error_tolerance)

    cpdef DeconvolutedPeak has_peak(self, double neutral_mass, double error_tolerance=1e-5, bint use_mz=False):
        """Find the peak that best matches ``neutral_mass`` within ``error_tolerance`` mass accuracy ppm.

        If ``use_mz`` is True, instead of matching neutral masses, match peaks using m/z instead.

        If :attr:`indexed` is not :const:`True`, then :meth:`reindex` will be called.

        Parameters
        ----------
        neutral_mass: double
            The mass to search for
        error_tolerance: double
            The PPM error tolerance to apply
        use_mz: bool
            Whether to search using m/z instead of neutral mass

        Returns
        -------
        DeconvolutedPeak
            The found peak, or None if no peak is found
        """
        if not self.indexed:
            self.reindex()
        return self._has_peak(neutral_mass, error_tolerance, use_mz)

    cdef DeconvolutedPeak getitem(self, size_t i):
        return <DeconvolutedPeak>PyTuple_GET_ITEM(self.peaks, i)

    cdef tuple getslice(self, size_t start, size_t end):
        return <tuple>PyTuple_GetSlice(self.peaks, start, end)

    cpdef tuple all_peaks_for(self, double neutral_mass, double tolerance=1e-5):
        """Find all peaks that match ``neutral_mass`` within ``error_tolerance`` mass accuracy ppm.

        Parameters
        ----------
        neutral_mass: double
            The mass to search for
        error_tolerance: double
            The PPM error tolerance to apply

        Returns
        -------
        tuple of DeconvolutedPeak
            The found peaks
        """
        cdef:
            double lo, hi, lo_err, hi_err
            int lo_ix, hi_ix
            DeconvolutedPeak lo_peak
            DeconvolutedPeak hi_peak
        if not self.indexed:
            self.reindex()
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
        return <tuple>PyTuple_GetSlice(self.peaks, lo_ix, hi_ix)

    cdef DeconvolutedPeak _get_nearest_peak(self, double neutral_mass, double* errout):
        """Find the peak nearest to ``neutral_mass``, regardless of error.

        Parameters
        ----------
        neutral_mass: double
            The mass to search for
        errout: double*
            A pointer to store the mass error in

        Returns
        -------
        DeconvolutedPeak
            The nearest peak
        """
        cdef DeconvolutedPeak peak
        peak = binary_search_nearest_neutral_mass(self.peaks, neutral_mass, errout)
        return peak

    def get_nearest_peak(self, double neutral_mass):
        """Find the peak nearest to ``neutral_mass``, regardless of error.

        Parameters
        ----------
        neutral_mass: double
            The mass to search for

        Returns
        -------
        DeconvolutedPeak
            The nearest peak
        double
            The error between ``neutral_mass`` and the found peak
        """
        cdef:
            DeconvolutedPeak peak
            double errout
        peak = self._get_nearest_peak(neutral_mass, &errout)
        return peak, errout

    def between(self, m1, m2, tolerance=1e-5, use_mz=False):
        """Retrieve a :class:`DeconvolutedPeakSet` containing all the peaks
        whose mass is between ``m1`` and ``m2``.

        If ``use_mz`` is :const:`True` then search by m/z instead of mass

        Parameters
        ----------
        m1 : float
            The lower mass limit
        m2 : float
            The upper mass limit
        use_mz: bool
            Whether to search for m/z instead of neutral mass

        Returns
        -------
        DeconvolutedPeakSet
        """
        acc = []
        collecting = False
        value = None
        if not use_mz:
            getter = operator.attrgetter("neutral_mass")
            iterator = self.peaks
        else:
            getter = operator.attrgetter("mz")
            iterator = self._mz_ordered
        for peak in iterator:
            value = getter(peak)
            if not collecting and value >= m1 and value < m2:
                collecting = True
            elif collecting and value > m2:
                break

            if collecting:
                acc.append(peak.clone())

        return self.__class__(acc).reindex()

    def select_top_n(self, n: int=200):
        top_peaks = sorted(self, key=operator.attrgetter('intensity'), reverse=True)[:n]
        dup = self.__class__((p.clone() for p in top_peaks))
        dup.reindex()
        return dup


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


@cython.cdivision(True)
cdef double _ppm_error(double x, double y) nogil:
    return (x - y) / y


cdef int _binary_search(double* array, double target, double error_tolerance, size_t n, size_t* out) nogil:
    cdef:
        size_t lo, hi, mid, i, j
        size_t best_index
        double err, found_mass
        double best_error, abs_error

    lo = 0
    hi = n

    while hi != lo:
        mid = (hi + lo) // 2
        found_mass = array[mid]
        if fabs(_ppm_error(found_mass, target)) < error_tolerance:
            best_index = mid
            best_error = INF
            i = 0
            while mid - i >= 0 and i <= mid:
                found_mass = array[mid - i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid - i
                        best_error = abs_error
                else:
                    break
                i += 1
            i = 1
            while (mid + i) < (n - 1):
                found_mass = array[mid + i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid + i
                        best_error = abs_error
                else:
                    break
                i += 1
            out[0] = best_index
            if best_error == INF:
                return 3
            return 0
        elif hi - lo == 1:
            out[0] = 0
            return 1
        elif found_mass > target:
            hi = mid
        else:
            lo = mid
    out[0] = 0
    return 2


cdef int _binary_search_interval(double* array, double target, double error_tolerance, size_t n, size_t* start, size_t* end) nogil:
    cdef:
        size_t lo, hi, mid, i, j
        size_t best_index
        double err, found_mass
        double best_error, abs_error

    lo = 0
    hi = n

    while hi != lo:
        mid = (hi + lo) // 2
        found_mass = array[mid]
        if fabs(_ppm_error(found_mass, target)) < error_tolerance:
            best_index = mid
            best_error = INF
            i = 0
            while mid - i >= 0 and i <= mid:
                found_mass = array[mid - i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid - i
                        best_error = abs_error
                else:
                    break
                i += 1
            start[0] = mid - i + 1
            i = 1
            while (mid + i) < (n - 1):
                found_mass = array[mid + i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid + i
                        best_error = abs_error
                else:
                    break
                i += 1
            end[0] = mid + i
            return 0
        elif hi - lo == 1:
            start[0] = 0
            end[0] = 0
            return 1
        elif found_mass > target:
            hi = mid
        else:
            lo = mid
    start[0] = 0
    end[0] = 0
    return 2


cdef DeconvolutedPeak binary_search_neutral_mass(tuple peak_set, double neutral_mass, double error_tolerance):
    cdef:
        size_t lo, hi, mid, i, j
        double err, found_mass
        DeconvolutedPeak found_peak

    lo = 0
    hi = PyTuple_GET_SIZE(peak_set)

    while hi != lo:
        mid = (hi + lo) // 2
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
        mid = (hi + lo) // 2
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
        mid = (hi + lo) // 2
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


@cython.boundscheck(False)
@cython.wraparound(False)
def decode_envelopes(np.ndarray[cython.floating, ndim=1, mode='c'] array):
    cdef:
        list envelope_list, current_envelope
        tuple converted
        Envelope envelope_obj
        size_t i, n
        cython.floating a, b
    envelope_list = []
    current_envelope = []

    i = 0
    n = array.shape[0]
    while i < n:
        a = array[i]
        b = array[i + 1]
        i += 2
        if a == 0 and b == 0:
            if current_envelope is not None:
                if current_envelope:
                    converted = tuple(current_envelope)
                    envelope_obj = Envelope._create(converted)
                    envelope_list.append(envelope_obj)
                current_envelope = []
        else:
            current_envelope.append(EnvelopePair._create(a, b))
    converted = tuple(current_envelope)
    envelope_obj = Envelope._create(converted)
    envelope_list.append(envelope_obj)
    return envelope_list


INTERVAL_INDEX_SIZE = 0

cdef class DeconvolutedPeakSetIndexed(DeconvolutedPeakSet):
    def __init__(self, peaks):
        self.neutral_mass_array = NULL
        self.mz_array = NULL
        self.interval_index = NULL
        super(DeconvolutedPeakSetIndexed, self).__init__(peaks)
        self._size = self.get_size()

    def __dealloc__(self):
        self._release_buffers()

    def __reduce__(self):
        return self.__class__, (self.peaks, ), self.__getstate__()

    def __getstate__(self):
        d = {"indexed": self.indexed}
        return d

    def __setstate__(self, d):
        if d.get("indexed", False):
            self.reindex()

    def _release_buffers(self):
        if self.neutral_mass_array != NULL:
            free(self.neutral_mass_array)
            self.neutral_mass_array = NULL
        if self.mz_array != NULL:
            free(self.mz_array)
            self.mz_array = NULL
        if self.interval_index != NULL:
            free_index_list(self.interval_index)
            self.interval_index = NULL

    def set_interval_index_size(self, index_size):
        global INTERVAL_INDEX_SIZE
        INTERVAL_INDEX_SIZE = index_size
        self.reindex()

    cdef void _build_index_arrays(self):
        cdef:
            size_t i, n
            DeconvolutedPeak peak
        n = PyTuple_GET_SIZE(self.peaks)
        self._size = n
        self._release_buffers()
        self.neutral_mass_array = <double*>malloc(sizeof(double) * n)
        self.mz_array = <double*>malloc(sizeof(double) * n)

        for i in range(n):
            peak = self.getitem(i)
            self.neutral_mass_array[i] = peak.neutral_mass
            self.mz_array[peak._index.mz] = peak.mz

        # The interpolating interval index is just a bit slower than
        # the normal binary search, likely due to cache friendliness.
        if n > 2 and INTERVAL_INDEX_SIZE > 0:
            if self.interval_index != NULL:
                free_index_list(self.interval_index)
                self.interval_index = NULL
            interval_index = <index_list*>malloc(sizeof(index_list))
            build_interval_index(self, interval_index, INTERVAL_INDEX_SIZE)
            if check_index(interval_index) != 0:
                free_index_list(interval_index)
            else:
                self.interval_index = interval_index

    cpdef reindex(self):
        super(DeconvolutedPeakSetIndexed, self).reindex()
        self._build_index_arrays()
        return self

    cdef DeconvolutedPeak _has_peak(self, double neutral_mass, double error_tolerance=1e-5, bint use_mz=False):
        """Find the peak that best matches ``neutral_mass`` within ``error_tolerance`` mass accuracy ppm.

        If ``use_mz`` is True, instead of matching neutral masses, match peaks using m/z instead.

        Parameters
        ----------
        neutral_mass: double
            The mass to search for
        error_tolerance: double
            The PPM error tolerance to apply
        use_mz: bool
            Whether to search using m/z instead of neutral mass

        Returns
        -------
        DeconvolutedPeak
            The found peak, or None if no peak is found
        """
        cdef:
            int status
            size_t i, n, s
            DeconvolutedPeak peak

        n = self._size
        if use_mz:
            status = _binary_search(self.mz_array, neutral_mass, error_tolerance, n, &i)
            if status != 0:
                return None
            else:
                return <DeconvolutedPeak>PyTuple_GET_ITEM(self._mz_ordered, i)
        if n == 0:
            return None
        if self.interval_index != NULL:
            find_search_interval(self.interval_index, neutral_mass, &s, &n)
            status = _binary_search_with_hint(self.neutral_mass_array, neutral_mass, error_tolerance, n, s, &i)
        else:
            status = _binary_search(self.neutral_mass_array, neutral_mass, error_tolerance, n, &i)
        if status != 0:
            return None
        peak = self.getitem(i)
        if abs((peak.neutral_mass - neutral_mass) / neutral_mass) < error_tolerance:
            return peak
        else:
            return None

    def _test_interval(self, double neutral_mass, double tolerance=1e-5):
        cdef:
            int status
            size_t n, start, end
        n = self._size
        status = _binary_search_interval(
            self.neutral_mass_array, neutral_mass, tolerance, n, &start, &end)
        return status, start, end

    cpdef tuple all_peaks_for(self, double neutral_mass, double tolerance=1e-5):
        """Find all peaks that match ``neutral_mass`` within ``error_tolerance`` mass accuracy ppm.

        Parameters
        ----------
        neutral_mass: double
            The mass to search for
        error_tolerance: double
            The PPM error tolerance to apply

        Returns
        -------
        tuple of DeconvolutedPeak
            The found peaks
        """
        cdef:
            int status
            size_t n, s, start, end
        if not self.indexed:
            self.reindex()
        n = self._size
        if n == 0:
            return ()
        if self.interval_index != NULL:
            find_search_interval(self.interval_index, neutral_mass, &s, &n)
            status = _binary_search_interval_with_hint(
                self.neutral_mass_array, neutral_mass, tolerance, n, s, &start, &end)
        else:
            status = _binary_search_interval(
                self.neutral_mass_array, neutral_mass, tolerance, n, &start, &end)
        if status != 0:
            return ()
        return <tuple>PyTuple_GetSlice(self.peaks, start, end)

    cdef int _interval_for(self, double neutral_mass, double tolerance, size_t* start, size_t* end) nogil:
        cdef:
            int status
            size_t n
        n = self._size
        status = _binary_search_interval(
            self.neutral_mass_array, neutral_mass, tolerance, n, start, end)
        return status

    def test_interval(self, double value):
        cdef:
            int status
            size_t start, end
            size_t i

        status = find_search_interval(self.interval_index, value, &start, &end)
        return start, end, interpolate_index(self.interval_index, value)

    def find_interval_for(self, double value):
        return interpolate_index(self.interval_index, value)

    def check_interval(self, size_t i):
        cdef:
            index_cell cell
        cell = self.interval_index.index[i]
        return cell

    cpdef DeconvolutedPeakSet copy(self):
        cdef DeconvolutedPeakSetIndexed dup = DeconvolutedPeakSet.copy(self)
        dup._build_index_arrays()
        return dup



cdef int _binary_search_with_hint(double* array, double target, double error_tolerance, size_t n, size_t hint, size_t* out) nogil:
    """Performs a best-matching binary search with a hint at the lower bound as well as the upper bound.

    The best-matching step occurs once a match is found, iterating forwards and backwards from the initial
    match until an index is found that minimizes the parts-per-million error with ``target``

    The binary search starts with the lower bound ``lo = hint`` instead of ``lo = 0``. The upper bound
    ``hi = n`` can be used as a hint to instruct the search to only consider a small segment of the array.

    Parameters
    ----------
    array: double*
        The array of doubles to search within
    target: double
        The value to search for
    error_tolerance: double
        The parts-per-million error tolerance to use to select solutions
    n: size_t
        The upper limit into ``array`` to search, usually the size of the array
        but may also be made lower to reduce the space to search
    hint: size_t
        The lower limit into ``array`` to search. Leaving this as 0
        will make this behave the same as :func:`_binary_search`.
    out: size_t*
        The best matching index
    """
    cdef:
        size_t lo, hi, mid, i, j
        size_t best_index
        double err, found_mass
        double best_error, abs_error

    lo = hint
    hi = n

    while hi != lo:
        mid = (hi + lo) // 2
        found_mass = array[mid]
        if fabs(_ppm_error(found_mass, target)) < error_tolerance:
            best_index = mid
            best_error = INF
            i = 0
            while mid - i >= 0 and i <= mid:
                found_mass = array[mid - i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid - i
                        best_error = abs_error
                else:
                    break
                i += 1
            i = 1
            while (mid + i) < (n - 1):
                found_mass = array[mid + i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid + i
                        best_error = abs_error
                else:
                    break
                i += 1
            out[0] = best_index
            if best_error == INF:
                return 3
            return 0
        elif hi - lo == 1:
            out[0] = 0
            return 1
        elif found_mass > target:
            hi = mid
        else:
            lo = mid
    out[0] = 0
    return 2


cdef int _binary_search_interval_with_hint(double* array, double target, double error_tolerance, size_t n, size_t hint, size_t* start, size_t* end) nogil:
    cdef:
        size_t lo, hi, mid, i, j
        size_t best_index
        double err, found_mass
        double best_error, abs_error

    lo = hint
    hi = n

    while hi != lo:
        mid = (hi + lo) // 2
        found_mass = array[mid]
        if fabs(_ppm_error(found_mass, target)) < error_tolerance:
            best_index = mid
            best_error = INF
            i = 0
            while mid - i >= 0 and i <= mid:
                found_mass = array[mid - i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid - i
                        best_error = abs_error
                else:
                    break
                i += 1
            start[0] = mid - i + 1
            i = 1
            while (mid + i) < (n - 1):
                found_mass = array[mid + i]
                abs_error = fabs(_ppm_error(target, found_mass))
                if abs_error < error_tolerance:
                    if abs_error < best_error:
                        best_index = mid + i
                        best_error = abs_error
                else:
                    break
                i += 1
            end[0] = mid + i
            return 0
        elif hi - lo == 1:
            start[0] = 0
            end[0] = 0
            return 1
        elif found_mass > target:
            hi = mid
        else:
            lo = mid
    start[0] = 0
    end[0] = 0
    return 2


cdef int build_interval_index(DeconvolutedPeakSet peaks, index_list* index, size_t index_size):
    cdef:
        double* linear_spacing
        double current_value, err, next_value
        size_t i, start_i, end_i, index_i, peaks_size
        DeconvolutedPeak peak
    peaks_size = peaks.get_size()
    if peaks_size > 0:
        index.low = peaks.getitem(0).neutral_mass
        index.high = peaks.getitem(peaks_size - 1).neutral_mass
    else:
        index.low = 0
        index.high = 1
    index.size = index_size
    linear_spacing = build_linear_spaced_array(
        index.low,
        index.high,
        index_size)

    index.index = <index_cell*>malloc(sizeof(index_cell) * index_size)

    for index_i in range(index_size):
        current_value = linear_spacing[index_i]
        peak = peaks._get_nearest_peak(current_value, &err)
        if peaks_size > 0:
            start_i = peak._index.neutral_mass
            if index_i > 0:
                start_i = index.index[index_i - 1].end - 1
            if index_i == index_size - 1:
                end_i = peaks_size - 1
            else:
                next_value = linear_spacing[index_i + 1]
                i = peak._index.neutral_mass
                while i < peaks_size:
                    peak = peaks.getitem(i)
                    if abs(current_value - peak.neutral_mass) > abs(next_value - peak.neutral_mass):
                        break
                    i += 1
                end_i = i
        else:
            start_i = 0
            end_i = 0
        index.index[index_i].center_value = current_value
        index.index[index_i].start = start_i
        index.index[index_i].end = end_i

    free(linear_spacing)
    return 0


cdef double* build_linear_spaced_array(double low, double high, size_t n):
    cdef:
        double* array
        double delta
        size_t i

    delta = (high - low) / (n - 1)

    array = <double*>malloc(sizeof(double) * n)

    for i in range(n):
        array[i] = low + i * delta
    return array


cdef void free_index_list(index_list* index):
    free(index.index)
    free(index)


cdef int check_index(index_list* index) nogil:
    if index.size == 0:
        return 1
    elif (index.high - index.low) == 0:
        return 2
    else:
        return 0


cdef size_t interpolate_index(index_list* index, double value):
    cdef:
        double v
        size_t i
    v = (((value - index.low) / (index.high - index.low)) * index.size)
    i = <size_t>v
    return i


cdef int find_search_interval(index_list* index, double value, size_t* start, size_t* end):
    cdef:
        size_t i
    if value > index.high:
        i = index.size - 1
    elif value < index.low:
        i = 0
    else:
        i = interpolate_index(index, value)
    if i > 0:
        if i < index.size:
            start[0] = index.index[i - 1].start
        else:
            # if we're at index.size or greater, act as if we're at the last
            # cell of the index
            start[0] = index.index[index.size - 2].start
    else:
        # if somehow the index were negative, this could end badly.
        start[0] = index.index[i].start
    if i >= (index.size - 1):
        end[0] = index.index[index.size - 1].end + 1
    else:
        end[0] = index.index[i + 1].end + 1
    return 0


cdef int create_deconvoluted_peak_set_t(DeconvolutedPeakSetIndexed peak_set, deconvoluted_peak_set_t* destination) except 1:
    cdef:
        size_t size, i
        deconvoluted_peak_t* peaks
        double* index
        DeconvolutedPeak peak


    size = peak_set.get_size()
    peaks = <deconvoluted_peak_t*>malloc(sizeof(deconvoluted_peak_t) * size)
    if peaks == NULL:
        PyErr_SetString(MemoryError, "Failed to allocate peak array for C peak structure")
        return 1

    index = <double*>malloc(sizeof(double) * size)
    if index == NULL:
        free(peaks)
        PyErr_SetString(MemoryError, "Failed to allocate mass index for C peak structure")
        return 1

    for i in range(size):
        peak = peak_set.getitem(i)
        peaks[i].neutral_mass = peak.neutral_mass
        peaks[i].charge = peak.charge
        peaks[i].intensity = peak.intensity
        peaks[i].index = i

        index[i] = peak.neutral_mass

    destination.peaks = peaks
    destination.mass_index = index
    destination.size = size
    destination.flags = CPeaksFlags.owns_peaks_and_index
    return 0


cdef int free_deconvoluted_peak_set_t(deconvoluted_peak_set_t* destination) nogil:
    if destination.flags & CPeaksFlags.owns_peaks:
        free(destination.peaks)
    if destination.flags & CPeaksFlags.owns_index:
        free(destination.mass_index)
    return 0


cdef deconvoluted_peak_set_t deconvoluted_peak_set_all_peaks_for(deconvoluted_peak_set_t* self, double neutral_mass, double error_tolerance=2e-5) nogil:
    cdef:
        size_t n, s, i
        size_t start, end
        int status
        deconvoluted_peak_set_t result

    n = self.size
    start = 0
    end = 0

    result.peaks = NULL
    result.mass_index = NULL
    result.size = 0
    result.flags = CPeaksFlags.borrowing

    status = _binary_search_interval(
        self.mass_index, neutral_mass, error_tolerance, n, &start, &end)

    if status != 0:
        return result

    result.peaks = &self.peaks[start]
    result.mass_index = &self.mass_index[start]
    result.size = end - start
    return result


cdef deconvoluted_peak_t* deconvoluted_peak_set_has_peak(deconvoluted_peak_set_t* self, double neutral_mass, double error_tolerance=2e-5) nogil:
    cdef:
        size_t n, s, i
        deconvoluted_peak_t* peak
    n = self.size
    s = 0
    i = 0
    status = _binary_search_with_hint(self.mass_index, neutral_mass, error_tolerance, n, s, &i)
    if i < n and i > 0:
        peak = &self.peaks[i]
        if fabs((peak.neutral_mass - neutral_mass) / neutral_mass) < error_tolerance:
            return peak
    return NULL


cdef int deconvoluted_peak_eq(deconvoluted_peak_t* self, deconvoluted_peak_t* other) nogil:
    if fabs(self.neutral_mass - other.neutral_mass) > 1e-6:
        return False
    if self.charge != other.charge:
        return False
    if fabs(self.intensity - other.intensity) > 1e-3:
        return False
    return True


cdef size_t deconvoluted_peak_hash(deconvoluted_peak_t* self) nogil:
    cdef:
        size_t value

    value = 0xBADFADA << self.charge
    value &= <long>self.neutral_mass
    return value


@cython.final
cdef class _CPeakSet:

    @classmethod
    def from_peak_list(cls, DeconvolutedPeakSetIndexed peaks):
        cdef:
            deconvoluted_peak_set_t* ptr
            _CPeakSet self

        ptr = <deconvoluted_peak_set_t*>malloc(sizeof(deconvoluted_peak_set_t))
        if ptr == NULL:
            raise MemoryError()
        if create_deconvoluted_peak_set_t(peaks, ptr) != 0:
            raise ValueError("Failed to convert peak set")
        self = cls.__new__(cls)
        self.ptr = ptr
        return self

    def __dealloc__(self):
        if self.ptr != NULL:
            free_deconvoluted_peak_set_t(self.ptr)
            self.ptr = NULL

    def __getitem__(self, int i):
        if i >= self.ptr.size or i < 0:
            raise IndexError(i)
        return self.ptr.peaks[i]

    def __len__(self):
        return self.ptr.size

    cdef deconvoluted_peak_t* getitem(self, size_t i) nogil:
        if i >= self.ptr.size or i < 0:
            return NULL
        return &self.ptr.peaks[i]

    cdef deconvoluted_peak_set_t _all_peaks_for(self, double neutral_mass, double error_tolerance=2e-5) nogil:
        cdef:
            deconvoluted_peak_set_t result

        result = deconvoluted_peak_set_all_peaks_for(self.ptr, neutral_mass, error_tolerance)
        return result

    cpdef tuple all_peaks_for(self, double neutral_mass, double error_tolerance=2e-5):
        cdef:
            size_t n, i
            tuple out
            deconvoluted_peak_set_t result
            # deconvoluted_peak_t peak
            object peak

        result = self._all_peaks_for(neutral_mass, error_tolerance)
        out = PyTuple_New(result.size)
        for i in range(result.size):
            peak = result.peaks[i]
            Py_INCREF(peak)
            PyTuple_SetItem(out, i, peak)
        return out

    @cython.cdivision(True)
    cdef deconvoluted_peak_t* _has_peak(self, double neutral_mass, double error_tolerance=2e-5) nogil:
        cdef:
            size_t n, s, i
            deconvoluted_peak_t* peak

        return deconvoluted_peak_set_has_peak(self.ptr, neutral_mass, error_tolerance)

    cpdef has_peak(self, double neutral_mass, double error_tolerance=2e-5):
        cdef:
            size_t n, s, i
            deconvoluted_peak_t* peak
        peak = self._has_peak(neutral_mass, error_tolerance)
        if peak == NULL:
            return None
        else:
            return peak[0]
