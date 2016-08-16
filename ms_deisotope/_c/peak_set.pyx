import operator

cimport cython
from ms_deisotope._c.averagine cimport mass_charge_ratio


@cython.cdivision
cdef double ppm_error(double x, double y):
    return (x - y) / y


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


cdef class EnvelopePair:
    pass


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
        return "[%s]" % (', '.join("(%0.4f, %0.2f)" % t for t in self),)

    def clone(self):
        return self.__class__(self)


cdef class DeconvolutedPeak:
    def __init__(self, neutral_mass, intensity, charge, signal_to_noise, index, full_width_at_half_max,
                 a_to_a2_ratio=None, most_abundant_mass=None, average_mass=None, score=None,
                 envelope=None, mz=None, fit=None, chosen_for_msms=False):
        if index is None:
            index = _Index()
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

    def __hash__(self):
        return hash((self.mz, self.intensity, self.charge))

    def clone(self):
        return DeconvolutedPeak(self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
                                self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
                                self.most_abundant_mass, self.average_mass, self.score,
                                self.envelope, self.mz, self.fit, self.chosen_for_msms)

    def __reduce__(self):
        return DeconvolutedPeak, (self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
                                  self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
                                  self.most_abundant_mass, self.average_mass, self.score,
                                  self.envelope, self.mz, self.fit, self.chosen_for_msms)

    cpdef bint _eq(self, DeconvolutedPeak other):
        return (abs(self.neutral_mass - other.neutral_mass) < 1e-5) and (
            abs(self.intensity - other.intensity) < 1e-5)

    def __richcmp__(self, other, int code):
        if code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))


cdef class DeconvolutedPeakSolution(DeconvolutedPeak):
    def __init__(self, solution, fit, *args, **kwargs):
        self.solution = solution
        self.fit = fit
        super(DeconvolutedPeakSolution, self).__init__(*args, **kwargs)

    def clone(self):
        return DeconvolutedPeakSolution(
            self.solution, self.fit, self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.chosen_for_msms)

    def __reduce__(self):
        return DeconvolutedPeakSolution, (
            self.solution, self.fit, self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
            self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
            self.most_abundant_mass, self.average_mass, self.score,
            self.envelope, self.mz, self.chosen_for_msms)

    def __iter__(self):
        yield self.solution
        yield self
        yield self.fit


cdef class DeconvolutedPeakSet:

    def __init__(self, peaks):
        self.peaks = peaks
        self._mz_ordered = None

    def _reindex(self):
        """
        Updates the :attr:`index` of each peak in `self` and updates the
        sorted order.

        Returns
        -------
        self: DeconvolutedPeakSet
        """
        self.peaks = tuple(sorted(self.peaks, key=operator.attrgetter("neutral_mass")))
        self._mz_ordered = tuple(sorted(self.peaks, key=operator.attrgetter("mz")))
        for i, peak in enumerate(self.peaks):
            peak.index = _Index()
            peak.index.neutral_mass = i
        for i, peak in enumerate(self._mz_ordered):
            peak.index.mz = i
        return self

    def __len__(self):
        return len(self.peaks)

    cdef DeconvolutedPeak _has_peak(self, double neutral_mass, double error_tolerance=2e-5, bint use_mz=False):
        if use_mz:
            return None
        else:
            return binary_search_neutral_mass(self.peaks, neutral_mass, error_tolerance)

    cpdef DeconvolutedPeak has_peak(self, double neutral_mass, double error_tolerance=2e-5, bint use_mz=False):
        return self._has_peak(neutral_mass, error_tolerance, use_mz)



cdef double INF
INF = float('inf')


cdef DeconvolutedPeak _sweep_solution_neutral_mass(tuple array, double value, size_t mid, double tolerance):
    cdef:
        size_t best_index, i, n
        double best_error, abs_error
        DeconvolutedPeak target

    best_index = mid
    best_error = INF
    n = len(array)

    i = 0
    while mid - i != 0:
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
    while (mid + i) != (n - 1):
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
    hi = len(peak_set)

    while hi != lo:
        mid = (hi + lo) / 2
        found_peak = <DeconvolutedPeak>peak_set[mid]
        found_mass = found_peak.neutral_mass

        if abs(ppm_error(found_mass, neutral_mass)) < error_tolerance:
            _sweep_solution_neutral_mass(peak_set, neutral_mass, mid, error_tolerance)
            return found_peak
        elif hi - lo == 1:
            return None
        elif found_mass > neutral_mass:
            hi = mid
        else:
            lo = mid
