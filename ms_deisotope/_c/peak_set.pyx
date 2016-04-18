# cython: embedsignature=True

import operator
from cpython.tuple cimport PyTuple_GET_ITEM, PyTuple_GetItem, PyTuple_GetSlice
from cpython cimport PyObject
from brainpy import calculate_mass, PROTON as _PROTON

cdef double PROTON
PROTON = _PROTON

mz_getter = operator.attrgetter('mz')
neutral_mass_getter = operator.attrgetter("neutral_mass")


cdef double ppm_error(double x, double y):
    return (x - y) / y


cdef double mass_charge_ratio(double neutral_mass, int charge, double charge_carrier=PROTON):
    return (neutral_mass + charge * charge_carrier) / abs(charge)


cdef class DeconvolutedPeak(object):
    cdef:
        public double neutral_mass
        public double intensity
        public int charge
        public double signal_to_noise
        public long index
        public double full_width_at_half_max
        public double a_to_a2_ratio
        public double most_abundant_mass
        public double average_mass
        public double score
        public list envelope
        public double _mz

    def __cinit__(self, neutral_mass, intensity, charge, signal_to_noise, index, full_width_at_half_max,
                 a_to_a2_ratio, most_abundant_mass, average_mass, score, envelope):
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
        self.envelope = envelope
        self._mz = self.mz

    def clone(self):
        return DeconvolutedPeak(self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
                                self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
                                self.most_abundant_mass, self.average_mass, self.score,
                                self.envelope)

    def __reduce__(self):
        return DeconvolutedPeak, (self.neutral_mass, self.intensity, self.charge, self.signal_to_noise,
                                  self.index, self.full_width_at_half_max, self.a_to_a2_ratio,
                                  self.most_abundant_mass, self.average_mass, self.score,
                                  self.envelope)

    @property
    def mz(self):
        return mass_charge_ratio(self.neutral_mass, self.charge)

    def contains_mz(self, mz, tolerance=2e-5):
        for i, m in enumerate(self.envelope):
            if abs(m[0] - mz) / mz <= tolerance:
                return True
        return False

    def __repr__(self):
        fields = [
            "neutral_mass", "intensity", "signal_to_noise",
            "index", "full_width_at_half_max", "charge",
            "a_to_a2_ratio", "most_abundant_mass", "average_mass",
            "score", "envelope"
        ]

        return "DeconvolutedPeak(%s)" % ', '.join("%s=%r" % (attr, getattr(self, attr)) for attr in fields)


cdef class DeconvolutedPeakSet(object):
    cdef:
        public tuple peaks
        public tuple _mz_ordered

    def __init__(self, peaks):
        self.peaks = peaks
        self._mz_ordered = None

    def _reindex(self):
        self.peaks = tuple(sorted(self.peaks, key=operator.attrgetter("neutral_mass")))
        self._mz_ordered = tuple(sorted(self.peaks, key=operator.attrgetter("mz")))
        for i, peak in enumerate(self.peaks):
            peak.index = i
        return self

    def __len__(self):
        return len(self.peaks)

    def get_nearest_peak(self, neutral_mass, use_mz=False):
        if use_mz:
            return _binary_search_nearest_match_mz(self._mz_ordered, neutral_mass, 0, len(self.peaks), 2e-4)
        else:
            return _binary_search_nearest_match_neutral_mass(self.peaks, neutral_mass, 0, len(self.peaks), 2e-4)

    def has_peak(self, neutral_mass, tolerance=1e-5, use_mz=False):
        if use_mz:
            return _binary_search_ppm_error_mz(self._mz_ordered, neutral_mass, 0, len(self.peaks), tolerance)
        return _binary_search_ppm_error_neutral_mass(self.peaks, neutral_mass, 0, len(self.peaks), tolerance)

    def contains_mz_in_envelope(self, mz, isolation_window=1.5, tolerance=1e-5):
        try:
            return max(
                [peak for peak in self.between(
                    mz - isolation_window,
                    mz + isolation_window,
                    use_mz=True)],
                key=operator.attrgetter("intensity"))
        except ValueError:
            return None
        # return _binary_search_envelope(self._mz_ordered, mz, 0, len(self), tolerance)

    def __repr__(self):
        return "<DeconvolutedPeakSet %d Peaks>" % (len(self))

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.__class__(tuple(p.clone() for p in self.peaks[item]))
        return self.peaks[item]

    def clone(self):
        return self.__class__(tuple(p.clone() for p in self))

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
            elif collecting:
                acc.append(peak)
        return self.__class__(acc)


cdef DeconvolutedPeak _null_peak
_null_peak = DeconvolutedPeak(0,0,0,0,0,0,0,0,0,0,[])


cdef DeconvolutedPeak _binary_search_ppm_error_neutral_mass(tuple array, double value, size_t lo, size_t hi, double tolerance):
    cdef:
        size_t mid, lower_edge
        DeconvolutedPeak target
        double target_value, error
    if (hi - lo) < 5:
        return _sweep_solution_ppm_error_neutral_mass(array, value, lo, hi, tolerance)
    else:
        mid = (hi + lo) / 2
        target = <DeconvolutedPeak>PyTuple_GetItem(array, mid)
        target_value = target.neutral_mass
        error = ppm_error(value, target_value)
        if abs(error) <= tolerance:
            return _sweep_solution_ppm_error_neutral_mass(array, value, max(mid - (mid if mid < 5 else 5), lo), min(mid + 5, hi), tolerance)
        elif target_value > value:
            return _binary_search_ppm_error_neutral_mass(array, value, lo, mid, tolerance)
        elif target_value < value:
            return _binary_search_ppm_error_neutral_mass(array, value, mid, hi, tolerance)
    return _null_peak


cdef DeconvolutedPeak _binary_search_ppm_error_mz(tuple array, double value, size_t lo, size_t hi, double tolerance):
    cdef:
        size_t mid, lower_edge
        DeconvolutedPeak target
        double target_value, error
    if (hi - lo) < 5:
        return _sweep_solution_ppm_error_mz(array, value, lo, hi, tolerance)
    else:
        mid = (hi + lo) / 2
        target = <DeconvolutedPeak>PyTuple_GetItem(array, mid)
        target_value = target._mz
        error = ppm_error(value, target_value)
        if abs(error) <= tolerance:
            return _sweep_solution_ppm_error_mz(array, value, max(mid - (mid if mid < 5 else 5), lo), min(mid + 5, hi), tolerance)
        elif target_value > value:
            return _binary_search_ppm_error_mz(array, value, lo, mid, tolerance)
        elif target_value < value:
            return _binary_search_ppm_error_mz(array, value, mid, hi, tolerance)
    return _null_peak


cdef DeconvolutedPeak _sweep_solution_ppm_error_neutral_mass(tuple array, double value, size_t lo, size_t hi, double tolerance):
    cdef:
        double best_error, error, abs_error
        size_t i, limit

    best_index = -1
    best_intensity = 0
    best_error = 1000000000000000
    for i in range(hi - lo):
        target = <DeconvolutedPeak>PyTuple_GetItem(array, lo + i)
        error = ppm_error(value, target.neutral_mass)
        abs_error = abs(error)
        if abs_error < tolerance and (abs_error < (best_error * 1.1)) and (target.intensity > best_intensity):
            best_index = lo + i
            best_error = abs_error
    if best_index == -1:
        return _null_peak
    else:
        return <DeconvolutedPeak>PyTuple_GetItem(array, best_index)


cdef DeconvolutedPeak _sweep_solution_ppm_error_mz(tuple array, double value, size_t lo, size_t hi, double tolerance):
    cdef:
        double best_error, error, abs_error
        size_t i, limit

    best_index = -1
    best_intensity = 0
    best_error = 1000000000000000
    for i in range(hi - lo):
        target = <DeconvolutedPeak>PyTuple_GetItem(array, lo + i)
        error = ppm_error(value, target._mz)
        abs_error = abs(error)
        if abs_error < tolerance and (abs_error < (best_error * 1.1)) and (target.intensity > best_intensity):
            best_index = lo + i
            best_error = abs_error
    if best_index == -1:
        return _null_peak
    else:
        return <DeconvolutedPeak>PyTuple_GetItem(array, best_index)


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


cdef DeconvolutedPeak _binary_search_nearest_match_neutral_mass(tuple array, double value, size_t lo, size_t hi, double* errout):
    cdef:
        size_t mid
        double v
        double err

    if (hi - lo) < 5:
        return _sweep_nearest_match_neutral_mass(array, value, lo, hi, errout)
    else:
        mid = (hi + lo) / 2
        v = array[mid].neutral_mass
        if abs(v - value) < 1.:
            return _sweep_nearest_match_neutral_mass(array, value, lo, hi, errout)
        elif v > value:
            return _binary_search_nearest_match_neutral_mass(array, value, lo, mid, errout)
        else:
            return _binary_search_nearest_match_neutral_mass(array, value, mid, hi, errout)


cdef DeconvolutedPeak _sweep_nearest_match_mz(tuple array, double value, size_t lo, size_t hi, double* errout):
    cdef:
        size_t i
        size_t best_index
        double best_error, err
        double v

    best_error = float('inf')
    best_index = -1
    for i in range(hi - lo):
        i += lo
        v = array[i]._mz
        err = abs(v - value)
        if err < best_error:
            best_error = err
            best_index = i
    errout[0] = best_error
    return array[best_index]


cdef DeconvolutedPeak _binary_search_nearest_match_mz(tuple array, double value, size_t lo, size_t hi, double* errout):
    cdef:
        size_t mid
        double v
        double err

    if (hi - lo) < 5:
        return _sweep_nearest_match_mz(array, value, lo, hi, errout)
    else:
        mid = (hi + lo) / 2
        v = array[mid]._mz
        if abs(v - value) < 1.:
            return _sweep_nearest_match_mz(array, value, lo, hi, errout)
        elif v > value:
            return _binary_search_nearest_match_mz(array, value, lo, mid, errout)
        else:
            return _binary_search_nearest_match_mz(array, value, mid, hi, errout)
