import operator
from .utils import Base, ppm_error
from brainpy import mass_charge_ratio


class DeconvolutedPeak(Base):
    __slots__ = [
        "neutral_mass", "intensity", "signal_to_noise",
        "index", "full_width_at_half_max", "charge",
        "a_to_a2_ratio", "most_abundant_mass", "average_mass",
        "score", "envelope"
    ]

    def __init__(self, neutral_mass, intensity, charge, signal_to_noise, index, full_width_at_half_max,
                 a_to_a2_ratio=None, most_abundant_mass=None, average_mass=None, score=None,
                 envelope=None):
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


class DeconvolutedPeakSet(Base):
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
            return _get_nearest_peak(self._mz_ordered, neutral_mass, use_mz=use_mz)
        else:
            return _get_nearest_peak(self.peaks, neutral_mass, use_mz=use_mz)

    def has_peak(self, neutral_mass, tolerance=1e-5, use_mz=False):
        if use_mz:
            return binary_search(self._mz_ordered, neutral_mass, tolerance, mz_getter)
        return binary_search(self.peaks, neutral_mass, tolerance, neutral_mass_getter)

    def contains_mz_in_envelope(self, mz, isolation_window=2., tolerance=1e-5):
        try:
            return max(
                [peak for peak in self.between(
                    mz - isolation_window,
                    mz + isolation_window,
                    use_mz=True)],
                key=operator.attrgetter("score"))
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


mz_getter = operator.attrgetter('mz')
neutral_mass_getter = operator.attrgetter("neutral_mass")


def _get_nearest_peak(peaklist, neutral_mass, use_mz=False):
    lo = 0
    hi = len(peaklist)

    tol = 1

    getter = operator.attrgetter("mz") if use_mz else operator.attrgetter("neutral_mass")

    def sweep(lo, hi):
        best_error = float('inf')
        best_index = None
        for i in range(hi - lo):
            i += lo
            v = getter(peaklist[i])
            err = abs(v - neutral_mass)
            if err < best_error:
                best_error = err
                best_index = i
        return peaklist[best_index], best_error

    def binsearch(lo, hi):
        if (hi - lo) < 5:
            return sweep(lo, hi)
        else:
            mid = (hi + lo) / 2
            v = getter(peaklist[mid])
            if abs(v - neutral_mass) < tol:
                return sweep(lo, hi)
            elif v > neutral_mass:
                return binsearch(lo, mid)
            else:
                return binsearch(mid, hi)
    return binsearch(lo, hi)


def _sweep_solution(array, value, lo, hi, tolerance, getter=neutral_mass_getter):
    best_index = -1
    best_error = float('inf')
    for i in range(hi - lo):
        target = array[lo + i]
        error = ppm_error(value, getter(target))
        abs_error = abs(error)
        if abs_error < tolerance and abs_error < best_error:
            best_index = lo + i
            best_error = abs_error
    if best_index == -1:
        return None
    else:
        return array[best_index]


def _binary_search(array, value, lo, hi, tolerance, getter=neutral_mass_getter):
    if (hi - lo) < 5:
        return _sweep_solution(array, value, lo, hi, tolerance, getter)
    else:
        mid = (hi + lo) / 2
        target = array[mid]
        target_value = getter(target)
        error = ppm_error(value, target_value)

        if abs(error) <= tolerance:
            return _sweep_solution(array, value, max(mid - 5, lo), min(mid + 5, hi), tolerance, getter)
        elif target_value > value:
            return _binary_search(array, value, lo, mid, tolerance, getter)
        elif target_value < value:
            return _binary_search(array, value, mid, hi, tolerance, getter)
    raise Exception("No recursion found!")


def binary_search(array, value, tolerance=2e-5, getter=neutral_mass_getter):
    size = len(array)
    if size == 0:
        return None
    return _binary_search(array, value, 0, size, tolerance, getter)


def _binary_search_envelope(array, value, lo, hi, tolerance):
    mid = (hi + lo) / 2
    target = array[mid]
    if target.contains_mz(value, tolerance):
        matches = [target]
        for i in range(max(lo - 100, 0), min(hi + 100, len(array))):
            target = array[i]
            if target.contains_mz(value, tolerance):
                matches.append(target)
        return max(matches, key=operator.attrgetter("intensity"))
    elif hi - lo == 1:
        return None
    elif target.mz < value:
        return _binary_search_envelope(array, value, mid, hi, tolerance)
    elif target.mz > value:
        return _binary_search_envelope(array, value, lo, mid, tolerance)
