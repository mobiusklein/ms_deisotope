import operator
from .utils import Base, ppm_error
from brainpy import mass_charge_ratio


class DeconvolutedPeak(Base):
    def __init__(self, neutral_mass, intensity, charge, signal_to_noise, index, full_width_at_half_max, score=None):
        self.neutral_mass = neutral_mass
        self.intensity = intensity
        self.signal_to_noise = signal_to_noise
        self.index = index
        self.full_width_at_half_max = full_width_at_half_max
        self.charge = charge
        self.score = score

    def clone(self):
        return DeconvolutedPeak(self.mz, self.intensity, self.charge, self.signal_to_noise,
                                self.index, self.full_width_at_half_max, self.score)

    @property
    def mz(self):
        return mass_charge_ratio(self.neutral_mass, self.charge)


class DeconvolutedPeakSet(Base):
    def __init__(self, peaks):
        self.peaks = tuple(sorted(peaks, key=operator.attrgetter("neutral_mass")))
        for i, peak in enumerate(peaks):
            peak.index = i

    def __len__(self):
        return len(self.peaks)

    def has_peak(self, neutral_mass, tolerance=1e-5):
        return binary_search(self.peaks, neutral_mass, tolerance)

    def __repr__(self):
        return "<DeconvolutedPeakSet %d Peaks>" % (len(self))

    def __getitem__(self, item):
        if isinstance(item, slice):
            return DeconvolutedPeakSet(self.peaks[item])
        return self.peaks[item]

    def clone(self):
        return DeconvolutedPeakSet(p.clone() for p in self)

    def between(self, m1, m2, tolerance=1e-5):
        p1 = self.has_peak(m1, tolerance)
        p2 = self.has_peak(m2, tolerance)
        return self[p1.index:p2.index]


def _sweep_solution(array, value, lo, hi, tolerance, verbose=False):
    best_index = -1
    best_error = float('inf')
    for i in range(hi - lo):
        target = array[lo + i]
        error = ppm_error(value, target.neutral_mass)
        abs_error = abs(error)
        if abs_error < tolerance and abs_error < best_error:
            best_index = lo + i
            best_error = abs_error
    if best_index == -1:
        return None
    else:
        return array[best_index]


def _binary_search(array, value, lo, hi, tolerance, verbose=False):
    if (hi - lo) < 5:
        return _sweep_solution(array, value, lo, hi, tolerance, verbose)
    else:
        mid = (hi + lo) / 2
        target = array[mid]
        target_value = target.neutral_mass
        error = ppm_error(value, target_value)

        if abs(error) <= tolerance:
            return _sweep_solution(array, value, max(mid - 5, lo), min(mid + 5, hi), tolerance, verbose)
        elif target_value > value:
            return _binary_search(array, value, lo, mid, tolerance, verbose)
        elif target_value < value:
            return _binary_search(array, value, mid, hi, tolerance, verbose)
    raise Exception("No recursion found!")


def binary_search(array, value, tolerance=2e-5, verbose=False):
    size = len(array)
    if size == 0:
        return None
    return _binary_search(array, value, 0, size, tolerance, verbose)
