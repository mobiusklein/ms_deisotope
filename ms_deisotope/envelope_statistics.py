import operator

from collections import namedtuple

from ms_deisotope.averagine import mass_charge_ratio

intensity_getter = operator.attrgetter("intensity")
mz_getter = operator.attrgetter("mz")
snr_getter = operator.attrgetter("signal_to_noise")


def a_to_a2_ratio(envelope):
    if len(envelope) < 3:
        return 0.
    a0 = envelope[0]
    a2 = envelope[2]
    if a0.mz < 0 or a2.mz < 0:
        return 0.
    return a0.intensity / a2.intensity


def total_area(envelope):
    return sum(p.area for p in envelope)


def most_abundant_mz(envelope):
    return max([p for p in envelope if p.mz > 1], key=intensity_getter).mz


def average_mz(envelope):
    envelope = [p for p in envelope if p.mz > 1]
    return weighted_average(map(mz_getter, envelope), map(intensity_getter, envelope))


def average_signal_to_noise(envelope):
    envelope = [p for p in envelope if p.mz > 1]
    return sum(map(snr_getter, envelope)) / float(len(envelope))


def weighted_average(values, weights):
    # We need to traverse weights twice, therefore convert it to a list
    # as a generator (i.e. py3's map) would always have sum(weights) == 0
    weights = list(weights)
    return sum(v * w for v, w in zip(values, weights)) / float(sum(weights))


_CoIsolation = namedtuple("CoIsolation", ("neutral_mass", "intensity", "charge"))


class CoIsolation(_CoIsolation):
    @property
    def mz(self):
        return mass_charge_ratio(self.neutral_mass, self.charge)


class PrecursorPurityEstimator(object):

    def __init__(self, lower_extension=1.5, default_width=1.5):
        self.lower_extension = lower_extension
        self.default_width = default_width

    def precursor_purity(self, scan, precursor_peak):
        peak_set = scan.peak_set
        mz = precursor_peak.mz
        isolation_window = scan.isolation_window
        if isolation_window is None:
            lower_bound = mz - self.default_width
            upper_bound = mz + self.default_width
        else:
            lower_bound = isolation_window.lower_bound
            upper_bound = isolation_window.upper_bound
        envelope = precursor_peak.envelope
        assigned = sum([p.intensity for p in envelope])
        total = sum([p.intensity for p in peak_set.between(lower_bound, upper_bound)])
        if total == 0:
            return 0
        purity = 1 - (assigned / total)
        return purity

    def coisolation(self, scan, precursor_peak, relative_intensity_threshold=0.1, ignore_singly_charged=True):
        peak_set = scan.deconvoluted_peak_set
        mz = precursor_peak.mz

        isolation_window = scan.isolation_window
        if isolation_window is None:
            lower_bound = mz - self.default_width
            upper_bound = mz + self.default_width
        else:
            lower_bound = isolation_window.lower_bound
            upper_bound = isolation_window.upper_bound
        extended_lower_bound = lower_bound - self.lower_extension

        peaks = peak_set.between(extended_lower_bound, upper_bound, use_mz=True)
        others = [
            CoIsolation(p.neutral_mass, p.intensity, p.charge)
            for p in peaks if p != precursor_peak and p.intensity > (
                precursor_peak.intensity * relative_intensity_threshold) and
            p.envelope[-1].mz > lower_bound and ((abs(p.charge) != 1 and ignore_singly_charged) or
                                                 not ignore_singly_charged)
        ]
        return others

    def __call__(self, scan, precursor_peak):
        purity = self.precursor_purity(scan, precursor_peak)
        coisolation = self.coisolation(scan, precursor_peak)
        return purity, coisolation
