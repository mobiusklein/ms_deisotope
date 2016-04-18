import operator


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


def most_abundant_mz(envelope):
    return max([p for p in envelope if p.mz > 1], key=intensity_getter).mz


def average_mz(envelope):
    envelope = [p for p in envelope if p.mz > 1]
    return sum(map(mz_getter, envelope)) / float(len(envelope))


def average_signal_to_noise(envelope):
    envelope = [p for p in envelope if p.mz > 1]
    return sum(map(mz_getter, envelope)) / float(len(envelope))


class InterferenceDetection(object):
    def __init__(self, peaklist):
        self.peaklist = peaklist

    def detect_interference(self, envelope):
        min_peak = envelope.experimental[0]
        max_peak = envelope.experimental[-1]

        region = self.peaklist.between(
            min_peak.mz - min_peak.full_width_at_half_max,
            max_peak.mz + max_peak.full_width_at_half_max)

        included_intensity = sum(p.intensity for p in envelope.experimental)
        region_intensity = sum(p.intensity for p in region)

        score = 1 - (included_intensity / region_intensity)
        return score
