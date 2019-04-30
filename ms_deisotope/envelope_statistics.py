'''A collection of simple statistics describing isotopic envelopes.
These functions are primarily used during deconvolution to calculate
properties of :class:`~.DeconvolutedPeak` objects.
'''
import operator

intensity_getter = operator.attrgetter("intensity")
mz_getter = operator.attrgetter("mz")
snr_getter = operator.attrgetter("signal_to_noise")


def a_to_a2_ratio(envelope):
    """Calculate the ratio of the abundance of the A+0 (monoisotopic)
    peak to the abundance of the A+2 (+2 neutron) peak.

    Parameters
    ----------
    envelope : :class:`list` of :class:`~.FittedPeak`
        The sequence of experimental peaks matched.

    Returns
    -------
    float
    """
    if len(envelope) < 3:
        return 0.
    a0 = envelope[0]
    a2 = envelope[2]
    if a0.mz < 0 or a2.mz < 0:
        return 0.
    return a0.intensity / a2.intensity


def total_area(envelope):
    """Sum the area of all peaks in the envelope.

    Parameters
    ----------
    envelope : :class:`list` of :class:`~.FittedPeak`
        The sequence of experimental peaks matched.

    Returns
    -------
    float
    """
    return sum(p.area for p in envelope)


def most_abundant_mz(envelope):
    """Find the most abundant m/z peak in the envelope

    Parameters
    ----------
    envelope : :class:`list` of :class:`~.FittedPeak`
        The sequence of experimental peaks matched.

    Returns
    -------
    float
    """
    return max([p for p in envelope if p.mz > 1], key=intensity_getter).mz


def average_mz(envelope):
    """An intensity weighted average the m/z of the peaks
    in the envelope

    Parameters
    ----------
    envelope : :class:`list` of :class:`~.FittedPeak`
        The sequence of experimental peaks matched.

    Returns
    -------
    float
    """
    envelope = [p for p in envelope if p.mz > 1]
    return weighted_average(map(mz_getter, envelope), map(intensity_getter, envelope))


def average_signal_to_noise(envelope):
    """Average the signal to noise ratio of the peaks in the envelope

    Parameters
    ----------
    envelope : :class:`list` of :class:`~.FittedPeak`
        The sequence of experimental peaks matched.

    Returns
    -------
    float
    """
    envelope = [p for p in envelope if p.mz > 1]
    return sum(map(snr_getter, envelope)) / float(len(envelope))


def weighted_average(values, weights):
    """Calculate a weighted average over `values` weighted by `weights`.

    Parameters
    ----------
    values : :class:`~.Iterable`
        The values
    weights : :class:`~.Iterable`
        The weights

    Returns
    -------
    float
    """
    # We need to traverse weights twice, therefore convert it to a list
    # as a generator (i.e. py3's map) would always have sum(weights) == 0
    weights = list(weights)
    return sum(v * w for v, w in zip(values, weights)) / float(sum(weights))
