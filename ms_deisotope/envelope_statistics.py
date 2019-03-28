'''A collection of simple statistics describing isotopic envelopes.
These functions are primarily used during deconvolution to calculate
properties of :class:`~.DeconvolutedPeak` objects.
'''
import operator

from collections import namedtuple

from ms_deisotope.averagine import mass_charge_ratio

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


_CoIsolation = namedtuple("CoIsolation", ("neutral_mass", "intensity", "charge"))


class CoIsolation(_CoIsolation):
    '''Records the properties of a co-isolating ion associated with a primary ion.

    Attributes
    ----------
    neutral_mass: float
        The neutral mass of the ion
    intensity: float
        The intensity of the ion
    charge: int
        The charge of the ion
    mz: float
        The calculated m/z of the ion
    '''

    @property
    def mz(self):
        '''The calculated m/z of the ion

        Returns
        -------
        float
        '''
        return mass_charge_ratio(self.neutral_mass, self.charge)


class PrecursorPurityEstimator(object):
    """Estimate the contamination of precursor ions.

    Attributes
    ----------
    lower_extension : float
        An extra interval below the isolation window to search for coisolating ions,
        intended to collect trailing envelope peaks of lower m/z ions.
    default_width : float
        A default isolation window radius to use when an actual isolation window was
        not available.
    """
    def __init__(self, lower_extension=1.5, default_width=1.5):
        self.lower_extension = lower_extension
        self.default_width = default_width

    def precursor_purity(self, scan, precursor_peak, isolation_window=None):
        """Calculate the purity of an isolation window.

        If the isolation window is not given, :attr:`default_width` will be used.

        Parameters
        ----------
        scan : :class:`Scan`
            The precursor scan
        precursor_peak : :class:`~.DeconvolutedPeak`
            The precursor peak that was deconvolved, whose
            :attr:`~.DeconvolutedPeak.envelope` contains the signal assigned
            to this precursor ion.
        isolation_window : :class:`~.IsolationWindow`
            The isolation window to search within

        Returns
        -------
        float:
            The precursor ion purity, a value of 1.0 indicates that all the signal
            within the isolation window belonged to the assigned peak, and as the value
            approaches 0, the less signal in the isolation window belonged to the precursor
            ion.
        """
        peak_set = scan.peak_set
        mz = precursor_peak.mz
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
        purity = max(min((assigned / total), 1.0), 0)
        return purity

    def coisolation(self, scan, precursor_peak, isolation_window=None, relative_intensity_threshold=0.1,
                    ignore_singly_charged=False):
        """Find any deconvolution solutions which may have partially overlapped the isolation window.

        If the isolation window is not given, :attr:`default_width` will be used.

        Parameters
        ----------
        scan : :class:`Scan`
            The precursor scan
        precursor_peak : :class:`~.DeconvolutedPeak`
            The precursor peak that was deconvolved, whose
            :attr:`~.DeconvolutedPeak.envelope` contains the signal assigned
            to this precursor ion.
        isolation_window : :class:`~.IsolationWindow`
            The isolation window to search within
        relative_intensity_threshold: :class:`float`
            The minimum intensity fraction of `precursor_peak` that another peak must exceed
            in order to be considered.
        ignore_singly_charged: :class:`bool`
            Whether or not to omit singly charged peaks around the isolation window

        Returns
        -------
        :class:`list` of :class:`CoIsolation`
            A list of :class:`CoIsolation` instances describing the deconvoluted peaks
            whose envelopes overlap the isolation window and pass the thresholds.
        """
        peak_set = scan.deconvoluted_peak_set
        mz = precursor_peak.mz

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

    def __call__(self, scan, precursor_peak, isolation_window):
        """A convenience wrapper that calls both :meth:`precursor_purity` and
        :meth:`coisolation`.

        Parameters
        ----------
        scan : :class:`Scan`
            The precursor scan
        precursor_peak : :class:`~.DeconvolutedPeak`
            The precursor peak that was deconvolved, whose
            :attr:`~.DeconvolutedPeak.envelope` contains the signal assigned
            to this precursor ion.
        isolation_window : :class:`~.IsolationWindow`
            The isolation window to search within

        Returns
        -------
        purity: :class:`float`
            The return value of :meth:`precursor_purity`
        coisolation: :class:`list` of :class:`CoIsolation`
            The return value of :meth:`coisolation`

        See Also
        --------
        :meth:`precursor_purity`
        :meth:`coisolation`
        """
        purity = self.precursor_purity(scan, precursor_peak, isolation_window)
        coisolation = self.coisolation(scan, precursor_peak, isolation_window)
        return purity, coisolation
