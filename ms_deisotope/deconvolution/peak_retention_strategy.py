'''A set of strategies for retaining peaks following deconvolution to
recover from low quality spectra. These methods are not recommended for
use on spectra where most ions are multiply charged.

These objects are meant to be provided to :func:`~.deconvolut_peaks`'s
`retention_strategy` argument, which will take care of calling them with
the appropriate arguments.
'''
import operator
import abc

from six import add_metaclass

from ms_deisotope.averagine import neutral_mass

from .utils import from_fitted_peak

intensity_getter = operator.attrgetter("intensity")


@add_metaclass(abc.ABCMeta)
class PeakRetentionStrategyBase(object):
    """A strategy for retaining peaks left behind by the deconvolution procedure.

    Deconvolution requires that the isotopic peak structure is present in the peak list,
    which may not be consistent for all types of data. This type of method is appropriate
    when the caller expects real peaks may have been left behind, as is the case with smaller
    peptides and lower quality tandem mass spectra.
    """

    def create_peak(self, fitted_peak, charge=1):
        """Create a :class:`~.DeconvolutedPeak` from a given :class:`~.FittedPeak` ``peak``
        and an optional charge state.

        Parameters
        ----------
        fitted_peak : :class:`~.FittedPeak`
            The peak to treat as the monoisotopic peak of the unfitted deconvoluted peak
        charge : int, optional
            The charge state to specify the peak at (the default is 1)

        Returns
        -------
        DeconvolutedPeak
        """
        return from_fitted_peak(fitted_peak, charge)

    @abc.abstractmethod
    def retain_peaks(self, peaklist, original_peaklist=None, charge_range=None, solutions=None):
        """Given a list of :class:`~.FittedPeak` objects left over after deconvolution,
        produce a list of :class:`~.DeconvolutedPeak` objects from them according to some
        criteria.

        Deconvolution requires that the isotopic peak structure is present in the peak list,
        which may not be consistent for all types of data. This type of method is appropriate
        when the caller expects real peaks may have been left behind, as is the case with smaller
        peptides and lower quality tandem mass spectra.

        Parameters
        ----------
        peaklist : :class:`~.PeakSet`
            The list of peaks that remain after deconvolution
        original_peaklist : :class:`~.PeakSet`, optional
            The list of peaks that were initially presented for deconvolution,
            which may be used to infer relative thresholds from. If not provided,
            these thresholds may be learned from the left over peaks, but the
            parameters may not be as good.
        charge_range : :class:`tuple`, optional
            The range of charge states considered during deconvolution. Used to infer
            the minimum charge state to assign to the peaks this method creates. If
            not provided, the default charge state is :const`1`.
        solutions: :class:`~.DeconvolutedPeakSet`, optional
            The accepted solutions which might be used to infer/reject additional peaks.

        Returns
        -------
        :class:`list`:
            The list of :class:`~.DeconvolutedPeak` objects.
        """
        raise NotImplementedError()

    def __call__(self, peaklist, original_peaklist=None, charge_range=None, solutions=None):
        """A wrapper for :meth:`retain_peaks`. Given a list of :class:`~.FittedPeak`
        objects left over after deconvolution, produce a list of :class:`~.DeconvolutedPeak`
        objects from them according to some criteria.

        Parameters
        ----------
        peaklist : :class:`~.PeakSet`
            The list of peaks that remain after deconvolution
        original_peaklist : :class:`~.PeakSet`, optional
            The list of peaks that were initially presented for deconvolution,
            which may be used to infer relative thresholds from. If not provided,
            these thresholds may be learned from the left over peaks, but the
            parameters may not be as good.
        charge_range : :class:`tuple`, optional
            The range of charge states considered during deconvolution. Used to infer
            the minimum charge state to assign to the peaks this method creates. If
            not provided, the default charge state is :const`1`.
        solutions: :class:`~.DeconvolutedPeakSet`, optional
            The accepted solutions which might be used to infer/reject additional peaks.

        Returns
        -------
        :class:`list`:
            The list of :class:`~.DeconvolutedPeak` objects.

        See Also
        --------
        :meth:`retain_peaks`
        """
        return self.retain_peaks(peaklist, original_peaklist, charge_range)

    def infer_minimum_charge(self, charge_range):
        """Given a charge range :class:`tuple`, return the smallest absolute magnitude
        charge state possible in that range.

        This method is polarity aware, and does its reasoning on the absolute magnitude
        of the charge state range, rather than on its signed value.

        Parameters
        ----------
        charge_range : :class:`tuple`
            The range of charge values that were used during the deconvolution process

        Returns
        -------
        :class:`int`:
            The minimum charge state.
        """
        if charge_range is None:
            return 1
        abs_range = [abs(c) for c in charge_range]
        if abs_range[0] < abs_range[1]:
            min_charge = charge_range[0]
        else:
            min_charge = charge_range[1]
        return min_charge


class TopNRetentionStrategy(PeakRetentionStrategyBase):
    """This strategy retains up to at most :attr:`n_peaks` peaks from the
    leftover signal, and requires that any peaks it retains are atleast
    :attr:`base_peak_coefficient` * the base peak of the original peak list.

    This strategy treats the leftover peak as if it were the monoisotopic peak
    of an unobserved isotopic distribution. Because the monoisotopic peak ceases
    to dominate an isotopic distribution above a certain mass (as implied by an
    averagine isotopic model), any peak selected must have a neutral mass below
    :attr:`max_mass`.

    Attributes
    ----------
    n_peaks : int
        The maximum number of peaks to retain.
    base_peak_coefficient : float
        The fraction of the base peak intensity to threshold on. Defaults to :const:`0.05`
    max_mass : float
        The largest neutral mass that a peak may have before it would no longer have a dominant
        monoisotopic peak and the peak will be discarded.

    """
    def __init__(self, n_peaks=50, base_peak_coefficient=0.05, max_mass=850.0):
        self.n_peaks = n_peaks
        self.base_peak_coefficient = base_peak_coefficient
        self.max_mass = max_mass

    def retain_peaks(self, peaklist, original_peaklist=None, charge_range=None, solutions=None):
        if original_peaklist is None:
            base_peak_sequence = peaklist
        else:
            base_peak_sequence = original_peaklist
        try:
            base_peak = max([peak.intensity for peak in base_peak_sequence])
        except ValueError:
            return []
        min_charge = self.infer_minimum_charge(charge_range)
        min_charge /= abs(min_charge)
        threshold = self.base_peak_coefficient * base_peak
        peaklist = sorted(peaklist, key=intensity_getter, reverse=True)
        result = []
        for peak in peaklist:
            if neutral_mass(peak.mz, min_charge) > self.max_mass:
                continue
            if peak.intensity >= threshold:
                result.append(self.create_peak(peak, min_charge))
                if len(result) == self.n_peaks:
                    break
            else:
                break

        return result


simple_peak_retention = TopNRetentionStrategy()
