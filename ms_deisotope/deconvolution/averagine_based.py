# -*- coding: utf-8 -*-
'''Implementations of averagine-based deconvoluters.

Averagine-based deconvoluters use an "average monomer" isotopic model
to interpolate the isotopic patterns for any peak in the experimental
spectrum. The term "averagine" comes from the name Senko gave to the
"average amino acid" when introducing the concept in [1].

References
----------
[1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic
    masses and ion populations for large biomolecules from resolved isotopic distributions.
    Journal of the American Society for Mass Spectrometry, 6(4), 229–233.
    http://doi.org/10.1016/1044-0305(95)00017-8
'''

from ms_deisotope.averagine import PROTON, AveragineCache, peptide, glycopeptide, glycan
from ms_deisotope.constants import IGNORE_BELOW, TRUNCATE_AFTER, SCALE_METHOD

from ms_deisotope.scoring import penalized_msdeconv

from .base import (
    DeconvoluterBase)

from .exhaustive import (
    ExhaustivePeakSearchDeconvoluterBase,
    PeakDependenceGraphDeconvoluterBase)

from .utils import count_placeholders, prepare_peaklist


class AveragineDeconvoluterBase(DeconvoluterBase):
    """A base class derived from :class:`DeconvoluterBase` which provides some common methods
    for fitting isotopic patterns using an Averagine model.
    """

    def __init__(self, use_subtraction=False, scale_method=SCALE_METHOD, merge_isobaric_peaks=True,
                 minimum_intensity=5., *args, **kwargs):
        super(AveragineDeconvoluterBase, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)

    def fit_theoretical_distribution(self, peak, error_tolerance, charge, charge_carrier=PROTON, truncate_after=0.8,
                                     ignore_below=IGNORE_BELOW):
        """Fit an isotopic pattern seeded at `peak` at `charge` charge.

        Generates a theoretical isotopic pattern using :attr:`averagine`, calls
        :meth:`match_theoretical_isotopic_distribution`
        to extract experimental peaks matching this theoretical pattern, scales the theoretical distribution using
        :meth:`scale_theoretical_distribution`, and evaluates the quality of the fit using :attr:`scorer`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The putative monoisotopic peak to use for interpolating an isotopic pattern
        error_tolerance : float
            Parts-per-million error tolerance for isotopic pattern matching
        charge : int
            The charge state to produce an isotopic pattern for
        charge_carrier : float, optional
            The charge carrier mass, defaults to |PROTON|

        Returns
        -------
        :class:`~.IsotopicFitRecord`
            The fitted isotopic pattern
        """
        tid = self.averagine.isotopic_cluster(
            peak.mz, charge, charge_carrier=charge_carrier,
            truncate_after=truncate_after, ignore_below=ignore_below)
        eid = self.match_theoretical_isotopic_distribution(
            tid, error_tolerance=error_tolerance)
        record = self._evaluate_theoretical_distribution(
            eid, tid, peak, charge)
        return record

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON, truncate_after=0.8,
                              ignore_below=IGNORE_BELOW):
        """Given a set of candidate monoisotopic peaks and charge states, and a PPM error tolerance,
        fit each putative isotopic pattern.

        Calls :meth:`fit_theoretical_distribution` on each candidate.

        If a fit does not satisfy :attr:`scorer` `.reject`, it is discarded. If a fit has only one real peak
        and has a charge state greater than 1, it will also be discarded.

        Parameters
        ----------
        peak_charge_set : set
            The set of candidate (:class:`~.FittedPeak`, charge) tuples to try to fit
        error_tolerance : float
            Matching error tolerance
        charge_carrier : float, optional
            The charge carrier to use. Defaults to |PROTON|

        Returns
        -------
        set
            The set of :class:`~.IsotopicFitRecord` instances produced
        """
        results = []
        for peak, charge in peak_charge_set:
            if peak.mz < 1:
                continue
            fit = self.fit_theoretical_distribution(
                peak, error_tolerance, charge, charge_carrier, truncate_after,
                ignore_below)
            fit.missed_peaks = count_placeholders(fit.experimental)
            if not self._check_fit(fit):
                continue
            results.append(fit)
            if self.incremental_truncation is not None:
                results.extend(self.fit_incremental_truncation(
                    fit, self.incremental_truncation))
        return set(results)


try:
    from ms_deisotope._c.deconvoluter_base import AveragineDeconvoluterBase
except ImportError:
    pass


class AveragineDeconvoluter(AveragineDeconvoluterBase, ExhaustivePeakSearchDeconvoluterBase):
    """A Deconvoluter which uses an :title-reference:`averagine` [1] model to generate theoretical
    isotopic patterns for each peak to consider. Combines :class:`AveragineDeconvoluterBase` and
    :class:`ExhaustivePeakSearchDeconvoluterBase` to create a working Deconvoluter type.


    Attributes
    ----------
    averagine : :class:`~.AveragineCache`
        The averagine model and associated theoretical isotopic pattern cache to use
        to build theoretical isotopic patterns.
    peaklist : :class:`~.PeakSet`
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    verbose : bool
        How much diagnostic information to provide

    References
    ----------
    [1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
        for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
        Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8
    """

    def __init__(self, peaklist, averagine=None, scorer=penalized_msdeconv,
                 use_subtraction=True, scale_method=SCALE_METHOD,
                 verbose=False, **kwargs):
        if averagine is None:
            averagine = AveragineCache(peptide, dict())
        else:
            if not isinstance(averagine, AveragineCache):
                averagine = AveragineCache(averagine, dict())
        self.peaklist = prepare_peaklist(peaklist)
        self.averagine = averagine
        self.scorer = scorer
        self._deconvoluted_peaks = []
        self.verbose = verbose

        super(AveragineDeconvoluter, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks=True, **kwargs)


class MultiAveragineDeconvoluterBase(DeconvoluterBase):
    """A base class derived from :class:`DeconvoluterBase` which provides some common methods
    for fitting isotopic patterns using multiple Averagine models.
    """

    def fit_theoretical_distribution(self, peak, error_tolerance, charge, averagine, charge_carrier=PROTON, truncate_after=0.8,
                                     ignore_below=IGNORE_BELOW):
        """Fit an isotopic pattern seeded at `peak` at `charge` charge.

        Generates a theoretical isotopic pattern using :attr:`averagine`, calls
        :meth:`match_theoretical_isotopic_distribution`
        to extract experimental peaks matching this theoretical pattern, scales the theoretical distribution using
        :meth:`scale_theoretical_distribution`, and evaluates the quality of the fit using :attr:`scorer`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The putative monoisotopic peak to use for interpolating an isotopic pattern
        error_tolerance : float
            Parts-per-million error tolerance for isotopic pattern matching
        charge : int
            The charge state to produce an isotopic pattern for
        averagine : :class:`~.AveragineCache`
            The isotopic model to use for this fitting
        charge_carrier : float, optional
            The charge carrier mass, defaults to |PROTON|

        Returns
        -------
        :class:`~.IsotopicFitRecord`
            The fitted isotopic pattern
        """
        tid = averagine.isotopic_cluster(
            peak.mz, charge, charge_carrier=charge_carrier,
            truncate_after=truncate_after, ignore_below=ignore_below)
        eid = self.match_theoretical_isotopic_distribution(
            tid, error_tolerance=error_tolerance)
        record = self._evaluate_theoretical_distribution(
            eid, tid, peak, charge)
        return record

    def _fit_peaks_at_charges(self, peak_charge_set, error_tolerance, charge_carrier=PROTON,
                              truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        results = []
        for peak, charge in peak_charge_set:
            for averagine in self.averagines:
                if peak.mz < 1:
                    continue
                fit = self.fit_theoretical_distribution(
                    peak, error_tolerance, charge, averagine, charge_carrier=charge_carrier,
                    truncate_after=truncate_after, ignore_below=ignore_below)
                fit.missed_peaks = count_placeholders(fit.experimental)
                fit.data = averagine
                if not self._check_fit(fit):
                    continue
                results.append(fit)
                if self.incremental_truncation is not None:
                    results.extend(self.fit_incremental_truncation(
                        fit, self.incremental_truncation))
        return set(results)


try:
    from ms_deisotope._c.deconvoluter_base import MultiAveragineDeconvoluterBase
except ImportError:
    pass


class MultiAveragineDeconvoluter(MultiAveragineDeconvoluterBase, ExhaustivePeakSearchDeconvoluterBase):
    """A Deconvoluter which uses multiple :title-reference:`averagine` [1] model to generate theoretical
    isotopic patterns for each peak to consider. Combines :class:`MultiAveragineDeconvoluterBase` and
    :class:`ExhaustivePeakSearchDeconvoluterBase` to create a working Deconvoluter type.

    This differs from :class:`AveragineDeconvoluter`, in that it will produce multiple isotopic fits for
    each (peak, charge) pair. This is advantageous when the isotopic patterns produced by different models
    are sufficiently different enough that they will favor different peak sets.

    Attributes
    ----------
    averagine : list of :class:`~.ms_deisotope.averagine.AveragineCache`
        The averagine models and associated theoretical isotopic pattern caches to use
        to build theoretical isotopic patterns.
    peaklist : :class:`~.ms_peak_picker.PeakSet`
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : :class:`~.ms_deisotope.scoring.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    verbose : bool
        How much diagnostic information to provide

    References
    ----------
    [1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
        for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
        Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8
    """

    def __init__(self, peaklist, averagines=None, scorer=penalized_msdeconv,
                 use_subtraction=True, scale_method=SCALE_METHOD,
                 merge_isobaric_peaks=True, minimum_intensity=5.,
                 verbose=False, *args, **kwargs):
        self.peaklist = prepare_peaklist(peaklist)
        self.scorer = scorer
        self.use_subtraction = use_subtraction
        self.scale_method = scale_method

        cache_backend = dict
        if averagines is None:
            averagines = [peptide, glycopeptide, glycan]
        averagines = [
            AveragineCache(avg, backend=cache_backend()) if not isinstance(
                avg, AveragineCache) else avg
            for avg in averagines]
        self.averagines = averagines
        self.verbose = verbose

        self._deconvoluted_peaks = []

        super(MultiAveragineDeconvoluter, self).__init__(
            use_subtraction, scale_method, merge_isobaric_peaks,
            minimum_intensity, *args, **kwargs)


class AveraginePeakDependenceGraphDeconvoluter(AveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    """A Deconvoluter which uses an :title-reference:`averagine` [1] model to generate theoretical
    isotopic patterns for each peak to consider, using a peak dependence graph to solve complex mass
    spectra.

    Extends :class:`AveragineDeconvoluter` to include features from
    :class:`PeakDependenceGraphDeconvoluterBase` making it suitable for deconvoluting complex spectra where
    peak overlaps are common.

    Attributes
    ----------
    peaklist : :class:`~.PeakSet`
        The centroided mass spectrum to deconvolute
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    averagine : :class:`~.AveragineCache`
        The averagine model and associated theoretical isotopic pattern cache to use
        to build theoretical isotopic patterns.
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : :class:`~.PeakDependenceGraph`
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
    merge_isobaric_peaks : bool
        If multiple passes produce peaks with identical mass values,
        should those peaks be summed
    minimum_intensity : float
        Experimental peaks whose intensity is below this level will be ignored
        by peak querying methods
    scale_method : str
        The name of the method to use to scale theoretical isotopic pattern intensities
        to match the experimental isotopic pattern. For a description of options, see
        :meth:`~.TheoreticalIsotopicPattern.scale`.
    use_subtraction : bool
        Whether or not to apply a subtraction procedure to experimental peaks after they
        have been fitted. This is only necessary if the same signal may be examined multiple
        times as in a multi-pass method or when peak dependence is not considered
    verbose : bool
        Produce extra logging information

    References
    ----------
    [1] Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
        for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
        Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8
    """

    def __init__(self, peaklist, *args, **kwargs):
        super(AveraginePeakDependenceGraphDeconvoluter,
              self).__init__(peaklist, *args, **kwargs)


class MultiAveraginePeakDependenceGraphDeconvoluter(MultiAveragineDeconvoluter, PeakDependenceGraphDeconvoluterBase):
    """Extends :class:`MultiAveragineDeconvoluter` to include features from
    :class:`PeakDependenceGraphDeconvoluterBase` making it suitable for deconvoluting complex spectra where
    peak overlaps are common.

    Attributes
    ----------
    peaklist : :class:`~.ms_peak_picker.PeakSet`
        The centroided mass spectrum to deconvolute
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    averagine : list of :class:`~.ms_deisotope.averagine.AveragineCache`
        The averagine model and associated theoretical isotopic pattern cache to use
        to build theoretical isotopic patterns.
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : :class:`~.PeakDependenceGraph`
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
    merge_isobaric_peaks : bool
        If multiple passes produce peaks with identical mass values,
        should those peaks be summed
    minimum_intensity : float
        Experimental peaks whose intensity is below this level will be ignored
        by peak querying methods
    scale_method : str
        The name of the method to use to scale theoretical isotopic pattern intensities
        to match the experimental isotopic pattern. For a description of options, see
        :meth:`~.TheoreticalIsotopicPattern.scale`.
    use_subtraction : bool
        Whether or not to apply a subtraction procedure to experimental peaks after they
        have been fitted. This is only necessary if the same signal may be examined multiple
        times as in a multi-pass method or when peak dependence is not considered
    verbose : bool
        Produce extra logging information

    """

    def __init__(self, peaklist, *args, **kwargs):
        super(MultiAveraginePeakDependenceGraphDeconvoluter,
              self).__init__(peaklist, *args, **kwargs)
