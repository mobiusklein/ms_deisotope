'''Deconvolution strategies using a list of compositions.
'''

from ms_deisotope.averagine import (
    PROTON, isotopic_variants,
    TheoreticalIsotopicPattern,
    neutral_mass)

from ms_deisotope.envelope_statistics import average_mz, a_to_a2_ratio, most_abundant_mz
from ms_deisotope.constants import (
    IGNORE_BELOW,
    TRUNCATE_AFTER,
    ERROR_TOLERANCE,
    MAX_ITERATION,
    CONVERGENCE)

from ms_deisotope.peak_dependency_network import PeakDependenceGraph
from ms_deisotope.peak_set import DeconvolutedPeakSolution, DeconvolutedPeakSet
from ms_deisotope.scoring import IsotopicFitRecord

from .base import (DeconvoluterBase)


from .utils import (
    count_placeholders, prepare_peaklist,
    drop_placeholders, first_peak, mean, charge_range_,
    )


class CompositionListDeconvoluterBase(DeconvoluterBase):
    """A mixin class to provide common features for deconvoluters which process spectra
    using a list of targeted compositions.

    Attributes
    ----------
    composition_list : list of :class:`~.Mapping`
        A series of objects which represent elemental compositions and support
        the :class:`~.Mapping` interface to access their individual elements.
    """

    def __init__(self, composition_list, *args, **kwargs):
        self.composition_list = list(composition_list)
        self.incremental_truncation = kwargs.get(
            "incremental_truncation", None)
        super(CompositionListDeconvoluterBase, self).__init__(*args, **kwargs)

    def generate_theoretical_isotopic_cluster(self, composition, charge, truncate_after=TRUNCATE_AFTER,
                                              mass_shift=None, charge_carrier=PROTON,
                                              ignore_below=IGNORE_BELOW):
        """Generate a theoretical isotopic pattern for ``composition``

        Parameters
        ----------
        composition : :class:`~.Mapping`
            An object representing an elemental composition
        charge : int
            The charge state to generate the isotopic pattern for
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift : float, optional
            An arbitrary mass shift to apply to the generated theoretical isotopic pattern,
            moving all peaks forward by that mass charge ratio transformed mass.
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to |PROTON|

        Returns
        -------
        :class:`~.TheoreticalIsotopicPattern`
            The theoretical isotopic pattern generated
        """
        tid = isotopic_variants(
            composition, charge=charge, charge_carrier=charge_carrier)
        tid = TheoreticalIsotopicPattern(tid, tid[0].mz)
        tid.truncate_after(truncate_after)
        tid.ignore_below(ignore_below)
        if mass_shift is not None:
            tid.shift(tid[0].mz + mass_shift / abs(charge))
        return tid

    def recalibrate_theoretical_mz(self, theoretical_distribution, experimental_mz):
        """Recalibrate the m/z of the theoretical isotopic pattern to start from the
        peak matching the experimental monoisotopic m/z

        Parameters
        ----------
        theoretical_distribution : :class:`TheoreticalIsotopicPattern`
            The theoretical isotopic pattern to adjust
        experimental_mz : float
            The experimental monoisotopic peak m/z

        Returns
        -------
        TheoreticalIsotopicPattern
        """
        theoretical_distribution.shift(experimental_mz)
        return theoretical_distribution

    def fit_composition_at_charge(self, composition, charge, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON,
                                  truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW, mass_shift=None):
        """Produce an isotopic fit for `composition` at `charge` against the experimental peak set.

        This method requires that the instance also possess a method named `match_theoretical_isotopic_distribution`
        such as the one implemented in :class:`DeconvoluterBase`.

        Parameters
        ----------
        composition : :class:`~.Mapping`
            An object representing an elemental composition
        charge : int
            The charge state to generate the isotopic pattern for
        error_tolerance : float
            The mass accuracy required to for peak matches
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift : float, optional
            An arbitrary mass shift to apply to the generated theoretical isotopic pattern,
            moving all peaks forward by that mass charge ratio transformed mass.
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to |PROTON|

        Returns
        -------
        :class:`~.IsotopicFitRecord`
        """
        tid = self.generate_theoretical_isotopic_cluster(composition, charge=charge, truncate_after=truncate_after,
                                                         charge_carrier=charge_carrier, ignore_below=ignore_below,
                                                         mass_shift=mass_shift)
        monoisotopic_peak = self.peaklist.has_peak(tid[0].mz, error_tolerance)
        if monoisotopic_peak is not None:
            tid = self.recalibrate_theoretical_mz(tid, monoisotopic_peak.mz)
        eid = self.match_theoretical_isotopic_distribution(
            tid.peaklist, error_tolerance)

        missed_peaks = count_placeholders(eid)

        if missed_peaks > len(eid) / 2:
            return None

        self.scale_theoretical_distribution(tid, eid)
        score = self.scorer.evaluate(self.peaklist, eid, tid.peaklist)
        fit = IsotopicFitRecord(None, score, charge, tid, eid)
        fit.missed_peaks = missed_peaks
        return fit

    def _make_deconvoluted_peak_solution(self, fit, composition, charge_carrier):
        eid = fit.experimental
        tid = fit.theoretical
        charge = fit.charge
        rep_eid = drop_placeholders(eid)
        total_abundance = sum(
            p.intensity for p in eid if p.intensity > 1)

        monoisotopic_mass = neutral_mass(
            tid.monoisotopic_mz, charge, charge_carrier)
        monoisotopic_mz = tid.monoisotopic_mz

        reference_peak = first_peak(eid)
        peak = DeconvolutedPeakSolution(
            composition, fit,
            monoisotopic_mass, total_abundance, charge,
            signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
            index=reference_peak.index,
            full_width_at_half_max=mean(
                p.full_width_at_half_max for p in rep_eid),
            a_to_a2_ratio=a_to_a2_ratio(tid),
            most_abundant_mass=neutral_mass(
                most_abundant_mz(eid), charge),
            average_mass=neutral_mass(average_mz(eid), charge),
            score=fit.score,
            envelope=[(p.mz, p.intensity) for p in rep_eid],
            mz=monoisotopic_mz, area=sum(e.area for e in eid))
        return peak

    def deconvolute_composition(self, composition, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW,
                                mass_shift=None):
        """For each charge state under consideration, fit the theoretical isotopic pattern for this composition,
        and if the fit is satisfactory, add it to the results set.

        Parameters
        ----------
        composition : :class:`~.Mapping`
            An object representing an elemental composition
        error_tolerance : float
            The mass accuracy required to for peak matches
        charge_range : tuple
            The charge state range to generate the isotopic patterns for
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift : float, optional
            An arbitrary mass shift to apply to the generated theoretical isotopic pattern,
            moving all peaks forward by that mass charge ratio transformed mass.
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to |PROTON|
        """
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(composition, charge=charge, error_tolerance=error_tolerance,
                                                 truncate_after=truncate_after, charge_carrier=charge_carrier,
                                                 mass_shift=mass_shift, ignore_below=ignore_below)
            if fit is None:
                continue
            if not self.scorer.reject(fit):
                eid = fit.experimental
                tid = fit.theoretical
                rep_eid = drop_placeholders(eid)
                if (len(rep_eid) < 2) or (len(rep_eid) < (len(tid) / 2.)) or (len(rep_eid) == 1 and fit.charge > 1):
                    continue
                if self.incremental_truncation is not None:
                    fits = [fit]
                    for case in self.fit_incremental_truncation(fit, self.incremental_truncation):
                        if not self.scorer.reject(case):
                            fits.append(case)
                    fit = self.scorer.select.best(fits)
                peak = self._make_deconvoluted_peak_solution(
                    fit, composition, charge_carrier)
                self._deconvoluted_peaks.append(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)


class CompositionListDeconvoluter(CompositionListDeconvoluterBase):
    '''Fit exact isotopic patterns from a list of compositions.

    Fits are accepted as they are made, making this algorithm unsuitable for
    complex spectra where isotopic patterns will share peaks.

    Attributes
    ----------
    composition_list : list of :class:`~.Mapping`
        A series of objects which represent elemental compositions and support
        the :class:`~.Mapping` interface to access their individual elements.
    peaklist : :class:`~ms_peak_picker.PeakSet`
        The collection of :class:`~.ms_peak_picker.FittedPeak` instances and possible associated
        data to deconvolute.
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
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
    '''

    def __init__(self, peaklist, composition_list, scorer,
                 use_subtraction=False, scale_method='sum',
                 verbose=False, use_quick_charge=False, **kwargs):
        self.peaklist = prepare_peaklist(peaklist)
        self.scorer = scorer
        self.verbose = verbose
        self._deconvoluted_peaks = []
        self.use_quick_charge = use_quick_charge
        super(CompositionListDeconvoluter, self).__init__(
            composition_list,
            use_subtraction=use_subtraction, scale_method=scale_method,
            merge_isobaric_peaks=True, **kwargs)

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), charge_carrier=PROTON,
                    truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW, mass_shift=None, **kwargs):
        """Deconvolute the spectrum, extracting isotopic patterns from the composition list.

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift: float, optional
            An optional mass shift to apply to each composition

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`
        """
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance=error_tolerance,
                                         charge_range=charge_range, charge_carrier=charge_carrier,
                                         truncate_after=truncate_after, ignore_below=ignore_below,
                                         mass_shift=mass_shift)
        return DeconvolutedPeakSet(self._deconvoluted_peaks).reindex()


class CompositionListPeakDependenceGraphDeconvoluter(CompositionListDeconvoluter):
    '''Fit exact isotopic patterns from a list of compositions.

    Fits are added to a peak dependence graph, and the best fit is chosen after
    all fits are calculated at each iteration.

    Attributes
    ----------
    composition_list : list of :class:`~.Mapping`
        A series of objects which represent elemental compositions and support
        the :class:`~.Mapping` interface to access their individual elements.
    peaklist : :class:`~ms_peak_picker.PeakSet`
        The collection of ms_peak_picker.FittedPeak instances and possible associated
        data to deconvolute.
    scorer : :class:`~.IsotopicFitterBase`
        The criterion for evaluating individual isotopic pattern fits
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : :class:`~PeakDependenceGraph`
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
    '''

    def __init__(self, peaklist, composition_list, scorer,
                 use_subtraction=False, scale_method='sum',
                 verbose=False, use_quick_charge=False, **kwargs):
        max_missed_peaks = kwargs.get("max_missed_peaks", 1)
        super(CompositionListPeakDependenceGraphDeconvoluter, self).__init__(
            peaklist, composition_list, scorer=scorer, use_subtraction=use_subtraction,
            scale_method=scale_method, verbose=verbose, use_quick_charge=use_quick_charge,
            **kwargs)

        self.peak_dependency_network = PeakDependenceGraph(
            self.peaklist, maximize=self.scorer.is_maximizing(), **kwargs)
        self.max_missed_peaks = max_missed_peaks

    @property
    def max_missed_peaks(self):
        """The maximum number of missed peaks per isotopic fit record permitted.

        This property directly mirrors :attr:`PeakDependenceGraph.max_missed_peaks`

        Returns
        -------
        int
        """
        return self.peak_dependency_network.max_missed_peaks

    @max_missed_peaks.setter
    def max_missed_peaks(self, value):
        self.peak_dependency_network.max_missed_peaks = value

    def _save_peak_solution(self, solution):
        self._deconvoluted_peaks.append(solution)

    def deconvolute_composition(self, composition, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW,
                                mass_shift=None):
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(
                composition, charge, error_tolerance, charge_carrier=charge_carrier,
                truncate_after=truncate_after, mass_shift=mass_shift, ignore_below=ignore_below)
            if fit is None:
                continue
            rep_eid = drop_placeholders(fit.experimental)
            if len(rep_eid) == 1 and fit.charge > 1:
                continue
            if not self.scorer.reject(fit):
                self.peak_dependency_network.add_fit_dependence(fit)
            if self.incremental_truncation is not None:
                for case in self.fit_incremental_truncation(fit, self.incremental_truncation):
                    if not self.scorer.reject(case):
                        self.peak_dependency_network.add_fit_dependence(case)


    def populate_graph(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), truncate_after=TRUNCATE_AFTER,
                       charge_carrier=PROTON, ignore_below=IGNORE_BELOW, mass_shift=None):
        """For each composition, for each charge state under consideration, fit the theoretical
        isotopic pattern for this composition, and if the fit is satisfactory, add it to the
        peak dependence graph for later selecting the optimal solution.

        Parameters
        ----------
        error_tolerance : float
            The mass accuracy required to for peak matches
        charge_range : tuple
            The charge state range to generate the isotopic patterns for
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to |PROTON|
        mass_shift : float, optional
            An arbitrary mass shift to apply to the generated theoretical isotopic pattern,
            moving all peaks forward by that mass charge ratio transformed mass.
        """
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance, charge_range,
                                         truncate_after=truncate_after, charge_carrier=charge_carrier,
                                         ignore_below=ignore_below, mass_shift=mass_shift)

    def select_best_disjoint_subgraphs(self, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Construct connected envelope graphs from :attr:`peak_dependency_network` and
        extract the best disjoint isotopic pattern fits in each envelope graph. This in
        turn produces one or more :class:`~.DeconvolutedPeakSolution` instances from each
        disjoint fit, which are processed and added to the results set.

        Parameters
        ----------
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result
        """
        disjoint_envelopes = self.peak_dependency_network.find_non_overlapping_intervals()

        for cluster in disjoint_envelopes:
            for fit in cluster.disjoint_best_fits():
                eid = fit.experimental
                tid = fit.theoretical
                composition = fit.data
                rep_eid = drop_placeholders(eid)
                if len(rep_eid) < 2 or len(rep_eid) < len(tid) / 2.:
                    continue

                peak = self._make_deconvoluted_peak_solution(
                    fit, composition, charge_carrier)
                self._save_peak_solution(peak)
                if self.use_subtraction:
                    self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), iterations=MAX_ITERATION,   # pylint: disable=arguments-differ
                    truncate_after=TRUNCATE_AFTER, charge_carrier=PROTON, ignore_below=IGNORE_BELOW,
                    mass_shift=None, convergence=CONVERGENCE, **kwargs):
        """Deconvolute the spectrum, extracting isotopic patterns from the composition list.

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        mass_shift: float, optional
            An optional mass shift to apply to each composition
        convergence : float, optional
            The threshold of the below which after the `(sum(intensity_before) - sum(
            intensity_after)) / sum(intensity_after)`

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`
        """
        if not self.use_subtraction:
            iterations = 1
        begin_signal = sum([p.intensity for p in self.peaklist])
        for _ in range(iterations):
            self.populate_graph(error_tolerance, charge_range, charge_carrier=charge_carrier,
                                truncate_after=truncate_after, ignore_below=ignore_below,
                                mass_shift=mass_shift)
            self.select_best_disjoint_subgraphs(error_tolerance)
            self._slice_cache.clear()
            end_signal = sum([p.intensity for p in self.peaklist]) + 1
            if (begin_signal - end_signal) / end_signal < convergence:
                break
            begin_signal = end_signal
        return DeconvolutedPeakSet(self._deconvoluted_peaks).reindex()
