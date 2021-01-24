'''A collection of base classes for "exhaustive" search strategies.
These strategies attempt to assign *every* peak, but they are not meant
to be used on their own.

For complete implementations see :class:`~.AveragineDeconvoluter` and
:class:`~.AveraginePeakDependenceGraphDeconvoluter`.
'''

import operator

from ms_deisotope.constants import (
    ERROR_TOLERANCE, TRUNCATE_AFTER, IGNORE_BELOW, MAX_ITERATION,
    CONVERGENCE)

from ms_deisotope.peak_set import (
    DeconvolutedPeak,
    DeconvolutedPeakSet)

from ms_deisotope.averagine import PROTON, neutral_mass

from ms_deisotope.envelope_statistics import (
    average_mz,
    a_to_a2_ratio,
    most_abundant_mz)

from ms_deisotope.utils import (
    TrivialTargetedDeconvolutionResult)

from ms_deisotope.peak_dependency_network import (
    PeakDependenceGraph,
    NetworkedTargetedDeconvolutionResult)


from .base import DeconvoluterBase
from .utils import (
    ChargeIterator,
    has_previous_peak_at_charge,
    has_successor_peak_at_charge,
    drop_placeholders,
    count_placeholders,
    first_peak,
    mean,
    info,
    debug)


class ExhaustivePeakSearchDeconvoluterBase(DeconvoluterBase):
    """Provides common methods for algorithms which attempt to find a deconvolution for every peak
    in a spectrum. This assumes no dependence between different peaks, instead it relies on subtraction,
    breadth of search, and order of encounter to avoid artefactual fits. This is usually not reasonable,
    so instead please use this class's extension, :class:`PeakDependenceGraphDeconvoluterBase` which can
    express dependence of fits on common resources.

    This class is not meant to be instantiated, but instead used as a mixin for classes that also
    inherit from :class:`DeconvoluterBase` and provide methods `fit_theoretical_distribution`
    and `_fit_peaks_at_charges`

    Attributes
    ----------
    use_quick_charge: bool
        Whether or not to use the :title-reference:`QuickCharge` algorithm when generating
        putative charge states.
    incremental_truncation: float or :const:`None`
        If not :const:`None`, the isotopic pattern truncation lower bound to contract each fit to
        incrementally, generating new fits for each dropped peak using
        :meth:`~.DeconvoluterBase.fit_incremental_truncation`

    """

    def __init__(self, peaklist, *args, **kwargs):
        super(ExhaustivePeakSearchDeconvoluterBase,
              self).__init__(peaklist, *args, **kwargs)
        self.use_quick_charge = kwargs.get("use_quick_charge", False)
        self.incremental_truncation = kwargs.get("incremental_truncation", None)

    def _get_all_peak_charge_pairs(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                   left_search_limit=3, right_search_limit=3,
                                   recalculate_starting_peak=True, use_quick_charge=False):
        """Construct the set of all unique candidate (monoisotopic peak, charge state) pairs using
        the provided search parameters.

        The search is performed using :func:`has_previous_peak_at_charge`, :func:`has_successor_peak_at_charge`,
        :meth:`_find_previous_putative_peak`, and :meth:`_find_next_putative_peak`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        recalculate_starting_peak : bool, optional
            Whether or not to re-calculate the putative starting peak m/z based upon nearby
            peaks close to where isotopic peaks for `peak` should be. Defaults to True

        Returns
        -------
        set
            The set of all unique candidate (monoisotopic peak, charge state)
        """

        target_peaks = set()

        charges = ChargeIterator(*charge_range)

        if use_quick_charge:
            charges.sequence_from_quickcharge(self.peaklist, peak)

        for charge in charges:
            target_peaks.add((peak, charge))

            # Look Left
            for i in range(1, left_search_limit):
                prev_peak = has_previous_peak_at_charge(
                    self, peak, charge, i)
                if prev_peak is None:
                    continue
                target_peaks.add((prev_peak, charge))

                if recalculate_starting_peak:
                    target_peaks.update(self._find_previous_putative_peak(
                        peak.mz, charge, i, 2 * error_tolerance))

            # Look Right
            for i in range(1, right_search_limit):
                nxt_peak = has_successor_peak_at_charge(
                    self, peak, charge, i)
                if nxt_peak is None:
                    continue
                target_peaks.add((nxt_peak, charge))

                if recalculate_starting_peak:
                    target_peaks.update(self._find_next_putative_peak(
                        peak.mz, charge, i, 2 * error_tolerance))

            if recalculate_starting_peak:
                for i in range(min(left_search_limit, 2)):
                    target_peaks.update(self._find_next_putative_peak(
                        peak.mz, charge, step=i, tolerance=2 * error_tolerance))

        return target_peaks

    def _fit_all_charge_states(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), left_search_limit=3,
                               right_search_limit=3, recalculate_starting_peak=True, charge_carrier=PROTON,
                               truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Carry out the fitting process for `peak`.

        This method calls :meth:`_get_all_peak_charge_pairs` to collect all hypothetical solutions
        for `peak`, and invokes :meth:`_fit_peaks_at_charges` to evaluate them.

        The method :meth:`_fit_peaks_at_charges` is required by this interface, but is not defined by
        it, as it depends upon the underlying isotopic pattern fitting algorithm. See one of the
        Averagine-based algorithms for an implementation, such as :class:`AveragineDeconvoluterBase`,
        a complementary ancestor with this class.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        recalculate_starting_peak : bool, optional
            Whether or not to re-calculate the putative starting peak m/z based upon nearby
            peaks close to where isotopic peaks for `peak` should be. Defaults to True
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        set
            The set of :class:`~.IsotopicFitRecord` instances produced
        """
        target_peaks = self._get_all_peak_charge_pairs(
            peak, error_tolerance=error_tolerance,
            charge_range=charge_range,
            left_search_limit=left_search_limit,
            right_search_limit=right_search_limit,
            recalculate_starting_peak=recalculate_starting_peak,
            use_quick_charge=self.use_quick_charge)

        results = self._fit_peaks_at_charges(
            target_peaks, error_tolerance, charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)
        return (results)

    def charge_state_determination(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                                   left_search_limit=3, right_search_limit=3,
                                   charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                                   ignore_below=IGNORE_BELOW):
        """Determine the optimal isotopic fit for `peak`, extracting it's charge state and monoisotopic peak.

        This method invokes :meth:`_fit_all_charge_states`, and then uses :attr:`scorer`'s `select` method to
        choose the optimal solution.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.IsotopicFitRecord`
            The best scoring isotopic fit
        """
        results = self._fit_all_charge_states(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)

        if self.verbose:
            info("Fits for %r", peak)
            for rec in sorted(results)[-10:]:
                info(rec)
        try:
            result = self.scorer.select(results)
            return result
        except ValueError:
            return None

    def _make_deconvoluted_peak(self, fit, charge_carrier):
        score, charge, eid, tid = fit
        rep_eid = drop_placeholders(eid)
        total_abundance = sum(p.intensity for p in rep_eid)
        monoisotopic_mass = neutral_mass(
            tid.monoisotopic_mz, charge, charge_carrier)
        reference_peak = first_peak(eid)

        dpeak = DeconvolutedPeak(
            neutral_mass=monoisotopic_mass, intensity=total_abundance, charge=charge,
            signal_to_noise=mean(p.signal_to_noise for p in rep_eid),
            index=reference_peak.index,
            full_width_at_half_max=mean(
                p.full_width_at_half_max for p in rep_eid),
            a_to_a2_ratio=a_to_a2_ratio(tid),
            most_abundant_mass=neutral_mass(
                most_abundant_mz(eid), charge),
            average_mass=neutral_mass(average_mz(eid), charge),
            score=score,
            envelope=[(p.mz, p.intensity) for p in eid],
            mz=tid.monoisotopic_mz, fit=fit,
            area=sum(e.area for e in eid))
        return dpeak

    def deconvolute_peak(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                         left_search_limit=3, right_search_limit=3, charge_carrier=PROTON,
                         truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Perform a deconvolution for `peak`, generating a new :class:`ms_deisotope.peak_set.DeconvolutedPeak` instance
        corresponding to the optimal solution.

        This new peak has an m/z matching the monoisotopic peak of the pattern containing `peak`, and its intensity
        is the sum of all the matched peaks in its isotopic pattern. Its charge, isotopic fit, and other qualities
        are derived from the :class:`ms_deisotope.scoring.IsotopicFitRecord` instance corresponding to its best
        solution.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier, or more specifically, the moiety which is added for
            each incremental change in charge state. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.DeconvolutedPeak`
        """
        fit = self.charge_state_determination(
            peak, error_tolerance=error_tolerance,
            charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)
        if fit is None:
            return
        tid = fit.theoretical
        dpeak = self._make_deconvoluted_peak(fit, charge_carrier)
        self._deconvoluted_peaks.append(dpeak)
        if self.use_subtraction:
            self.subtraction(tid, error_tolerance)
        return dpeak

    def targeted_deconvolution(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                               left_search_limit=3, right_search_limit=3,
                               charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                               ignore_below=IGNORE_BELOW):
        """Express the intent that this peak's deconvolution solution will be retrieved at a later point in the process
        and that it should be deconvoluted, and return a handle to retrieve the results with.

        This algorithm's implementation is simple enough that this is equivalent to just performing the deconvolution
        now and storing the result in a :class:`~.TrivialTargetedDeconvolutionResult` instance.

        Otherwise identical to :meth:`deconvolute_peak`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.TrivialTargetedDeconvolutionResult`
        """
        dpeak = self.deconvolute_peak(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after,
            ignore_below=ignore_below)
        result = TrivialTargetedDeconvolutionResult(self, dpeak, peak)
        return result

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                    order_chooser=operator.attrgetter('index'),
                    left_search_limit=3, right_search_limit=3, charge_carrier=PROTON,
                    truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW,
                    **kwargs):
        """Completely deconvolute the spectrum.

        Visit each peak in the order chosen by `order_chooser`, and call :meth:`deconvolute_peak`
        on it with the provided arguments. This assumes all overlaps in isotopic pattern are captured
        by the search limits. This is usually not the case. For an alternative see
        :class:`PeakDependenceGraphDeconvoluterBase`

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        order_chooser : callable, optional:
            A callable used as a key function for sorting peaks into the order they will
            be visited during deconvolution. Defaults to :obj:`operator.attrgetter("index")`
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`
        """
        i = 0
        for peak in sorted(self.peaklist, key=order_chooser, reverse=True):
            if peak.mz < 2 or peak.intensity < self.minimum_intensity:
                continue
            self.deconvolute_peak(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit,
                right_search_limit=right_search_limit, charge_carrier=charge_carrier,
                truncate_after=truncate_after, ignore_below=ignore_below)
            i += 1

        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(
                self._deconvoluted_peaks)

        return DeconvolutedPeakSet(self._deconvoluted_peaks).reindex()


try:
    from ms_deisotope._c.deconvoluter_base import (
        _get_all_peak_charge_pairs as _c_get_all_peak_charge_pairs,
        _make_deconvoluted_peak as _c_make_deconvoluted_peak)
    ExhaustivePeakSearchDeconvoluterBase._get_all_peak_charge_pairs = _c_get_all_peak_charge_pairs
    ExhaustivePeakSearchDeconvoluterBase._make_deconvoluted_peak = _c_make_deconvoluted_peak
except ImportError as e:
    print(e)


class PeakDependenceGraphDeconvoluterBase(ExhaustivePeakSearchDeconvoluterBase):
    """Extends the concept of :class:`ExhaustivePeakSearchDeconvoluterBase` to include a way to handle
    conflicting solutions which claim the same experimental peak.

    This lets the Deconvoluter assign a single peak only once, and to the "best" solution to use it. To
    do this, the Deconvoluter constructs a graph where peaks are nodes, and isotopic fits are hyperedges
    connecting multiple nodes. Rather than deconvoluting the spectrum step by step, assigning signal as
    it explores the spectrum, the Deconvoluter instead inserts each considered isotopic fit into the graph.
    After completely traversing the spectrum, the Deconvoluter solves the dependence graph attempting to
    maximize some criterion and produce a set of disjoint isotopic fits. These fits are then assigned signal
    and added to the deconvoluted spectrum as normal.

    The criterion used is currently a greedy maximization (or minimization) of each connected component of
    the peak dependence graph.

    Attributes
    ----------
    max_missed_peaks : int
        The maximum number of missing peaks to tolerate in an isotopic fit
    peak_dependency_network : :class:`~.PeakDependenceGraph`
        The peak dependence graph onto which isotopic fit dependences on peaks
        are constructed and solved.
    """

    def __init__(self, peaklist, *args, **kwargs):
        max_missed_peaks = kwargs.get("max_missed_peaks", 1)
        self.subgraph_solver_type = kwargs.get("subgraph_solver", 'top')
        super(PeakDependenceGraphDeconvoluterBase,
              self).__init__(peaklist, *args, **kwargs)
        self.peak_dependency_network = PeakDependenceGraph(
            self.peaklist, maximize=self.scorer.is_maximizing())
        self.max_missed_peaks = max_missed_peaks
        self.fit_postprocessor = kwargs.get("fit_postprocessor", None)
        self._priority_map = {}

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

    def _explore_local(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER,
                       ignore_below=IGNORE_BELOW):
        """Given a peak, explore the local neighborhood for candidate isotopic fits and add each
        fit above a threshold to the peak dependence graph.

        The threshold assumes that a single peak's neighborhood will contain many, many fits, but
        that only the top `n` scoring fits are worth considering. For now, `n` is fixed at 100 or
        the half number of fits returned, whichever is larger. This is to prevent the fit graph
        from growing out of control and wasting time storing impractical fits. Any fit added to
        the graph will have to pass :attr:`scorer.select` as well, so weak fits will never be added,
        regardless of how many fits are allowed to be inserted.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 1
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 0
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.

        Returns
        -------
        int
            The number of fits added to the graph
        """
        results = self._fit_all_charge_states(
            peak, error_tolerance=error_tolerance, charge_range=charge_range, left_search_limit=left_search_limit,
            charge_carrier=charge_carrier, truncate_after=truncate_after, ignore_below=ignore_below)

        hold = set()
        for fit in results:
            if fit.charge > 1 and len(drop_placeholders(fit.experimental)) == 1:
                continue
            hold.add(fit)

        results = hold
        n = len(results)

        stop = max(min(n // 2, 100), 10)
        if n == 0:
            return 0
        i = 0
        for i in range(stop):
            if len(results) == 0:
                break
            candidate = self.scorer.select(results)
            if candidate is None:
                break
            self.peak_dependency_network.add_fit_dependence(candidate)
            results.discard(candidate)

        return i

    def populate_graph(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8), left_search_limit=1,
                       right_search_limit=0, charge_carrier=PROTON,
                       truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Visit each experimental peak and execute :meth:`_explore_local` on it with the provided
        parameters, populating the peak dependence graph with all viable candidates.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 1
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 0
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        """
        for peak in self.peaklist:
            if peak in self._priority_map or peak.intensity < self.minimum_intensity:
                if self.verbose:
                    debug("Skipping %r (%r)", peak,
                          peak.intensity < self.minimum_intensity)
                continue
            n = self._explore_local(
                peak, error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                charge_carrier=charge_carrier,
                truncate_after=truncate_after, ignore_below=ignore_below)
            if self.verbose:
                debug("Exlporing Area Around %r Yielded %d Fits", peak, n)

    def postprocess_fits(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                         charge_carrier=PROTON, *args, **kwargs):
        """Postprocesses fits before solving the peak dependence graph.

        Currently a no-op.

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        """
        if self.fit_postprocessor is None:
            return

    def select_best_disjoint_subgraphs(self, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Construct connected envelope graphs from :attr:`peak_dependency_network` and
        extract the best disjoint isotopic pattern fits in each envelope graph. This in turn
        produces one or more :class:`DeconvolutedPeak` instances from each disjoint fit,
        which are processed and added to the results set.

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
        i = 0
        if self.subgraph_solver_type == 'disjoint':
            solver = self._solve_subgraph_disjoint
        elif self.subgraph_solver_type == 'iterative':
            solver = self._solve_subgraph_iterative
        elif self.subgraph_solver_type == 'top':
            solver = self._solve_subgraph_top
        else:
            raise ValueError("Unknown solver type %r" %
                             (self.subgraph_solver_type, ))

        for cluster in disjoint_envelopes:
            solutions = solver(cluster, error_tolerance, charge_carrier)
            for dpeak in solutions:
                self.peak_dependency_network.add_solution(dpeak.fit, dpeak)
                self._deconvoluted_peaks.append(dpeak)
                i += 1

    def _solve_subgraph_top(self, cluster, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Given a :class:`~.DependenceCluster`, return the single best fit from the collection of
        co-dependent fits.

        Parameters
        ----------
        cluster : :class:`~.DependenceCluster`
            The connected subgraph whose nodes will be searched
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result

        Returns
        -------
        list of :class:`~DeconvolutedPeak`
            The solved deconvolution solutions
        """
        fit = cluster[0]
        _, _, eid, tid = fit
        rep_eid = drop_placeholders(eid)
        if len(rep_eid) == 0:
            return []
        dpeak = self._make_deconvoluted_peak(fit, charge_carrier)
        if self.use_subtraction:
            self.subtraction(tid, error_tolerance)
        return [dpeak]

    def _solve_subgraph_disjoint(self, cluster, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Given a :class:`~.DependenceCluster`, find a greedy disjoint set of isotopic fits.

        Parameters
        ----------
        cluster : :class:`~.DependenceCluster`
            The connected subgraph whose nodes will be searched
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result

        Returns
        -------
        list of :class:`~DeconvolutedPeak`
            The solved deconvolution solutions

        """
        disjoint_best_fits = cluster.disjoint_best_fits()
        i = 0
        solutions = []
        for fit in disjoint_best_fits:
            _, _, eid, tid = fit
            rep_eid = drop_placeholders(eid)
            if len(rep_eid) == 0:
                continue
            dpeak = self._make_deconvoluted_peak(fit, charge_carrier)
            solutions.append(dpeak)
            i += 1
            if self.use_subtraction:
                self.subtraction(tid, error_tolerance)
        return solutions

    def _solve_subgraph_iterative(self, cluster, error_tolerance=ERROR_TOLERANCE, charge_carrier=PROTON):
        """Given a :class:`~.DependenceCluster`, build a :class:`~.ConnectedSubgraph` and incrementally
        subtract the best fitting solution and update its overlapping envelopes.

        Parameters
        ----------
        cluster : :class:`~.DependencyCluster`
            The connected subgraph whose nodes will be searched
        error_tolerance : float, optional
            The error tolerance to use when performing subtraction, if subtraction is
            being performed.
        charge_carrier : float, optional
            The mass of the charge carrier as used for the deconvolution. Required to
            back-out the neutral mass of the deconvoluted result

        Returns
        -------
        list of :class:`~DeconvolutedPeak`
            The solved deconvolution solutions
        """
        print(cluster)
        subgraph = cluster.build_graph()
        solutions = []
        mask = set()
        best_node = subgraph[0]
        peak = self._make_deconvoluted_peak(best_node.fit, charge_carrier)
        solutions.append(peak)
        if self.use_subtraction:
            self.subtraction(best_node.fit.theoretical, error_tolerance)
        mask.add(best_node.index)
        print("Masking %d (%0.3f, %d)" % (
            best_node.index, best_node.fit.monoisotopic_peak.mz, best_node.fit.charge))

        overlapped_nodes = list(best_node.overlap_edges)
        maximize = subgraph.maximize

        n = len(subgraph)
        while len(mask) != n:
            best_node = None
            best_score = 0 if maximize else float('inf')
            retained = []
            for node in overlapped_nodes:
                if node.index in mask:
                    continue
                missed_peaks = node.fit.missed_peaks = count_placeholders(
                    node.fit.experimental)
                total_peaks = len(node.fit.experimental)
                invalid_peak_count = (missed_peaks >= total_peaks - 1 and abs(node.fit.charge) > 1
                                     ) or missed_peaks == total_peaks
                if invalid_peak_count or missed_peaks > self.max_missed_peaks:
                    mask.add(node.index)
                    print("Masking %d (%0.3f, %d). Invalidated Peak Count %r | Missed Peaks %d" % (
                        node.index, node.fit.monoisotopic_peak.mz,
                        node.fit.charge, invalid_peak_count, missed_peaks))

                    continue
                fit = node.fit
                fit.theoretical.normalize().scale(fit.experimental, self.scale_method)
                fit.score = self.scorer.evaluate(
                    self.peaklist, fit.experimental, fit.theoretical.peaklist)
                if self.scorer.reject(fit):
                    mask.add(node.index)
                    continue
                else:
                    retained.append(node)
                if maximize:
                    if fit.score > best_score:
                        best_node = node
                        best_score = fit.score
                else:
                    if fit.score < best_score:
                        best_node = node
                        best_score = fit.score

            if best_node is not None:
                peak = self._make_deconvoluted_peak(
                    best_node.fit, charge_carrier)
                solutions.append(peak)
                if self.use_subtraction:
                    self.subtraction(
                        best_node.fit.theoretical, error_tolerance)
                mask.add(best_node.index)
                print("Masking %d (%0.3f, %d)" % (
                    best_node.index, best_node.fit.monoisotopic_peak.mz, best_node.fit.charge))

                overlapped_nodes = [
                    node for node in best_node.overlap_edges if node.index not in mask]
            else:
                overlapped_nodes = []
            if not overlapped_nodes and len(mask) != n:
                overlapped_nodes = [
                    node for node in subgraph if node.index not in mask]
        return solutions

    def targeted_deconvolution(self, peak, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),
                               left_search_limit=3, right_search_limit=3, charge_carrier=PROTON,
                               truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW):
        """Express the intent that this peak's deconvolution solution will be retrieved at a later point in the process
        and that it should be deconvoluted, and return a handle to retrieve the results with.

        As the operation does not immediately result in a deconvoluted peak but just adds the resulting fits to
        :attr:`peak_dependency_network`, this method constructs an instance of
        :class:`~.NetworkedTargetedDeconvolutionResult` which holds all
        the required information for recovering the best fit containing `peak`.

        Parameters
        ----------
        peak : :class:`~.FittedPeak`
            The peak to start the search from
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        left_search_limit : int, optional
            The number of steps to search to the left of `peak`. Defaults to 3
        right_search_limit : int, optional
            The number of steps to search to the right of `peak`. Defaults to 3
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|

        Returns
        -------
        :class:`~.NetworkedTargetedDeconvolutionResult`
        """
        self._explore_local(
            peak, error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit,
            right_search_limit=right_search_limit, charge_carrier=charge_carrier,
            truncate_after=truncate_after, ignore_below=ignore_below)
        result = NetworkedTargetedDeconvolutionResult(self, peak)
        self._priority_map[peak] = result
        return result

    def _deconvolution_step(self, iteration_step, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),   # pylint: disable=arguments-differ
                            left_search_limit=1, right_search_limit=0, charge_carrier=PROTON,
                            truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW, **kwargs):
        if iteration_step != 0:
            self.peak_dependency_network.reset()
        self.populate_graph(
            error_tolerance=error_tolerance, charge_range=charge_range,
            left_search_limit=left_search_limit, right_search_limit=right_search_limit,
            charge_carrier=charge_carrier,
            truncate_after=truncate_after, ignore_below=ignore_below)
        self.postprocess_fits(
            charge_range=charge_range,
            charge_carrier=charge_carrier,
            error_tolerance=error_tolerance)
        self.select_best_disjoint_subgraphs(error_tolerance, charge_carrier)
        self._slice_cache.clear()

    def deconvolute(self, error_tolerance=ERROR_TOLERANCE, charge_range=(1, 8),   # pylint: disable=arguments-differ
                    left_search_limit=1, right_search_limit=0, iterations=MAX_ITERATION,
                    charge_carrier=PROTON, truncate_after=TRUNCATE_AFTER, ignore_below=IGNORE_BELOW,
                    convergence=CONVERGENCE, **kwargs):
        """Completely deconvolute the spectrum.

        For each iteration, clear :attr:`peak_depencency_network`, then invoke :meth:`populate_graph`
        followed by :meth:`select_best_disjoint_subgraphs` to populate the resulting
        :class:`~.DeconvolutedPeakSet`

        Parameters
        ----------
        error_tolerance : float, optional
            The parts-per-million error tolerance in m/z to search with. Defaults to |ERROR_TOLERANCE|
        charge_range : tuple, optional
            The range of charge states to consider. Defaults to (1, 8)
        order_chooser : callable, optional:
            A callable used as a key function for sorting peaks into the order they will
            be visited during deconvolution. Defaults to :obj:`operator.attrgetter("index")`
        left_search_limit : int, optional
            The number of steps to search to the left of :obj:`peak`. Defaults to 1
        right_search_limit : int, optional
            The number of steps to search to the right of :obj:`peak`. Defaults to 0
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to |PROTON|
        truncate_after : float, optional
            The percent of intensity to ensure is included in a theoretical isotopic pattern
            starting from the monoisotopic peak. This will cause theoretical isotopic patterns
            to be truncated, excluding trailing peaks which do not contribute substantially to
            the overall shape of the isotopic pattern.
        ignore_below : float, optional
            The minimum relative abundance to consider a peak in a theoretical isotopic
            pattern
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
        i = 0
        for i in range(iterations):
            if self.verbose:
                info("<== Starting Iteration %d ===================>", (i, ))
            self._deconvolution_step(
                i, error_tolerance=error_tolerance, charge_range=charge_range,
                left_search_limit=left_search_limit, right_search_limit=right_search_limit,
                charge_carrier=charge_carrier, truncate_after=truncate_after,
                ignore_below=ignore_below, **kwargs)
            end_signal = sum([p.intensity for p in self.peaklist]) + 1
            if (begin_signal - end_signal) / end_signal < convergence:
                if self.verbose:
                    info("(%0.4e - %0.4e) / %0.4e < %0.2g, Converged!",
                         begin_signal, end_signal, end_signal, convergence)
                break
            begin_signal = end_signal
        else:
            if self.verbose:
                info("Did Not Converge.")
        if self.verbose:
            info("Finished Deconvolution in %d Iterations", (i + 1, ))
        if self.merge_isobaric_peaks:
            self._deconvoluted_peaks = self._merge_peaks(
                self._deconvoluted_peaks)

        return DeconvolutedPeakSet(list(self._deconvoluted_peaks)).reindex()


try:
    _has_c = True
    from ms_deisotope._c.deconvoluter_base import (
        _explore_local as _c_explore_local,
        populate_graph as cpopulate_graph)
    PeakDependenceGraphDeconvoluterBase._explore_local = _c_explore_local
    PeakDependenceGraphDeconvoluterBase.populate_graph = cpopulate_graph
except ImportError as e:
    _has_c = False
