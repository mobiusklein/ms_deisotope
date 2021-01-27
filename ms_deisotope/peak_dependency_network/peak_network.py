import operator
import logging

from collections import defaultdict

from .subgraph import ConnectedSubgraph
from .intervals import SpanningMixin, IntervalTreeNode
from ..utils import Base, TargetedDeconvolutionResultBase
from ..task import LogUtilsMixin

logger = logging.getLogger("deconvolution_scan_processor")


def ident(x):
    return x


class PeakNode(Base):
    """
    Represent a single FittedPeak and the junction of multiple
    IsotopicFitRecords which depend uponit

    Attributes
    ----------
    links : dict
        Mapping IsotopicFitRecords which depend upon this peak to those fit's scores
    peak : ms_peak_picker.FittedPeak
        The peak being depended upon
    """

    def __init__(self, peak, links=None):
        if links is None:
            links = {}
        self.peak = peak
        self.links = links
        self._hash = hash(self.peak)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        try:
            return self.peak == other.peak
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __contains__(self, fit):
        return fit in self.links

    def __repr__(self):
        return "PeakNode(%s, %s)" % (self.peak, self.links)


class DependenceCluster(SpanningMixin):
    """
    Represent a cluster of isotopic fits which are overlapping

    Attributes
    ----------
    dependencies : set
        The set of all fits which are codependent
    start: float
        The smallest mz of all elements of this cluster
    end: float
        The largest mz of all elements of this cluster
    maximize: bool
        Whether the objective is to maximize or minimize the
        isotopic fit score
    """

    def __init__(self, parent=None, dependencies=None, maximize=True):
        if parent is None:
            parent = self
        if dependencies is None:
            dependencies = []
        else:
            dependencies = sorted(dependencies, key=lambda x: x.score, reverse=maximize)
        self.parent = parent
        self.dependencies = dependencies
        self.rank = 0
        self.maximize = maximize
        self._reset()

    def _reset(self):
        self.start = self._start()
        self.end = self._end()
        self.best_fit = self._best_fit()

    def add(self, fit):
        """
        Adds a new IsotopicFitRecord to this cluster, and ensures the sorted
        property still holds

        Parameters
        ----------
        fit : IsotopicFitRecord
            The fit being added to the record

        """
        self.dependencies.append(fit)
        self.dependencies.sort(key=lambda x: x.score, reverse=self.maximize)
        self._reset()

    def disjoint_subset(self):
        graph = self.build_graph()
        return graph.find_heaviest_path()

    def build_graph(self):
        graph = ConnectedSubgraph(self.dependencies, maximize=self.maximize)
        return graph

    def _best_fit(self):
        """
        Retrieve the absolute best single fit for the cluster

        Returns
        -------
        name : IsotopicFitRecord
            The best scoring IsotopicFitRecord in this cluster
        """
        return self.dependencies[0]

    def disjoint_best_fits(self):
        """
        Compute the best set of disjoint isotopic fits spanning this cluster

        Returns
        -------
        list of IsotopicFitRecord
        """
        fit_sets = tuple(self.disjoint_subset())
        best_fits = fit_sets
        return [node.fit for node in best_fits]

    def _start(self):
        """
        Determine the first mz coordinate for members of this cluster

        Returns
        -------
        float
        """
        return min([member.experimental[0].mz for member in self.dependencies])

    def _end(self):
        """
        Determines the last mz coordinate for members of this cluster

        Returns
        -------
        float
        """
        return max([member.experimental[-1].mz for member in self.dependencies])

    def __len__(self):
        return len(self.dependencies)

    def __iter__(self):
        return iter(self.dependencies)

    def __contains__(self, fit):
        return fit in self.dependencies

    def __eq__(self, other):
        return self.dependencies == other.dependencies

    def __repr__(self):
        return "DependenceCluster(dependencies=[%s], start=%0.4f, end=%0.4f)" % (
            '\n'.join(map(str, self.dependencies[:10])), self.start, self.end)

    def __getitem__(self, i):
        return self.dependencies[i]

    def uses_mz(self, mz):
        for fit in self:
            for peak in fit.experimental:
                if peak.mz == mz:
                    return True
        return False

    def fits_using_mz(self, mz):
        fits = []
        for fit in self:
            for peak in fit.experimental:
                if peak.mz == mz:
                    fits.append(fit)
                    break
        return fits

try:
    _has_c = True
    from ms_deisotope._c.peak_dependency_network.peak_network import (
        PeakNode, DependenceCluster, PeakDependenceGraphBase)
except ImportError:
    _has_c = False

    class PeakDependenceGraphBase(object):

        def reset(self):
            # Keep a record of all clusters from previous iterations
            self._all_clusters.extend(
                self.clusters if self.clusters is not None else [])
            self.nodes = dict()
            self.dependencies = set()
            self._interval_tree = None
            self._populate_initial_graph()

        def _populate_initial_graph(self):
            for peak in self.peaklist:
                self.nodes[peak.index] = PeakNode(peak)

        def add_fit_dependence(self, fit_record):
            '''Add the relatoinship between the experimental peaks
            in `fit_record` to the graph, expressed as a hyper-edge
            denoted by `fit_record`.

            This adds `fit_record` to :attr:`PeakNode.links` for each
            node corresponding to the :class:`~.FittedPeak` instance in
            :attr:`IsotopicFitRecord.experimental` of `fit_record`. It also
            adds `fit_record to :attr:`dependencies`

            Parameters
            ----------
            fit_record: :class:`~.IsotopicFitRecord`
            '''
            for peak in fit_record.experimental:
                # Check for null peak
                if peak.index == 0:
                    continue
                self.nodes[peak.index].links[fit_record] = fit_record.score
            self.dependencies.add(fit_record)

        def nodes_for(self, fit_record, cache=None):
            if cache is None:
                return [self.nodes[p.index] for p in fit_record.experimental if p.peak_count >= 0]
            else:
                try:
                    return cache[fit_record]
                except KeyError:
                    cache[fit_record] = value = [
                        self.nodes[p.index] for p in fit_record.experimental if p.peak_count >= 0]
                    return value

        def drop_fit_dependence(self, fit_record):
            '''Remove this fit from the graph, deleting all
            hyper-edges.
            '''
            for node in self.nodes_for(fit_record):
                try:
                    del node.links[fit_record]
                except KeyError:
                    pass

        def best_exact_fits(self):
            '''For each distinct group of experimental peaks, retain only
            the best scoring fit using exactly those peaks.
            '''
            by_peaks = defaultdict(list)
            best_fits = []
            for fit in self.dependencies:
                by_peaks[tuple(fit.experimental)].append(fit)
            for _peak_tuple, fits in by_peaks.items():
                fits = sorted(fits, key=lambda x: x.score,
                              reverse=not self.maximize)
                for fit in fits[:-1]:
                    self.drop_fit_dependence(fit)
                best_fits.append(fits[-1])
            self.dependencies = set(best_fits)

        def _gather_independent_clusters(self, nodes_for_cache=None):
            clusters = defaultdict(set)
            if nodes_for_cache is None:
                nodes_for_cache = dict()
            seed_nodes = self.nodes.values()
            for seed_node in seed_nodes:
                # This peak is depended upon by each fit in `dependencies`
                dependencies = set(seed_node.links.keys())

                if len(dependencies) == 0:
                    continue

                # These fits also depend upon these other peaks, and those peaks are depended upon
                # for other fits in turn, which depend upon the assignment of this peak.
                for dep in list(dependencies):
                    for node in self.nodes_for(dep, nodes_for_cache):
                        dependencies |= clusters[node]

                # Create a fresh copy to share out again to avoid eliminate possible errors of shared storage.
                # This is likely unecessary under the current implementation.
                dependencies = set(dependencies)

                # Update all co-depended nodes with the full set of all fits which depend upon them
                for dep in dependencies:
                    for node in self.nodes_for(dep, nodes_for_cache):
                        clusters[node] = dependencies
            return clusters


def dependency_cluster_uses_mz(self, mz):
    for fit in self:
        for peak in fit.experimental:
            if peak.mz == mz:
                return True
    return False


def dependency_cluster_fits_using_mz(self, mz):
    fits = []
    for fit in self:
        for peak in fit.experimental:
            if peak.mz == mz:
                fits.append(fit)
                break
    return fits


class PeakDependenceGraph(PeakDependenceGraphBase, LogUtilsMixin):
    def __init__(self, peaklist, nodes=None, dependencies=None, max_missed_peaks=1,
                 use_monoisotopic_superceded_filtering=True, maximize=True):
        if nodes is None:
            nodes = {}
        if dependencies is None:
            dependencies = set()
        self.peaklist = peaklist
        self.nodes = nodes
        self.dependencies = dependencies
        self.max_missed_peaks = max_missed_peaks
        self.use_monoisotopic_superceded_filtering = use_monoisotopic_superceded_filtering
        self.maximize = maximize
        if len(self.nodes) == 0:
            self._populate_initial_graph()
        self.clusters = None
        self._interval_tree = None
        self._solution_map = {}
        self._all_clusters = []

    def add_solution(self, key, solution):
        self._solution_map[key] = solution

    @property
    def interval_tree(self):
        '''An :class:`~.IntervalTreeNode` built over all dependency clusters that
        the :class:`PeakDependenceGraph` has seen over the course of all iterations.

        Raises
        ------
        NoIsotopicClustersError:
            When no isotopic clusters, and in turn no :class:`DependenceCluster`s are found

        Returns
        -------
        :class:`~.IntervalTreeNode`
        '''
        if self._interval_tree is None:
            # Build tree from both the current set of clusters
            # and all clusters from previous iterations
            self._interval_tree = IntervalTreeNode.build(
                self.clusters + self._all_clusters)
            if self._interval_tree is None:
                raise NoIsotopicClustersError(
                    "Could not build intervals for peak retrieval with %d clusters" % len(
                        self.clusters), self)
        return self._interval_tree

    def _find_fuzzy_solution_for(self, peak, shift=0.5):
        tree = self.interval_tree
        clusters = tree.contains_point(peak.mz + shift)
        if len(clusters) == 0:
            # No dependency clusters span this point. No signal?
            return None
        else:
            best_fits = [cluster.disjoint_best_fits() for cluster in clusters]

            acc = []
            for fits in best_fits:
                acc.extend(fits)
            best_fits = acc

            index = 0
            error = float('inf')

            for i, fit in enumerate(best_fits):
                err = abs(fit.monoisotopic_peak.mz - peak.mz)
                if err < error:
                    error = err
                    index = i
            fit = best_fits[index]
            try:
                return self._solution_map[fit]
            except KeyError:
                best_fits.sort(key=lambda f: abs(f.monoisotopic_peak.mz - peak.mz))
                # The first index is the minimum error, which the loop above found,
                # but by virtue of the KeyError getting us here, we know it is not
                # present in the solution map.
                for f in  best_fits[1:]:
                    try:
                        return self._solution_map[f]
                    except KeyError:
                        continue
                # No solution could be found
                return None

    def find_cluster_for(self, peak):
        tree = self.interval_tree
        if tree is None:
            raise TypeError(
                "This network hasn't constructed clusters yet. Call find_non_overlapping_intervals first!")
        clusters = tree.contains_point(peak.mz)
        for cluster in clusters:
            # by construction, only one cluster contains a particular peak
            if dependency_cluster_uses_mz(cluster, peak.mz):
                return cluster

    def find_solution_for(self, peak):
        '''Find the best isotopic pattern fit which includes ``peak``
        or is close enough to ``peak`` that it might have been selected
        instead.

        Parameters
        ----------
        peak: :class:`~.FittedPeak`
            The peak to find the best fit for.

        Returns
        -------
        :class:`~.DeconvolutedPeak`
        '''
        peak_node = self.nodes[peak.index]
        tree = self.interval_tree
        clusters = tree.contains_point(peak.mz)
        if not clusters:
            return self._find_fuzzy_solution_for(peak)
        best_fits = [cluster.disjoint_best_fits() for cluster in clusters]

        acc = []
        for fits in best_fits:
            acc.extend(fits)
        best_fits = acc

        self.debug("For query m/z %f has candidates %r" % (peak.mz, best_fits))

        # Extract only fits that use the query peak
        common = tuple(set(best_fits) & set(peak_node.links))

        used_common = False
        if len(common) > 1 or not common:
            if len(common) > 1:
                self.debug("There is not exactly one solution for %r (%d)" % (peak, len(common)))
            # If there were no fits for this peak, then it may be that this peak
            # was not included in a fit. Try to find the nearest solution.
            fit = self._find_best_fit_by_weighted_distance(peak, best_fits)
        else:
            used_common = True
            fit = common[0]
        try:
            return self._solution_map[fit]
        except KeyError:
            self.debug(
                "The closest common fit was not chosen. The closest fit was common? %r" % (
                    used_common))
            candidates = list(set(best_fits) & set(self._solution_map))
            fit = self._find_best_fit_by_weighted_distance(peak, candidates)
            return self._solution_map[fit]

    def _find_best_fit_by_weighted_distance(self, peak, candidates):
        '''Locate the fit that is nearest to ``peak`` in m/z space, weighted by
        the score of the fit such that given two fits with the same m/z, the better
        scoring fit is chosen.

        Parameters
        ----------
        peak: :class:`~.FittedPeak`
            The peak to find solutions with reference to.
        candidates: :class:`list`
            A set of candidate :class:`~.IsotopicFitRecord` to compare between

        Returns
        -------
        :class:`~.IsotopicFitRecord`
        '''
        i = 0
        err = float('inf')
        for j, case in enumerate(candidates):
            case_err = abs(case.monoisotopic_peak.mz - peak.mz)
            if self.maximize:
                case_err /= case.score
            else:
                case_err /= (1. / (case.score + 1e-6))
            if case_err < err:
                i = j
                err = case_err
        fit = candidates[i]
        return fit

    def claimed_nodes(self):
        '''Return all nodes with at least one
        fit connecting it.

        Returns
        -------
        :class:`set`
        '''
        peaks = set()
        for fit in self.dependencies:
            for peak in fit.experimental:
                if peak.peak_count < 0:
                    continue
                peaks.add(self.nodes[peak.index])
        return peaks

    def drop_superceded_fits(self):
        '''Drop all fits where there is a better fit
        using the first fit's monoisotopic peak at the
        same charge state.
        '''
        suppressed = []
        keep = []

        if self.maximize:
            comparator = operator.lt
        else:
            comparator = operator.gt

        for dep in self.dependencies:
            mono_peak = dep.monoisotopic_peak
            if mono_peak.peak_count < 0:
                continue
            mono_peak_node = self.nodes[mono_peak.index]
            suppress = False
            for candidate in mono_peak_node.links:
                if dep.charge == candidate.charge:
                    if mono_peak in candidate.experimental:
                        if comparator(dep.score, candidate.score):
                            suppress = True
                            break
            if suppress:
                suppressed.append(dep)
            else:
                keep.append(dep)
        for dropped in suppressed:
            self.drop_fit_dependence(dropped)
        self.dependencies = set(keep)

    def drop_gapped_fits(self, n=None):
        '''Discard any fit with more missed peaks than
        :attr:`max_missed_peaks`.
        '''
        if n is None:
            n = self.max_missed_peaks
        keep = []
        for dep in self.dependencies:
            if dep.missed_peaks > n:
                self.drop_fit_dependence(dep)
            else:
                keep.append(dep)
        self.dependencies = set(keep)

    def find_non_overlapping_intervals(self):
        self.drop_gapped_fits()
        self.best_exact_fits()
        if self.use_monoisotopic_superceded_filtering:
            self.drop_superceded_fits()
        clusters = self._gather_independent_clusters()
        # Use an `id` keyed dictionary over these copied sets to ensure we have exactly one reference
        # to each set of inter-dependent fits, and then convert each set into an instance of `DependenceCluster`
        clusters = [DependenceCluster(dependencies=c, maximize=self.maximize) for c in {
            id(v): v for v in clusters.values() if v}.values()]
        clusters = sorted(clusters, key=operator.attrgetter("start"))
        self.clusters = clusters
        return clusters

    def __iter__(self):
        return iter(self.clusters)

    def __getitem__(self, i):
        return self.clusters[i]

    def __repr__(self):
        return "PeakDependenceGraph(%s, %d)" % (self.peaklist, len(self.dependencies))


PeakDependenceGraph.log_with_logger(logger)


class NetworkedTargetedDeconvolutionResult(TargetedDeconvolutionResultBase):
    '''Stores the necessary information to retrieve the local optimal solution for
    a single peak for a deconvolution algorithm from the optimal fit containing
    :attr:`query_peak` in the set of disjoint best fits for the enclosing connected
    component

    Attributes
    ----------
    query_peak : FittedPeak
        The peak queried with
    solution_peak : DeconvolutedPeak
        The optimal solution peak
    '''
    def __init__(self, deconvoluter, peak, *args, **kwargs):
        super(NetworkedTargetedDeconvolutionResult, self).__init__(deconvoluter, *args, **kwargs)
        self.query_peak = peak
        self.solution_peak = None

    def _get_solution(self):
        try:
            self.solution_peak = self.deconvoluter.peak_dependency_network.find_solution_for(
                self.query_peak)
        except IndexError:
            pass

    def get(self):
        """Fetch the optimal solution after the computation has finished.

        Returns
        -------
        DeconvolutedPeak
        """
        self._get_solution()
        return self.solution_peak


class NoIsotopicClustersError(ValueError):
    def __init__(self, msg, graph=None, scan_id=None):
        ValueError.__init__(self, msg)
        self.graph = graph
        self.scan_id = scan_id
