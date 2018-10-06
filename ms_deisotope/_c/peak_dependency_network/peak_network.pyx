import operator
import warnings

from collections import defaultdict

cimport cython

from ms_peak_picker._c.peak_set cimport FittedPeak, PeakSet

from ms_deisotope._c.scoring cimport IsotopicFitRecord
from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin, IntervalTreeNode
from ms_deisotope._c.peak_dependency_network.subgraph cimport ConnectedSubgraph


@cython.freelist(1000000)
cdef class PeakNode(object):
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

    cdef:
        public FittedPeak peak
        public dict links
        long _hash

    @staticmethod
    cdef PeakNode _create(FittedPeak peak):
        cdef PeakNode inst = PeakNode.__new__(PeakNode)
        inst.peak = peak
        inst.links = dict()
        inst._hash = hash(peak)
        return inst

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


cdef class DependenceCluster(SpanningMixin):
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

    cdef:
        public object parent
        public list dependencies
        public bint maximize
        public IsotopicFitRecord best_fit

    def __init__(self, parent=None, dependencies=None, maximize=True):
        if parent is None:
            parent = self
        if dependencies is None:
            dependencies = []
        else:
            dependencies = sorted(dependencies, key=lambda x: x.score, reverse=maximize)
        self.parent = parent
        self.dependencies = dependencies
        self.maximize = maximize
        self._reset()

    cpdef _reset(self):
        self.start = self._start()
        self.end = self._end()
        self.best_fit = self._best_fit()

    cpdef add(self, IsotopicFitRecord fit):
        """
        Adds a new IsotopicFitRecord to this cluster, and ensures the sorted
        property still holds

        Parameters
        ----------
        fit : IsotopicFitRecord
            The fit being added to the record

        """
        self.dependencies.append(fit)
        self.dependencies.sort(reverse=self.maximize)
        self._reset()

    cpdef list disjoint_subset(self):
        graph = self.build_graph()
        return graph.find_heaviest_path()

    cpdef ConnectedSubgraph build_graph(self):
        graph = ConnectedSubgraph(self.dependencies, maximize=self.maximize)
        return graph

    cdef IsotopicFitRecord _best_fit(self):
        """
        Retrieve the absolute best single fit for the cluster

        Returns
        -------
        name : IsotopicFitRecord
            The best scoring IsotopicFitRecord in this cluster
        """
        return self.dependencies[0]

    cpdef list disjoint_best_fits(self):
        """
        Compute the best set of disjoint isotopic fits spanning this cluster

        Returns
        -------
        list of IsotopicFitRecord
        """
        fit_sets = tuple(self.disjoint_subset())
        best_fits = fit_sets
        return [node.fit for node in best_fits]

    cdef double _start(self):
        """
        Determine the first mz coordinate for members of this cluster

        Returns
        -------
        float
        """
        return min([member.experimental[0].mz for member in self.dependencies])

    cdef double _end(self):
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


cdef class PeakDependenceGraph(object):
    cdef:
        public dict nodes
        public set dependencies
        public PeakSet peaklist
        public int max_missed_peaks
        public bint use_monoisotopic_superceded_filtering
        public bint maximize

        public list clusters
        public dict _solution_map
        public list _all_clusters

        public IntervalTreeNode _interval_tree

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

    def reset(self):
        # Keep a record of all clusters from previous iterations
        self._all_clusters.extend(
            self.clusters if self.clusters is not None else [])
        self.nodes = dict()
        self.dependencies = set()
        self._interval_tree = None
        self._populate_initial_graph()

    def add_solution(self, key, solution):
        self._solution_map[key] = solution

    @property
    def interval_tree(self):
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

    def _deep_fuzzy_solution_for(self, peak, shift=0.5):
        pass

    def _find_fuzzy_solution_for(self, peak, shift=0.5):
        tree = self.interval_tree
        clusters = tree.contains_point(peak.mz + shift)
        if len(clusters) == 0:
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
            return self._solution_map[fit]

    def find_solution_for(self, peak):
        peak_node = self.nodes[peak.index]
        tree = self.interval_tree
        clusters = tree.contains_point(peak.mz)
        if len(clusters) == 0:
            return self._find_fuzzy_solution_for(peak)
        best_fits = [cluster.disjoint_best_fits() for cluster in clusters]

        acc = []
        for fits in best_fits:
            acc.extend(fits)
        best_fits = acc

        # Extract only fits that use the query peak
        common = tuple(set(best_fits) & set(peak_node.links))

        if len(common) > 1 or len(common) == 0:
            if len(common) > 1:
                warnings.warn("Too many solutions exist for %r" % peak)
            # If there were no fits for this peak, then it may be that this peak
            # was not included in a fit. Try to find the nearest solution.
            i = 0
            err = float('inf')
            for j, case in enumerate(best_fits):
                case_err = abs(case.monoisotopic_peak.mz - peak.mz)
                if self.maximize:
                    case_err /= case.score
                else:
                    case_err /= (1. / case.score)
                if case_err < err:
                    i = j
                    err = case_err
            fit = best_fits[i]
        else:
            fit = common[0]
        return self._solution_map[fit]

    cpdef _populate_initial_graph(self):
        cdef:
            size_t i, n
            FittedPeak peak
        n = self.peaklist.get_size()
        for i in range(n):
            peak = self.peaklist.getitem(i)

        for peak in self.peaklist:
            self.nodes[peak.index] = PeakNode(peak)

    def add_fit_dependence(self, fit_record):
        for peak in fit_record.experimental:
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
        for node in self.nodes_for(fit_record):
            try:
                del node.links[fit_record]
            except KeyError:
                pass

    def claimed_nodes(self):
        peaks = set()
        for fit in self.dependencies:
            for peak in fit.experimental:
                if peak.peak_count < 0:
                    continue
                peaks.add(self.nodes[peak.index])
        return peaks

    def drop_superceded_fits(self):
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

    def best_exact_fits(self):
        by_peaks = defaultdict(list)
        best_fits = []
        for fit in self.dependencies:
            by_peaks[tuple(fit.experimental)].append(fit)
        for peak_tuple, fits in by_peaks.items():
            fits = sorted(fits, key=lambda x: x.score, reverse=not self.maximize)
            for fit in fits[:-1]:
                self.drop_fit_dependence(fit)
            best_fits.append(fits[-1])
        self.dependencies = set(best_fits)

    def drop_gapped_fits(self, n=None):
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
        clusters = defaultdict(set)
        self.drop_gapped_fits()
        self.best_exact_fits()
        if self.use_monoisotopic_superceded_filtering:
            self.drop_superceded_fits()

        nodes_for_cache = {}

        for node in self.nodes.values():
            # This peak is depended upon by each fit in `dependencies`
            dependencies = set(node.links.keys())

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
        return "PeakDependenceNetwork(%s, %d)" % (self.peaklist, len(self.dependencies))


class NoIsotopicClustersError(ValueError):
    def __init__(self, msg, graph=None, scan_id=None):
        ValueError.__init__(self, msg)
        self.graph = graph
        self.scan_id = scan_id
