import operator
from collections import defaultdict

from .subgraph import ConnectedSubgraph
from .intervals import SpanningMixin
from ..utils import Base


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
        return "PeakNode(%s, %d)" % (self.peak, self.links)


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
        self.dependencies.sort(key=lambda x: x.score)
        self._reset()

    def disjoint_subset(self):
        graph = ConnectedSubgraph(self.dependencies, maximize=self.maximize)
        return graph.find_heaviest_path()

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

    def interval_contains_point(self, point):
        return self.start <= point < self.end

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


class PeakDependenceGraph(object):
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
        self._clusters = None

    def _populate_initial_graph(self):
        for peak in self.peaklist:
            self.nodes[peak.index] = PeakNode(peak)

    def add_fit_dependence(self, fit_record):
        for peak in fit_record.experimental:
            if peak.index == 0:
                continue
            self.nodes[peak.index].links[fit_record] = fit_record.score
        self.dependencies.add(fit_record)

    def nodes_for(self, fit_record):
        return [self.nodes[p.index] for p in fit_record.experimental if p.index != 0]

    def drop_fit_dependence(self, fit_record):
        for node in self.nodes_for(fit_record):
            del node.links[fit_record]

    def claimed_nodes(self):
        peaks = set()
        for fit in self.dependencies:
            for peak in fit.experimental:
                if peak.index == 0:
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
            if mono_peak.index == 0:
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

        for node in self.nodes.values():
            # This peak is depended upon by each fit in `dependencies`
            dependencies = set(node.links.keys())

            if len(dependencies) == 0:
                continue

            # These fits also depend upon these other peaks, and those peaks are depended upon
            # for other fits in turn, which depend upon the assignment of this peak.
            for dep in list(dependencies):
                for node in self.nodes_for(dep):
                    dependencies |= clusters[node]

            # Create a fresh copy to share out again to avoid eliminate possible errors of shared storage.
            # This is likely unecessary under the current implementation.
            dependencies = set(dependencies)

            # Update all co-depended nodes with the full set of all fits which depend upon them
            for dep in dependencies:
                for node in self.nodes_for(dep):
                    clusters[node] = dependencies

        # Use an `id` keyed dictionary over these copied sets to ensure we have exactly one reference
        # to each set of inter-dependent fits, and then convert each set into an instance of `DependenceCluster`
        clusters = [DependenceCluster(dependencies=c, maximize=self.maximize) for c in {
            id(v): v for v in clusters.values() if v}.values()]
        clusters = sorted(clusters, key=operator.attrgetter("start"))
        self._clusters = clusters
        return clusters

    def __repr__(self):
        return "PeakDependenceNetwork(%s, %d)" % (self.peaklist, len(self.dependencies))
