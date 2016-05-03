from collections import defaultdict

from .utils import Base
from .scoring import IsotopicFitRecord


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
        self._hash = None

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.peak)
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


def isdisjoint_list(a, b):
    for pa, pb in zip(a, b):
        if pa == pb:
            return False
    return True


class DependenceCluster(Base):
    """
    Represent a cluster of peak fits which are overlapping

    Attributes
    ----------
    dependencies : set
        The set of all fits which are codependent
    """
    def __init__(self, parent=None, dependencies=None):
        if parent is None:
            parent = self
        if dependencies is None:
            dependencies = []
        else:
            dependencies = sorted(dependencies, key=lambda x: x.score, reverse=True)
        self.parent = parent
        self.dependencies = dependencies
        self.rank = 0

    def add(self, fit):
        self.dependencies.append(fit)
        self.dependencies.sort(key=lambda x: x.score)

    @property
    def best_fit(self):
        return self.dependencies[0]

    @property
    def start(self):
        return min([member.experimental[0].mz for member in self.dependencies])

    @property
    def end(self):
        return max([member.experimental[-1].mz for member in self.dependencies])

    def __contains__(self, fit):
        return fit in self.dependencies

    def __eq__(self, other):
        return self.dependencies == other.dependencies

    def __repr__(self):
        return "DependenceCluster(dependencies=[%s], start=%0.4f, end=%0.4f)" % (
            '\n'.join(map(str, self.dependencies[:10])), self.start, self.end)


class PeakDependenceGraph(object):
    def __init__(self, peaklist, nodes=None, dependencies=None):
        if nodes is None:
            nodes = {}
        if dependencies is None:
            dependencies = set()
        self.peaklist = peaklist
        self.nodes = nodes
        self.dependencies = dependencies
        if len(self.nodes) == 0:
            self._populate_initial_graph()

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

    def claimed_nodes(self):
        peaks = set()
        for fit in self.dependencies:
            for peak in fit.experimental:
                if peak.index == 0:
                    continue
                peaks.add(self.nodes[peak.index])
        return peaks

    def find_non_overlapping_intervals(self, max_size=20):
        clusters = defaultdict(set)

        for node in self.nodes.values():
            dependencies = set(node.links.keys())

            if len(dependencies) == 0:
                continue

            for dep in list(dependencies):
                for node in self.nodes_for(dep):
                    dependencies |= clusters[node]

            dependencies = sorted(dependencies, key=lambda x: x.score, reverse=True)
            dependencies = set(dependencies[:max_size])

            for dep in dependencies:
                for node in self.nodes_for(dep):
                    clusters[node] = dependencies

        clusters = [DependenceCluster(dependencies=c) for c in {
            id(v): v for v in clusters.values() if v}.values()]
        return clusters
