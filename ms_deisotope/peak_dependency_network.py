from collections import defaultdict

from .utils import Base
from .scoring import IsotopicFitRecord


class PeakNode(Base):
    def __init__(self, peak, links=None):
        if links is None:
            links = {}
        self.peak = peak
        self.links = links

    def __hash__(self):
        return hash(self.peak)

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
        if dependencies is None:
            dependencies = set()
        else:
            dependencies = set(dependencies)
        self.parent = parent
        self.dependencies = dependencies

    def add(self, fit):
        self.dependencies.add(fit)

    def __contains__(self, fit):
        return fit in self.dependencies

    def __eq__(self, other):
        return self.dependencies == other.dependencies


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
            self.nodes[peak.index].links[fit_record] = fit_record.score
        self.dependencies.add(fit_record)

    def claimed_nodes(self):
        peaks = set()
        for fit in self.dependencies:
            for peak in fit.experimental:
                peaks.add(self.nodes[peak.index])
        return peaks

    def find_highest_scoring_active_set(self):
        active_set = set()

        claimed_nodes = defaultdict(DependenceCluster)
        for fit_record in self.dependencies:
            for peak in fit_record.experimental:
                selected_cluster = claimed_nodes[peak.index]
                selected_cluster.add(fit_record)
            


