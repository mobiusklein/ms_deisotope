import warnings
import operator
from collections import defaultdict

from ms_deisotope.peak_dependency_network.intervals import (
    # Use pure Python implementations to allow non-primitive floats as coordinates
    _SpanningMixin as SpanningMixin, _IntervalTreeNode as IntervalTreeNode)
from ms_deisotope.peak_dependency_network.subgraph import GreedySubgraphSelection

from .lcms_feature import EmptyFeature
from .feature_fit import map_coord


def is_valid(feature):
    return feature is not None and not isinstance(feature, EmptyFeature)


class FeatureNode(object):
    """Holds all the information about a single `LCMSFeature` instance
    and all the `LCMSFeatureFit` instances that depend upon it.

    Attributes
    ----------
    feature : LCMSFeature
    links : dict
        Mapping from feature fit to score
    """
    def __init__(self, feature, links=None):
        if links is None:
            links = {}
        self.feature = feature
        self.links = links
        self._hash = hash(self.feature)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        try:
            return self.feature == other.feature
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __contains__(self, fit):
        return fit in self.links

    def __repr__(self):
        return "FeatureNode(%s, %s)" % (self.feature, self.links)

try:
    _FeatureNode = FeatureNode
    has_c = True
    from ms_deisotope._c.feature_map.dependence_network import FeatureNode
except ImportError:
    has_c = False

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
        self.best_fit = None
        self.maximize = maximize
        self._reset()

    def _reset(self):
        self.start = self._start()
        self.end = self._end()
        self.best_fit = self._best_fit()

    def add(self, fit):
        """
        Adds a new LCMSFeatureSetFit to this cluster, and ensures the sorted
        property still holds

        Parameters
        ----------
        fit : LCMSFeatureSetFit
            The fit being added to the record

        """
        self.dependencies.append(fit)
        self.dependencies.sort(key=lambda x: x.score, reverse=self.maximize)
        self._reset()

    def disjoint_subset(self, max_size=1000):
        if len(self.dependencies) > max_size:
            dependencies = self.dependencies[:max_size]
        else:
            dependencies = self.dependencies
        graph = ConnectedSubgraph(dependencies, maximize=self.maximize)
        return graph.find_heaviest_path()

    def _best_fit(self):
        """
        Retrieve the absolute best single fit for the cluster

        Returns
        -------
        name : LCMSFeatureSetFit
            The best scoring LCMSFeatureSetFit in this cluster
        """
        return self.dependencies[0]

    def disjoint_best_fits(self, max_size=1000):
        """
        Compute the best set of disjoint isotopic fits spanning this cluster

        Returns
        -------
        list of LCMSFeatureSetFit
        """
        fit_sets = tuple(self.disjoint_subset(max_size=max_size))
        best_fits = fit_sets
        return [node.fit for node in best_fits]

    def _start(self):
        """
        Determine the first mz coordinate for members of this cluster

        Returns
        -------
        float
        """
        return map_coord(*min([member.start for member in self.dependencies]))

    def _end(self):
        """
        Determines the last mz coordinate for members of this cluster

        Returns
        -------
        float
        """
        return map_coord(*max([member.end for member in self.dependencies]))

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
        return "DependenceCluster(dependencies=[%s], start=%r, end=%r)" % (
            '\n'.join(map(str, self.dependencies[:10])), self.start, self.end)

    def __getitem__(self, i):
        return self.dependencies[i]


try:
    _DependenceCluster = DependenceCluster
    has_c = True
    from ms_deisotope._c.feature_map.dependence_network import DependenceCluster
except ImportError:
    has_c = False


class FeatureSetFitNode(SpanningMixin):
    def __init__(self, fit, index=None):
        self.fit = fit
        self.edges = set()
        self.overlap_edges = set()

        self._hash = None
        self.score = fit.score

        self.feature_indices = {(f.mz, f.start_time, f.end_time) for f in fit.features if is_valid(f)}

        self.start = (fit.start[0], max(f.start_time for f in fit.features if is_valid(f)))
        self.end = (fit.end[0], min(f.end_time for f in fit.features if is_valid(f)))

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.fit)
        return self._hash

    def __eq__(self, other):
        return self.fit is other.fit

    def __gt__(self, other):
        return self.fit > other.fit

    def __lt__(self, other):
        return self.fit < other.fit

    def __ne__(self, other):
        return self.fit is not other.fit

    def visit(self, other):
        if self.isdisjoint(other):
            self.edges.add(other)
            other.edges.add(self)
        else:
            self.overlap_edges.add(other)
            other.overlap_edges.add(self)

    def isdisjoint(self, other):
        return self.feature_indices.isdisjoint(other.feature_indices)

    def __repr__(self):
        return "FeatureSetFitNode(%r)" % self.fit

    @staticmethod
    def overlap(a, b):
        return len(a.feature_indices & b.feature_indices) > 0


try:
    _FeatureSetFitNode = FeatureSetFitNode
    has_c = True
    from ms_deisotope._c.feature_map.dependence_network import FeatureSetFitNode
except ImportError:
    has_c = False


class ConnectedSubgraph(object):
    def __init__(self, fits, maximize=True):
        self.nodes = tuple(map(FeatureSetFitNode, fits))
        self.maximize = maximize
        i = 0
        for node in self.nodes:
            node.index = i
            i += 1
        self.populate_edges()

    def __getitem__(self, i):
        return self.nodes[i]

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def populate_edges(self):
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            for j in range(i + 1, len(self.nodes)):
                other = self.nodes[j]
                node.visit(other)

    def find_heaviest_path(self, method="greedy"):
        if len(self) == 1:
            return set(self.nodes)
        if method == "greedy":
            solution = GreedySubgraphSelection.solve(
                tuple(self), maximize=self.maximize, overlap_fn=FeatureSetFitNode.overlap)
            return solution
        else:
            raise NotImplementedError(method)


try:
    _ConnectedSubgraph = ConnectedSubgraph
    has_c = True
    from ms_deisotope._c.feature_map.dependence_network import ConnectedSubgraph
except ImportError:
    has_c = False


class FeatureDependenceGraphBase(object):

    def nodes_for(self, fit_record, cache=None):
        if cache is None:
            return [self.nodes[p] for p in fit_record.features if is_valid(p)]
        else:
            try:
                return cache[fit_record]
            except KeyError:
                cache[fit_record] = value = [
                    self.nodes[p] for p in fit_record.features if is_valid(p)]
                return value

    def add_fit_dependence(self, fit_record):
        for feature in fit_record.features:
            if not is_valid(feature):
                continue
            self.nodes[feature].links[fit_record] = fit_record.score
        self.dependencies.add(fit_record)

    def drop_superceded_fits(self):
        suppressed = []
        keep = []

        if self.maximize:
            comparator = operator.lt
        else:
            comparator = operator.gt

        for dep in self.dependencies:
            monoisotopic_feature = dep.monoisotopic_feature
            if not is_valid(monoisotopic_feature):
                continue
            monoisotopic_feature_node = self.nodes[monoisotopic_feature]
            suppress = False
            for candidate in monoisotopic_feature_node.links:
                if dep.charge == candidate.charge:
                    if monoisotopic_feature in candidate.features:
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

    def drop_fit_dependence(self, fit_record, cache=None):
        for node in self.nodes_for(fit_record, cache=cache):
            try:
                del node.links[fit_record]
            except KeyError:
                pass

    def drop_gapped_fits(self, n=None):
        if n is None:
            n = self.max_missing_features
        keep = []
        for dep in self.dependencies:
            if dep.missing_features > n:
                self.drop_fit_dependence(dep)
            else:
                keep.append(dep)
        self.dependencies = set(keep)

    def best_exact_fits(self):
        by_feature_set = defaultdict(list)
        best_fits = []
        for fit in self.dependencies:
            by_feature_set[tuple(fit.features)].append(fit)
        for peak_tuple, fits in by_feature_set.items():
            fits = sorted(fits, key=lambda x: x.score,
                          reverse=not self.maximize)
            for fit in fits[:-1]:
                self.drop_fit_dependence(fit)
            best_fits.append(fits[-1])
        self.dependencies = set(best_fits)

    def _gather_independent_clusters(self, nodes_for_cache=None):
        clusters = defaultdict(set)
        if nodes_for_cache is None:
            nodes_for_cache = {}

        for node in self.nodes.values():
            # This feature is depended upon by each fit in `dependencies`
            dependencies = set(node.links.keys())

            if len(dependencies) == 0:
                continue

            # These fits also depend upon these other features, and those features are depended upon
            # for other fits in turn, which depend upon the assignment of this feature.
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


try:
    _FeatureDependenceGraphBase = FeatureDependenceGraphBase
    has_c = True
    from ms_deisotope._c.feature_map.dependence_network import FeatureDependenceGraphBase
except ImportError:
    has_c = False


class FeatureDependenceGraph(FeatureDependenceGraphBase):
    def __init__(self, feature_map, nodes=None, dependencies=None, max_missing_features=1,
                 use_monoisotopic_superceded_filtering=True, maximize=True):
        if nodes is None:
            nodes = {}
        if dependencies is None:
            dependencies = set()
        self.feature_map = feature_map
        self.nodes = nodes
        self.dependencies = dependencies
        self.max_missing_features = max_missing_features
        self.use_monoisotopic_superceded_filtering = use_monoisotopic_superceded_filtering
        self.maximize = maximize
        if len(self.nodes) == 0:
            self._populate_initial_graph()
        self.clusters = None
        self._interval_tree = None
        self._solution_map = {}

    def reset(self):
        self.nodes = dict()
        self.dependencies = set()
        self._interval_tree = None
        self._solution_map = dict()
        self._populate_initial_graph()

    def add_solution(self, key, solution):
        self._solution_map[key] = solution

    @property
    def interval_tree(self):
        if self._interval_tree is None:
            self._interval_tree = IntervalTreeNode.build(self.clusters)
            if self._interval_tree is None:
                raise ValueError(
                    "Could not build intervals for peak retrieval with %d clusters" % len(
                        self.clusters), self)
        return self._interval_tree

    def _deep_fuzzy_solution_for(self, peak, shift=0.5):
        pass

    def _find_fuzzy_solution_for(self, feature, shift=0.5):
        tree = self.interval_tree
        clusters = tree.contains_point(feature.mz + shift)
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
                err = abs(fit.monoisotopic_feature.mz - feature.mz)
                if err < error:
                    error = err
                    index = i
            fit = best_fits[index]
            return self._solution_map[fit]

    def find_solution_for(self, feature):
        feature_node = self.nodes[feature]
        tree = self.interval_tree
        if tree is None:
            tree = IntervalTreeNode.build(self.clusters)
            if tree is None:
                raise ValueError(
                    "Could not build intervals for peak retrieval with %d clusters" % len(
                        self.clusters))
        clusters = tree.contains_point(feature.mz)
        if len(clusters) == 0:
            return self._find_fuzzy_solution_for(feature)
        best_fits = [cluster.disjoint_best_fits() for cluster in clusters]

        acc = []
        for fits in best_fits:
            acc.extend(fits)
        best_fits = acc

        common = tuple(set(best_fits) & set(feature_node.links))

        if len(common) > 1 or len(common) == 0:
            if len(common) > 1:
                warnings.warn("Too many solutions exist for %r" % feature)
            # If there were no fits for this peak, then it may be that this peak
            # was not included in a fit. Try to find the nearest solution.
            i = 0
            err = float('inf')
            for j, case in enumerate(best_fits):
                case_err = abs(case.monoisotopic_feature.mz - feature.mz)
                if case_err < err:
                    i = j
                    err = case_err
            fit = best_fits[i]
        else:
            fit = common[0]
        return self._solution_map[fit]

    def _populate_initial_graph(self):
        for feature in self.feature_map:
            self.nodes[feature] = FeatureNode(feature)

    def find_non_overlapping_intervals(self):
        self.drop_gapped_fits()
        self.best_exact_fits()
        if self.use_monoisotopic_superceded_filtering:
            self.drop_superceded_fits()

        nodes_for_cache = {}

        clusters = self._gather_independent_clusters(nodes_for_cache)

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
        return "%s(%s, %d)" % (self.__class__.__name__, self.feature_map, len(self.dependencies))
