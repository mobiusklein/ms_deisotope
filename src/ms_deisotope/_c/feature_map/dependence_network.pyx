import operator

from cpython cimport PyObject
from cpython.list cimport PyList_GetItem, PyList_Size, PyList_GET_SIZE, PyList_GET_ITEM, PyList_Append
from cpython.exc cimport PyErr_Occurred, PyErr_Clear, PyErr_ExceptionMatches
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_Contains, PyDict_DelItem, PyDict_Next, PyDict_Values
from cpython.tuple cimport PyTuple_GetItem, PyTuple_Size
from cpython.set cimport PySet_Size

from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin, IntervalTreeNode

from ms_deisotope._c.feature_map.lcms_feature cimport FeatureBase, EmptyFeature
from ms_deisotope._c.feature_map.feature_fit cimport LCMSFeatureSetFit, map_coord

from ms_deisotope.peak_dependency_network.subgraph import GreedySubgraphSelection


cdef class SpanningMixin2D(object):
    """Provides methods for checking whether an entity
    which has a defined start and end point over a single
    dimension contains or overlaps with another entity in
    that same dimension.
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __contains__(self, i):
        """Tests for point inclusion, `start <= i <= end`

        Parameters
        ----------
        i : Number
            The point to be tested

        Returns
        -------
        bool
        """
        return self._contains(i)

    cdef bint _contains(self, map_coord i):
        """Tests for point inclusion, `start <= i <= end`

        Parameters
        ----------
        i : map_coord
            The point to be tested

        Returns
        -------
        bool
        """
        return self.start.ge(i) and i.le(self.end)

    cpdef bint contains(self, map_coord i):
        """Tests for point inclusion, `start <= i <= end`

        Parameters
        ----------
        i : map_coord
            The point to be tested

        Returns
        -------
        bool
        """
        return self._contains(i)

    cpdef bint overlaps(self, SpanningMixin2D interval):
        """Tests whether another spanning entity with
        a defined start and end point overlaps with
        this spanning entity

        Parameters
        ----------
        interval : SpanningMixin

        Returns
        -------
        bool
        """
        cdef bint cond
        cond = (
            (self.start.le(interval.start) and self.end.ge(interval.start)) or \
            (self.start.ge(interval.start) and self.end.le(interval.end)) or \
            (self.start.ge(interval.start) and self.end.ge(interval.end) and self.start.ge(interval.end)) or \
            (self.start.le(interval.start) and self.end.ge(interval.start)) or \
            (self.start.le(interval.end) and self.end.ge(interval.end)))
        return cond

    cpdef bint is_contained_in_interval(self, SpanningMixin2D interval):
        """Tests whether this spanning entity is
        completely contained inside the other entity

        Parameters
        ----------
        interval : SpanningMixin

        Returns
        -------
        bool
        """
        return self.start.ge(interval.start) and self.end.le(interval.end)

    cpdef bint contains_interval(self, SpanningMixin2D interval):
        """Tests whether the other spanning entity is
        completely contained inside this entity

        Parameters
        ----------
        interval : SpanningMixin

        Returns
        -------
        bool
        """
        return self.start.le(interval.start) and self.end.ge(interval.end)



cpdef double score_getter(LCMSFeatureSetFit x):
    return x.score


cdef class FeatureNode(object):
    """Holds all the information about a single `LCMSFeature` instance
    and all the `LCMSFeatureSetFit` instances that depend upon it.

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

    def __reduce__(self):
        return self.__class__, (self.feature, self.links)

    def __eq__(self, FeatureNode other):
        if other is None:
            return False
        try:
            return self.feature == other.feature
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __contains__(self, fit):
        return fit in self.links

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__, self.feature, self.links)


cdef class DependenceCluster(SpanningMixin2D):
    def __init__(self, parent=None, dependencies=None, maximize=True):
        if parent is None:
            parent = self
        if dependencies is None:
            dependencies = []
        else:
            dependencies = sorted(dependencies, key=score_getter, reverse=maximize)
        self.parent = parent
        self.dependencies = dependencies
        self.best_fit = None
        self.maximize = maximize
        self._reset()

    cpdef _reset(self):
        self.start = self._start()
        self.end = self._end()
        self.best_fit = self._best_fit()

    cpdef add(self, LCMSFeatureSetFit fit):
        """
        Adds a new LCMSFeatureSetFit to this cluster, and ensures the sorted
        property still holds

        Parameters
        ----------
        fit : LCMSFeatureSetFit
            The fit being added to the record

        """
        self.dependencies.append(fit)
        self.dependencies.sort(key=score_getter, reverse=self.maximize)
        self._reset()

    cpdef disjoint_subset(self, int max_size=1000):
        cdef:
            list dependencies
        if len(self.dependencies) > max_size:
            dependencies = self.dependencies[:max_size]
        else:
            dependencies = self.dependencies

        graph = ConnectedSubgraph(dependencies, maximize=self.maximize)
        return graph.find_heaviest_path()

    cpdef list disjoint_best_fits(self, int max_size=1000):
        """
        Compute the best set of disjoint isotopic fits spanning this cluster

        Returns
        -------
        list of LCMSFeatureSetFit
        """
        fit_sets = tuple(self.disjoint_subset(max_size=max_size))
        best_fits = fit_sets
        return [node.fit for node in best_fits]

    cdef LCMSFeatureSetFit _best_fit(self):
        """
        Retrieve the absolute best single fit for the cluster

        Returns
        -------
        name : LCMSFeatureSetFit
            The best scoring LCMSFeatureSetFit in this cluster
        """
        return self.dependencies[0]

    cdef map_coord _start(self):
        """
        Determine the first mz-time coordinate for members of this cluster

        Returns
        -------
        map_coord
        """
        cdef:
            size_t i, n
            map_coord best, cur

        n = PyList_Size(self.dependencies)
        if n == 0:
            return None
        best = (<LCMSFeatureSetFit>PyList_GetItem(self.dependencies, 0)).get_start()
        for i in range(1, n):
            cur = (<LCMSFeatureSetFit>PyList_GetItem(self.dependencies, i)).get_start()
            if cur < best:
                best = cur
        return best.copy()

    cdef map_coord _end(self):
        """
        Determines the last mz-time coordinate for members of this cluster

        Returns
        -------
        map_coord
        """
        cdef:
            size_t i, n
            map_coord best, cur

        n = PyList_Size(self.dependencies)
        if n == 0:
            return None
        best = (<LCMSFeatureSetFit>PyList_GetItem(self.dependencies, 0)).get_end()
        for i in range(1, n):
            cur = (<LCMSFeatureSetFit>PyList_GetItem(self.dependencies, i)).get_end()

            if cur > best:
                best = cur
        return best.copy()

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
        return "%s(dependencies=[%s], start=%r, end=%r)" % (
            self.__class__.__name__,
            '\n'.join(map(str, self.dependencies[:10])), self.start, self.end)

    def __getitem__(self, i):
        return self.dependencies[i]


cdef bint is_valid(FeatureBase feature):
    return feature is not None and not isinstance(feature, EmptyFeature)


cdef double INF = float('infinity')


cdef class FeatureSetFitNode(SpanningMixin2D):
    def __init__(self, fit, index=None):
        self.fit = fit
        self.edges = set()
        self.overlap_edges = set()
        self.index = -1
        if index is not None:
            self.index = index
        self._hash = -1
        self.score = fit.score
        self._init_indices()

    cdef void _init_indices(self):
        cdef:
            map_coord start, end
            set feature_indices
            size_t i, n
            double st, et
            FeatureBase feat

        feature_indices = set()
        start = map_coord._create(self.fit.get_start().mz, 0)
        end = map_coord._create(self.fit.get_end().mz, INF)

        # Build a minimal bounding box over features, one that allows overflow.
        n = PyList_Size(self.fit.features)
        for i in range(n):
            feat = <FeatureBase>PyList_GetItem(self.fit.features, i)
            if not is_valid(feat):
                continue
            st = feat.get_start_time()
            et = feat.get_end_time()
            feature_indices.add((feat.get_mz(), st, et))
            if start.time < st:
                start.time = st
            if end.time > et:
                end.time = et
        self.feature_indices = feature_indices
        self.start = start
        self.end = end

    def __hash__(self):
        if self._hash == -1:
            self._hash = hash(self.fit)
        return self._hash

    def __eq__(self, FeatureSetFitNode other):
        return self.fit is other.fit

    def __gt__(self, FeatureSetFitNode other):
        return self.fit > other.fit

    def __lt__(self, FeatureSetFitNode other):
        return self.fit < other.fit

    def __ne__(self, FeatureSetFitNode other):
        return self.fit is not other.fit

    cpdef visit(self, FeatureSetFitNode other):
        if self.isdisjoint(other):
            self.edges.add(other)
            other.edges.add(self)
        else:
            self.overlap_edges.add(other)
            other.overlap_edges.add(self)

    cpdef bint isdisjoint(self, FeatureSetFitNode other):
        return PySet_Size(self.feature_indices & other.feature_indices)

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self.fit)

    @staticmethod
    def overlap(FeatureSetFitNode a, FeatureSetFitNode b):
        return (a.feature_indices & b.feature_indices)


cdef class ConnectedSubgraph(object):
    def __init__(self, fits, maximize=True):
        self.nodes = tuple(map(FeatureSetFitNode, fits))
        self.maximize = maximize
        self._init()

    cdef void _init(self):
        cdef:
            int i, n
            FeatureSetFitNode node

        n = PyTuple_Size(self.nodes)
        for i in range(n):
            node = <FeatureSetFitNode>PyTuple_GetItem(self.nodes, i)
            node.index = i
        self.populate_edges()

    def __getitem__(self, i):
        return self.nodes[i]

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)

    cpdef populate_edges(self):
        cdef:
            size_t i, j, n
            FeatureSetFitNode node, other
        n = PyTuple_Size(self.nodes)
        for i in range(n):
            node = <FeatureSetFitNode>PyTuple_GetItem(self.nodes, i)
            for j in range(i + 1, n):
                other = <FeatureSetFitNode>PyTuple_GetItem(self.nodes, j)
                node.visit(other)

    cpdef find_heaviest_path(self, method="greedy"):
        if len(self) == 1:
            return set(self.nodes)
        if method == "greedy":
            solution = GreedySubgraphSelection.solve(
                self.nodes, maximize=self.maximize, overlap_fn=FeatureSetFitNode.overlap)
            return solution
        else:
            raise NotImplementedError(method)


cdef class FeatureDependenceGraphBase(object):
    cpdef list nodes_for(self, LCMSFeatureSetFit fit_record, dict cache=None):
        cdef:
            size_t i, n
            FeatureBase p
            list result
            PyObject* ptmp

        if cache is None:
            result = []
            n = PyList_Size(fit_record.features)
            for i in range(n):
                p = <FeatureBase>PyList_GET_ITEM(fit_record.features, i)
                if is_valid(p):
                    result.append(<object>PyDict_GetItem(self.nodes, p))

            return result
        else:
            ptmp = PyDict_GetItem(cache, fit_record)
            if ptmp != NULL:
                return <list>ptmp
            result = []
            n = PyList_Size(fit_record.features)
            for i in range(n):
                p = <FeatureBase>PyList_GET_ITEM(fit_record.features, i)
                if is_valid(p):
                    result.append(<object>PyDict_GetItem(self.nodes, p))
            PyDict_SetItem(cache, fit_record, result)
            return result

    cpdef add_fit_dependence(self, LCMSFeatureSetFit fit_record):
        cdef:
            size_t i, n
            FeatureBase p
            FeatureNode node
            PyObject* ptmp

        n = PyList_Size(fit_record.features)
        for i in range(n):
            p = <FeatureBase>PyList_GET_ITEM(fit_record.features, i)
            if not is_valid(p):
                continue
            node = <FeatureNode>PyDict_GetItem(self.nodes, p)
            PyDict_SetItem(node.links, fit_record, fit_record.score)
        self.dependencies.add(fit_record)

    cpdef drop_fit_dependence(self, LCMSFeatureSetFit fit_record, dict cache=None):
        cdef:
            list nodes
            size_t i, n
            int code
            FeatureNode node
        nodes = self.nodes_for(fit_record, cache=cache)
        n = PyList_Size(nodes)
        for i in range(n):
            node = <FeatureNode>PyList_GET_ITEM(nodes, i)
            code = PyDict_DelItem(node.links, fit_record)
            if code < 0:
                if PyErr_Occurred():
                    if PyErr_ExceptionMatches(KeyError):
                        PyErr_Clear()

    cpdef drop_superceded_fits(self):
        cdef:
            dict cache
            set suppressed, keep
            LCMSFeatureSetFit dep, candidate
            bint suppress, maximize
            FeatureBase monoisotopic_feature
            FeatureNode monoisotopic_feature_node
            object tmp
            PyObject* ptmp

        suppressed = set()
        keep = set()

        maximize = self.maximize
        for tmp in self.dependencies:
            dep = <LCMSFeatureSetFit>tmp
            monoisotopic_feature = dep.monoisotopic_feature
            if not is_valid(monoisotopic_feature):
                continue
            ptmp = PyDict_GetItem(self.nodes, monoisotopic_feature)
            if ptmp == NULL:
                raise KeyError(monoisotopic_feature)
            else:
                monoisotopic_feature_node = <FeatureNode>ptmp
            suppress = False
            for tmp in monoisotopic_feature_node.links:
                candidate = <LCMSFeatureSetFit>tmp
                if dep.charge == candidate.charge:
                    if monoisotopic_feature in candidate.features:
                        if (maximize and (dep.score < candidate.score)) or (not maximize and (dep.score > candidate.score)):
                            suppress = True
                            break
            if suppress:
                suppressed.add(dep)
            else:
                keep.add(dep)
        cache = dict()
        for dropped in suppressed:
            self.drop_fit_dependence(dropped, cache)
        self.dependencies = keep

    cpdef drop_gapped_fits(self, n=None):
        cdef:
            int max_missed
            LCMSFeatureSetFit dep
            set keep
        if n is None:
            max_missed = self.max_missing_features
        else:
            max_missed = n
        keep = set()
        for dep in self.dependencies:
            if dep.missing_features > max_missed:
                self.drop_fit_dependence(dep)
            else:
                keep.add(dep)
        self.dependencies = keep

    cpdef best_exact_fits(self):
        cdef:
            dict by_feature_set
            list bucket, fits
            set best_fits
            LCMSFeatureSetFit fit, best_fit
            tuple key
            PyObject* ptmp
            PyObject* pval
            Py_ssize_t pos
            size_t i, n
            bint maximize

        by_feature_set = {}
        best_fits = set()
        for tmp in self.dependencies:
            fit = <LCMSFeatureSetFit>tmp
            key = tuple(fit.features)
            ptmp = PyDict_GetItem(by_feature_set, key)
            if ptmp == NULL:
                bucket = []
                PyDict_SetItem(by_feature_set, key, bucket)
            else:
                bucket = <list>ptmp
            bucket.append(fit)

        maximize = self.maximize
        pos = 0
        while PyDict_Next(by_feature_set, &pos, &ptmp, &pval):
            fits = <list>pval
            n = PyList_Size(fits)
            if n == 1:
                best_fit = <LCMSFeatureSetFit>PyList_GET_ITEM(fits, 0)
                best_fits.add(best_fit)
                continue
            best_fit = <LCMSFeatureSetFit>PyList_GET_ITEM(fits, 0)
            for i in range(1, n):
                fit = <LCMSFeatureSetFit>PyList_GET_ITEM(fits, i)
                if (maximize and (best_fit.score < fit.score)) or (not maximize and (best_fit.score > fit.score)):
                    self.drop_fit_dependence(best_fit)
                    best_fit = fit
                else:
                    self.drop_fit_dependence(fit)
            best_fits.add(best_fit)
        self.dependencies = best_fits

    cpdef _gather_independent_clusters(self, dict nodes_for_cache=None):
        cdef:
            dict clusters
            set dependencies
            list dependencies_snapshot, dependent_nodes
            FeatureNode seed_node, dep_node
            list nodes, nodes_of
            LCMSFeatureSetFit dep
            Py_ssize_t i, n, j, m, k, p
            PyObject* ptemp

        clusters = dict()
        if nodes_for_cache is None:
            nodes_for_cache = dict()

        nodes = PyDict_Values(self.nodes)
        n = PyList_GET_SIZE(nodes)
        for i in range(n):
            seed_node = <FeatureNode>PyList_GET_ITEM(nodes, i)

            # This feature is depended upon by each fit in `dependencies`
            dependencies = set(seed_node.links)

            m = PySet_Size(dependencies)
            if m == 0:
                continue
            # These fits also depend upon these other peaks, and those peaks are depended upon
            # for other fits in turn, which depend upon the assignment of this peak.
            dependencies_snapshot = list(dependencies)
            dependent_nodes = []
            for j in range(m):
                dep = <LCMSFeatureSetFit>PyList_GET_ITEM(dependencies_snapshot, j)
                nodes_of = self.nodes_for(dep, nodes_for_cache)
                p = PyList_GET_SIZE(nodes_of)
                for k in range(p):
                    dep_node = <FeatureNode>PyList_GET_ITEM(nodes_of, k)
                    ptemp = PyDict_GetItem(clusters, dep_node)
                    if ptemp != NULL:
                        dependencies |= (<set>ptemp)
                    # Update all co-depended nodes with the full set of all fits which depend upon them
                    # PyDict_SetItem(clusters, dep_node, dependencies)
            for _dep in dependencies:
                dep = <LCMSFeatureSetFit>_dep
                nodes_of = self.nodes_for(dep, nodes_for_cache)
                p = PyList_GET_SIZE(nodes_of)
                for k in range(p):
                    dep_node = <FeatureNode>PyList_GET_ITEM(nodes_of, k)
                    PyDict_SetItem(clusters, dep_node, dependencies)

        return clusters
