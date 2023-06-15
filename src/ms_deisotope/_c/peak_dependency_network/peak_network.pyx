import warnings

from collections import defaultdict

cimport cython

from cpython cimport Py_INCREF
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE, PyList_New, PyList_SET_ITEM
from cpython.tuple cimport PyTuple_GET_ITEM
from cpython.sequence cimport PySequence_Fast, PySequence_Fast_ITEMS
from cpython.int cimport PyInt_AsLong
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem, PyDict_DelItem, PyDict_Next, PyDict_Values
from cpython.object cimport PyObject
from cpython.set cimport PySet_Add

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

    def __init__(self, parent=None, dependencies=None, maximize=True):
        if parent is None:
            parent = self
        if dependencies is None:
            dependencies = []
        else:
            dependencies = sorted(dependencies, reverse=maximize)
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

    cpdef disjoint_subset(self):
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

    def uses_mz(self, double mz):
        cdef:
            FittedPeak peak
            IsotopicFitRecord fit
        for fit in self:
            for peak in fit.experimental:
                if peak.mz == mz:
                    return True
        return False

    def fits_using_mz(self, double mz):
        cdef:
            FittedPeak peak
            IsotopicFitRecord fit
        fits = []
        for fit in self:
            for peak in fit.experimental:
                if peak.mz == mz:
                    fits.append(fit)
                    break
        return fits


cdef class PeakDependenceGraphBase(object):

    cpdef reset(self):
        # Keep a record of all clusters from previous iterations
        self._all_clusters.extend(
            self.clusters if self.clusters is not None else [])
        self.nodes = dict()
        self.dependencies = set()
        self._interval_tree = None
        self._populate_initial_graph()

    cpdef _populate_initial_graph(self):
        cdef:
            size_t i, n
            FittedPeak peak
        n = self.peaklist.get_size()
        for i in range(n):
            peak = self.peaklist.getitem(i)

        for peak in self.peaklist:
            self.nodes[peak.index] = PeakNode(peak)

    cpdef add_fit_dependence(self, IsotopicFitRecord fit_record):
        """Add the relatoinship between the experimental peaks
        in `fit_record` to the graph, expressed as a hyper-edge
        denoted by `fit_record`.

        This adds `fit_record` to :attr:`PeakNode.links` for each
        node corresponding to the :class:`~.FittedPeak` instance in
        :attr:`IsotopicFitRecord.experimental` of `fit_record`. It also
        adds `fit_record to :attr:`dependencies`

        Parameters
        ----------
        fit_record: :class:`~.IsotopicFitRecord`
        """
        cdef:
            size_t i, n
            FittedPeak peak
            PeakNode node
        n = PyList_GET_SIZE(fit_record.experimental)
        for i in range(n):
            peak = <FittedPeak>PyList_GET_ITEM(fit_record.experimental, i)
            if peak.index == 0 and peak.intensity <= 1:
                continue
            node = <PeakNode>PyDict_GetItem(self.nodes, peak.index)
            PyDict_SetItem(node.links, fit_record, fit_record.score)
        self.dependencies.add(fit_record)

    cpdef list nodes_for(self, IsotopicFitRecord fit_record, dict cache=None):
        cdef:
            list result
            size_t i, n
            FittedPeak p
            PeakNode node
        if cache is None:
            result = []
            n = PyList_GET_SIZE(fit_record.experimental)
            for i in range(n):
                p = <FittedPeak>PyList_GET_ITEM(fit_record.experimental, i)
                if p.peak_count >= 0:
                    node = <PeakNode>PyDict_GetItem(self.nodes, p.index)
                    result.append(node)
            return result
        else:
            try:
                return cache[fit_record]
            except KeyError:
                result = []
                n = PyList_GET_SIZE(fit_record.experimental)
                for i in range(n):
                    p = <FittedPeak>PyList_GET_ITEM(fit_record.experimental, i)
                    if p.peak_count >= 0:
                        node = <PeakNode>PyDict_GetItem(self.nodes, p.index)
                        result.append(node)
                PyDict_SetItem(cache, fit_record, result)
                return result

    cpdef drop_fit_dependence(self, IsotopicFitRecord fit_record):
        """Remove this fit from the graph, deleting all
        hyper-edges.
        """
        cdef:
            list nodes
            size_t i, n

        nodes = self.nodes_for(fit_record)
        n = PyList_GET_SIZE(nodes)
        for i in range(n):
            node = <PeakNode>PyList_GET_ITEM(nodes, i)
            try:
                (PyDict_DelItem(node.links, fit_record))
            except KeyError:
                pass

    cpdef best_exact_fits(self):
        """For each distinct group of experimental peaks, retain only
        the best scoring fit using exactly those peaks.
        """
        cdef:
            dict by_peaks
            list best_fits, fits_for
            IsotopicFitRecord fit
            object _fit
            tuple peaks_key
            Py_ssize_t i, j, n
            size_t best_index
            PyObject* ptemp
            PyObject* pvalue

        by_peaks = {}
        best_fits = []
        for _fit in self.dependencies:
            fit = <IsotopicFitRecord>_fit
            peaks_key = tuple(fit.experimental)
            ptemp = PyDict_GetItem(by_peaks, peaks_key)
            if ptemp == NULL:
                fits_for = [fit]
                PyDict_SetItem(by_peaks, peaks_key, fits_for)
            else:
                fits_for = <list>ptemp
                fits_for.append(fit)
        j = 0
        while PyDict_Next(by_peaks, &j, &ptemp, &pvalue):
            fits_for = <list>pvalue
            n = PyList_GET_SIZE(fits_for)
            if n == 0:
                continue
            best_index = best_fit_index(fits_for, n, self.maximize)
            for i in range(0, n):
                if i != best_index:
                    fit = <IsotopicFitRecord>PyList_GET_ITEM(fits_for, i)
                    self.drop_fit_dependence(fit)
            fit = <IsotopicFitRecord>PyList_GET_ITEM(fits_for, best_index)
            best_fits.append(fit)
        self.dependencies = set(best_fits)

    cpdef _gather_independent_clusters(self, dict nodes_for_cache=None):
        cdef:
            dict clusters
            set dependencies
            list dependencies_snapshot, dependent_nodes
            PeakNode seed_node, dep_node
            list nodes, nodes_of
            IsotopicFitRecord dep
            Py_ssize_t i, n, j, m, k, p
            PyObject* ptemp

        clusters = dict()
        if nodes_for_cache is None:
            nodes_for_cache = dict()

        nodes = PyDict_Values(self.nodes)
        n = PyList_GET_SIZE(nodes)
        for i in range(n):
            seed_node = <PeakNode>PyList_GET_ITEM(nodes, i)

            # This peak is depended upon by each fit in `dependencies`
            dependencies = set(seed_node.links)

            if len(dependencies) == 0:
                continue

            # These fits also depend upon these other peaks, and those peaks are depended upon
            # for other fits in turn, which depend upon the assignment of this peak.
            dependencies_snapshot = list(dependencies)
            dependent_nodes = []
            m = PyList_GET_SIZE(dependencies_snapshot)
            for j in range(m):
                dep = <IsotopicFitRecord>PyList_GET_ITEM(dependencies_snapshot, j)
                nodes_of = self.nodes_for(dep, nodes_for_cache)
                p = PyList_GET_SIZE(nodes_of)
                for k in range(p):
                    dep_node = <PeakNode>PyList_GET_ITEM(nodes_of, k)
                    ptemp = PyDict_GetItem(clusters, dep_node)
                    if ptemp != NULL:
                        dependencies |= (<set>ptemp)
                    # Update all co-depended nodes with the full set of all fits which depend upon them
                    # PyDict_SetItem(clusters, dep_node, dependencies)

            # dependencies = set(dependencies)

            for _dep in dependencies:
                dep = <IsotopicFitRecord>_dep
                nodes_of = self.nodes_for(dep, nodes_for_cache)
                p = PyList_GET_SIZE(nodes_of)
                for k in range(p):
                    dep_node = <PeakNode>PyList_GET_ITEM(nodes_of, k)
                    PyDict_SetItem(clusters, dep_node, dependencies)

        return clusters


cdef double INFTY = float('inf')

cdef size_t best_fit_index(list fits_for, size_t n, bint maximize):
    cdef:
        double best_score
        size_t i, best_index
        PyObject* value

    if maximize:
        best_score = -INFTY
    else:
        best_score = INFTY
    best_index = 0
    for i in range(n):
        value = <PyObject*>PyList_GET_ITEM(fits_for, i)
        if (maximize and (<IsotopicFitRecord>value).score > best_score) or\
           (not maximize and (<IsotopicFitRecord>value).score < best_score):
            best_score = (<IsotopicFitRecord>value).score
            best_index = i
    return best_index