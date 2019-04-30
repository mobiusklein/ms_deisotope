import operator
import warnings

from collections import defaultdict

cimport cython

from cpython cimport Py_INCREF
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE, PyList_New, PyList_SET_ITEM
from cpython.tuple cimport PyTuple_GET_ITEM
from cpython.int cimport PyInt_AsLong
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.object cimport PyObject
from cpython.set cimport PySet_Add

from ms_peak_picker._c.peak_set cimport FittedPeak, PeakSet

from ms_deisotope._c.scoring cimport IsotopicFitRecord
from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin, IntervalTreeNode
from ms_deisotope._c.peak_dependency_network.subgraph cimport ConnectedSubgraph


cdef object score_getter = operator.attrgetter("score")


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
        cdef:
            size_t i, n
            FittedPeak peak
            PeakNode node
        n = PyList_GET_SIZE(fit_record.experimental)
        for i in range(n):
            peak = <FittedPeak>PyList_GET_ITEM(fit_record.experimental, i)
            if peak.index == 0:
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
