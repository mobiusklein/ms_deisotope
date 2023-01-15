cimport cython
from cpython.list cimport PyList_GetItem, PyList_Size, PyList_AsTuple
from cpython.tuple cimport PyTuple_GetItem, PyTuple_Size
from cpython.set cimport PySet_Add

from ms_deisotope._c.peak_dependency_network.intervals cimport (SpanningMixin, Interval)
from ms_deisotope._c.scoring cimport IsotopicFitRecord

from ms_peak_picker._c.peak_set cimport FittedPeak


@cython.final
cdef class FitNode(SpanningMixin):
    def __init__(self, IsotopicFitRecord fit, size_t index=-1):
        self.index = index
        self.fit = fit
        self.edges = set()
        self.overlap_edges = set()
        self.peak_indices = set()
        self._hash = hash(fit)
        self.score = fit.score
        self._init_fields(fit.experimental)

    @staticmethod
    cdef FitNode _create(IsotopicFitRecord fit, size_t index):
        cdef:
            FitNode node
        node = FitNode.__new__(FitNode)
        node.index = index
        node.fit = fit
        node.edges = set()
        node.overlap_edges = set()
        node.peak_indices = set()
        node._hash = hash(fit)
        node.score = fit.score
        node._init_fields(fit.experimental)
        return node

    cdef void _init_fields(self, list experimental):
        cdef:
            size_t i, n
            FittedPeak peak
        n = PyList_Size(experimental)
        for i in range(n):
            # Placeholder Peaks have a peak_count of -1
            peak = <FittedPeak>PyList_GetItem(experimental, i)
            if peak.peak_count >= 0:
                PySet_Add(self.peak_indices, peak)
        if n > 0:
            self.end = peak.mz
            self.start = (<FittedPeak>PyList_GetItem(experimental, 0)).mz
        else:
            self.end = 0
            self.start = 0

    def __hash__(self):
        return self._hash

    def __eq__(self, FitNode other):
        return self.fit is other.fit

    def __gt__(self, FitNode other):
        return self.score > other.score

    def __lt__(self, FitNode other):
        return self.score < other.score

    def __ne__(self, FitNode other):
        return self.fit is not other.fit

    cpdef bint isdisjoint(self, FitNode other):
        return self.peak_indices.isdisjoint(other.peak_indices)

    cpdef visit(self, FitNode other):
        if self.isdisjoint(other):
            PySet_Add(self.edges, other)
            PySet_Add(other.edges, self)
        else:
            PySet_Add(self.overlap_edges, other)
            PySet_Add(other.overlap_edges, self)

    def __repr__(self):
        return "FitNode(%r)" % (self.fit, )


cpdef bint span_overlap(FitNode a, FitNode b):
    return a.__contains__(b.start) or a.__contains__(b.end - 1) or\
        b.__contains__(a.start) or b.__contains__(a.end - 1)


cpdef bint peak_overlap(FitNode a, FitNode b):
    return len(a.peak_indices & b.peak_indices) > 0


@cython.nonecheck(False)
cpdef list layout_layers(list envelopes, bint maximize=True):
    """
    Produce a non-overlapping stacked layout of individual envelopes.
    """
    cdef:
        list layers, layer
        size_t env_i, env_n, layer_i, layer_n, member_i, member_n
        bint placed, collision
        FitNode envelope, member

    layers = [[]]
    envelopes.sort(reverse=maximize)
    env_n = PyList_Size(envelopes)
    for env_i in range(env_n):
        envelope = <FitNode>PyList_GetItem(envelopes, env_i)
        placed = False
        layer_n = PyList_Size(layers)
        for layer_i in range(layer_n):
            layer = <list>PyList_GetItem(layers, layer_i)
            collision = False
            member_n = PyList_Size(layer)
            for member_i in range(member_n):
                member = <FitNode>PyList_GetItem(layer, member_i)
                if peak_overlap(envelope, member):
                    collision = True
                    break
            if not collision:
                layer.append(envelope)
                placed = True
                break
        if not placed:
            layers.append([envelope])
    return layers


@cython.final
@cython.freelist(5)
cdef class GreedySubgraphSelection(object):
    cdef:
        public list intervals
        public list _layers
        public bint maximize

    def __init__(self, object subgraph, bint maximize=True):
        self.intervals = list(subgraph)
        self._layers = [[]]
        self.maximize = maximize
        self._build_layers()

    cdef void _build_layers(self):
        self._layers = layout_layers(self.intervals, maximize=self.maximize)

    cdef list _select_best_subset(self):
        return self._layers[0]

    cpdef select(self):
        return self._select_best_subset()


@cython.freelist(5)
cdef class ConnectedSubgraph(object):
    def __init__(self, object fits, bint maximize=True):
        self.maximize = maximize
        self._init_nodes(list(fits))
        self.populate_edges()

    cdef void _init_nodes(self, list fits):
        cdef:
            list nodes
            size_t i, n
            IsotopicFitRecord fit
        nodes = list()
        n = PyList_Size(fits)
        for i in range(n):
            fit = <IsotopicFitRecord>PyList_GetItem(fits, i)
            nodes.append(FitNode._create(fit, i))
        self.nodes = PyList_AsTuple(nodes)

    cdef void populate_edges(self):
        cdef:
            size_t i, j, n
            FitNode node, other

        n = PyTuple_Size(self.nodes)

        for i in range(n):
            node = <FitNode>PyTuple_GetItem(self.nodes, i)
            for j in range(n):
                if i == j:
                    continue
                other = <FitNode>PyTuple_GetItem(self.nodes, j)
                node.visit(other)

    def __getitem__(self, i):
        return self.nodes[i]

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)

    cpdef find_heaviest_path(self, method="greedy"):
        cdef:
            GreedySubgraphSelection solution
        if len(self) == 1:
            return set(self.nodes)
        if method == "greedy":
            solution = GreedySubgraphSelection(self, maximize=self.maximize)
            return solution.select()
        else:
            raise NotImplementedError(method)
