# cython: embedsignature=True
from collections import defaultdict

cimport cython
from cpython cimport PyObject
from cpython cimport PyList_Append, PyList_Size, PyList_GetItem
from cpython cimport PySet_Add
from cpython cimport PyDict_GetItem, PyDict_SetItem
from cpython cimport PyInt_AsLong

import numpy as np
cimport numpy as np

from ms_deisotope._c.peak_set cimport DeconvolutedPeak

np.import_array()


cdef class NodeBase(object):
    def __init__(self):
        self.key = self.make_key()

    cdef double get_neutral_mass(self):
        return 0

    cdef double get_intensity(self):
        return 0

    cdef long make_key(self):
        return 1

    cdef bint _eq(self, NodeBase other):
        if other is None:
            return False
        elif abs(self.get_neutral_mass() - other.get_neutral_mass()) >= 1e-4:
            return False
        elif abs(self.get_intensity() - other.get_intensity()) >= 1e-4:
            return False
        elif self.key != other.key:
            return False
        return True

    @property
    def neutral_mass(self):
        return self.get_neutral_mass()

    @property
    def intensity(self):
        return self.get_intensity()

    def __repr__(self):
        return "{self.__class__.__name__}({mass}, {intensity})".format(
            self=self, mass=self.neutral_mass, intensity=self.intensity)


cdef class PeakNode(NodeBase):

    def __init__(self, peak):
        self.peak = peak
        super(PeakNode, self).__init__()

    cdef double get_neutral_mass(self):
        return self.peak.neutral_mass

    cdef double get_intensity(self):
        return self.peak.intensity

    cdef long make_key(self):
        return self.peak._index.neutral_mass

    def __iter__(self):
        yield self.peak


cdef class PeakGroupNode(NodeBase):
    def __init__(self, peaks):
        self.peaks = sorted(peaks, key=lambda x: x.neutral_mass)
        super(PeakGroupNode, self).__init__()

    cdef size_t get_size(self):
        return len(self.peaks)

    cdef double get_neutral_mass(self):
        cdef:
            size_t i, n
            double acc, norm
            DeconvolutedPeak peak
        n = self.get_size()
        acc = 0
        norm = 0
        for i in range(n):
            peak = self.peaks[i]
            acc += peak.neutral_mass * (peak.intensity + 1)
            norm += peak.intensity + 1
        if n == 0 or norm == 0:
            return 0
        return acc / norm

    cdef double get_intensity(self):
        cdef:
            size_t i, n
            double acc
            DeconvolutedPeak peak
        n = self.get_size()
        acc = 0
        for i in range(n):
            peak = self.peaks[i]
            acc += (peak.intensity)
        return acc

    cdef long make_key(self):
        cdef:
            size_t i, n
            long acc
            DeconvolutedPeak peak
        n = self.get_size()
        acc = 123823903
        for i in range(n):
            peak = self.peaks[i]
            acc += (peak._index.neutral_mass)
            acc = acc << 1
        return acc

    def __iter__(self):
        return iter(self.peaks)


@cython.final
@cython.freelist(1000000)
cdef class Edge(object):

    @staticmethod
    cdef Edge _create(NodeBase start, NodeBase end, object annotation):
        cdef Edge self = Edge.__new__(Edge)
        self.start = start
        self.end = end
        self.annotation = annotation
        self.key = (self.start.key, self.end.key)
        self._hash = hash(self.key)
        return self

    def __init__(self, start, end, annotation):
        self.start = start
        self.end = end
        self.annotation = annotation
        self.key = (self.start.key, self.end.key)
        self._hash = hash(self.key)

    def __reduce__(self):
        return self.__class__, (self.start, self.end, self.annotation)

    def __eq__(self, other):
        cdef:
            Edge typed_other
        if other is None:
            return False
        else:
            typed_other = <Edge?>other
        # if self.key != typed_other.key:
        #     return False
        if not self.start._eq(typed_other.start):
            return False
        elif not self.end._eq(typed_other.end):
            return False
        elif self.annotation != typed_other.annotation:
            return False
        return True

    def __iter__(self):
        yield self.start
        yield self.end
        yield self.annotation

    def __hash__(self):
        return self._hash

    def __repr__(self):
        return ("{self.__class__.__name__}({self.start.neutral_mass:0.3f} "
                "-{self.annotation}-> {self.end.neutral_mass:0.3f})").format(self=self)


cdef class SpectrumGraph(object):
    def __init__(self):
        self.transitions = set()
        self.by_first = defaultdict(list)
        self.by_second = defaultdict(list)

    cpdef add(self, DeconvolutedPeak p1, DeconvolutedPeak p2, object annotation):
        if p1._index.neutral_mass > p2._index.neutral_mass:
            temp = p2
            p2 = p1
            p1 = temp
        trans = Edge._create(PeakNode(p1), PeakNode(p2), annotation)
        PySet_Add(self.transitions, trans)
        self.by_first[trans.key[0]].append(trans)
        self.by_second[trans.key[1]].append(trans)

    def __iter__(self):
        return iter(self.transitions)

    def __len__(self):
        return len(self.transitions)

    def __repr__(self):
        return "{self.__class__.__name__}({size})".format(
            self=self, size=len(self.transitions))

    def _get_maximum_index(self):
        try:
            trans = max(self.transitions, key=lambda x: x.key[1])
            return trans.key[1]
        except ValueError:
            return 0

    def adjacency_matrix(self):
        n = self._get_maximum_index() + 1
        A = np.zeros((n, n))
        for trans in self.transitions:
            A[trans.key] = 1
        return A

    cpdef list topological_sort(self, adjacency_matrix=None):
        if adjacency_matrix is None:
            adjacency_matrix = self.adjacency_matrix()
        else:
            adjacency_matrix = adjacency_matrix.copy()

        waiting = set()
        for i in range(adjacency_matrix.shape[0]):
            # Check for incoming edges. If no incoming
            # edges, add to the waiting set
            if adjacency_matrix[:, i].sum() == 0:
                waiting.add(i)
        ordered = list()
        while waiting:
            ix = waiting.pop()
            ordered.append(ix)
            outgoing_edges = adjacency_matrix[ix, :]
            indices_of_neighbors = np.nonzero(outgoing_edges)[0]
            # For each outgoing edge
            for neighbor_ix in indices_of_neighbors:
                # Remove the edge
                adjacency_matrix[ix, neighbor_ix] = 0
                # Test for incoming edges
                if (adjacency_matrix[:, neighbor_ix] == 0).all():
                    waiting.add(neighbor_ix)
        if adjacency_matrix.sum() > 0:
            raise ValueError("%d edges left over" % (adjacency_matrix.sum(),))
        else:
            return ordered

    def path_lengths(self):
        adjacency_matrix = self.adjacency_matrix()
        distances = np.zeros(self._get_maximum_index() + 1)
        for ix in self.topological_sort():
            incoming_edges = np.nonzero(adjacency_matrix[:, ix])[0]
            if incoming_edges.shape[0] == 0:
                distances[ix] = 0
            else:
                longest_path_so_far = distances[incoming_edges].max()
                distances[ix] = longest_path_so_far + 1
        return distances

    cpdef list paths_starting_at(self, ix, int limit=-1):
        paths = []
        for trans in self.by_first[ix]:
            paths.append([trans])
        finished_paths = []
        while True:
            extended_paths = []
            for path in paths:
                terminal = path[-1]
                edges = self.by_first[terminal.key[1]]
                if not edges:
                    finished_paths.append(path)
                for trans in edges:
                    extended_paths.append(path + [trans])
            paths = extended_paths
            if len(paths) == 0:
                break
        return self.transitive_closure(finished_paths, limit)

    cpdef list paths_ending_at(self, ix, int limit=-1):
        paths = []
        for trans in self.by_second[ix]:
            paths.append([trans])
        finished_paths = []
        while True:
            extended_paths = []
            for path in paths:
                terminal = path[0]
                edges = self.by_second[terminal.key[0]]
                if not edges:
                    finished_paths.append(path)
                for trans in edges:
                    new_path = [trans]
                    new_path.extend(path)
                    extended_paths.append(new_path)
            paths = extended_paths
            if len(paths) == 0:
                break
        return self.transitive_closure(finished_paths, limit)

    cpdef list transitive_closure(self, list paths, int limit=-1):
        cdef:
            list path, index_list
            list node_sets, keep, node_sets_items
            set node_set, other_node_set
            dict by_node
            PyObject* tmp
            object conv
            Edge node
            size_t i, j, n, m, k, q, v
            bint is_enclosed
        paths = sorted(paths, key=len, reverse=True)
        # precompute node_sets for each path
        node_sets = []
        by_node = {}
        for i in range(PyList_Size(paths)):
            path = <list>PyList_GetItem(paths, i)
            # track all nodes by index in this path in node_set
            node_set = set()
            for j in range(PyList_Size(path)):
                node = <Edge>PyList_GetItem(path, j)
                # add the start node index and end node index to the
                # set of nodes on this path
                # node_set.update(node.key)
                if j == 0:
                    conv = node.start._index.neutral_mass
                    node_set.add(conv)
                    tmp = PyDict_GetItem(by_node, conv)
                    if tmp == NULL:
                        index_list = []
                        PyDict_SetItem(by_node, conv, index_list)
                    else:
                        index_list = <list>tmp
                    index_list.append(i)
                conv = node.end._index.neutral_mass
                node_set.add(conv)
                tmp = PyDict_GetItem(by_node, conv)
                if tmp == NULL:
                    index_list = []
                    PyDict_SetItem(by_node, conv, index_list)
                else:
                    index_list = <list>tmp
                index_list.append(i)
            node_sets.append(node_set)

        keep = []
        k = 0
        for i in range(PyList_Size(paths)):
            path = <list>PyList_GetItem(paths, i)
            is_enclosed = False
            node_set = <set>PyList_GetItem(node_sets, i)
            n = len(node_set)
            node = <Edge>PyList_GetItem(path, 0)
            conv = node.start._index.neutral_mass
            index_list = <list>PyDict_GetItem(by_node, conv)
            for q in range(PyList_Size(index_list)):
                j = PyInt_AsLong(<object>PyList_GetItem(index_list, q))
                if i == j:
                    continue
                other_node_set = <set>PyList_GetItem(node_sets, j)
                m = len(other_node_set)
                if m < n:
                    break
                if node_set < other_node_set:
                    is_enclosed = True
                    break
            if not is_enclosed:
                keep.append(path)
                k += 1
            if limit > 0 and limit == k:
                break
        return keep

    def longest_paths(self, int limit=-1):
        cdef:
            list paths, segment
        # get all distinct paths
        paths = []
        for ix in np.argsort(self.path_lengths())[::-1]:
            segment = self.paths_ending_at(ix, limit)
            paths.extend(segment)
            if PyList_Size(segment) == 0:
                break
            # remove redundant paths
            paths = self.transitive_closure(paths, limit)
            if limit > 0:
                if PyList_Size(paths) > limit:
                    break
        return paths

@cython.final
@cython.freelist(1000)
cdef class MassWrapper(object):
    '''An adapter class to make types whose mass calculation is a method
    (:mod:`glypy` dynamic graph components) compatible with code where the
    mass calculation is an  attribute (:mod:`glycopeptidepy` objects and
    most things here)

    Hashes and compares as :attr:`obj`

    Attributes
    ----------
    obj: object
        The wrapped object
    mass: float
        The mass of :attr:`obj`
    '''

    def __init__(self, obj):
        self.obj = obj
        try:
            # object's mass is a method
            self.mass = obj.mass()
        except TypeError:
            # object's mass is a plain attribute
            self.mass = obj.mass

    def __repr__(self):
        return "{self.__class__.__name__}({self.obj})".format(self=self)

    def __eq__(self, other):
        return self.obj == other

    def __hash__(self):
        return hash(self.obj)


cdef class Path(object):
    def __init__(self, edge_list):
        self.transitions = edge_list
        self.total_signal = self._total_signal()
        self.start_mass = self[0].start.neutral_mass
        self.end_mass = self[-1].end.neutral_mass
        self._peaks_used = None
        self._edges_used = None

    @property
    def peaks(self):
        if self._peaks_used is None:
            self._peaks_used = self._build_peaks_set()
        return self._peaks_used

    def _build_peaks_set(self):
        peaks = set()
        for edge in self:
            peaks.update(edge.end)
        peaks.update(self[0].start)
        return peaks

    def _build_edges_used(self):
        mapping = defaultdict(list)
        for edge in self:
            mapping[edge.start, edge.end].append(edge)
        return mapping

    def __contains__(self, edge):
        return self.has_edge(edge, True)

    def has_edge(self, edge, match_annotation=False):
        if self._edges_used is None:
            self._edges_used = self._build_edges_used()
        key = (edge.start, edge.end)
        if key in self._edges_used:
            edges = self._edges_used[key]
            if match_annotation:
                for e in edges:
                    if e == edge:
                        return True
                return False
            else:
                for e in edges:
                    if e.key != edge.key:
                        continue
                    if e.start != edge.start:
                        continue
                    if e.end != edge.end:
                        continue
                    return True
                return False

    def has_peak(self, peak):
        if self._peaks_used is None:
            self._peaks_used = self._build_peaks_set()
        return peak in self._peaks_used

    def has_peaks(self, peaks_set):
        if self._peaks_used is None:
            self._peaks_used = self._build_peaks_set()
        return peaks_set & self._peaks_used

    def __iter__(self):
        return iter(self.transitions)

    def __getitem__(self, i):
        return self.transitions[i]

    def __len__(self):
        return self.get_size()

    cdef Edge get(self, Py_ssize_t i):
        if i < 0:
            return self.transitions[i]
        return <Edge>PyList_GetItem(self.transitions, i)

    cdef size_t get_size(self):
        return PyList_Size(self.transitions)

    cdef double _total_signal(self):
        total = 0
        for edge in self:
            total += edge.end.intensity
        total += self[0].start.intensity
        return total

    def __repr__(self):
        return "%s(%s, %0.4e, %f, %f)" % (
            self.__class__.__name__,
            '->'.join(str(e.annotation) for e in self),
            self.total_signal, self.start_mass, self.end_mass
        )
