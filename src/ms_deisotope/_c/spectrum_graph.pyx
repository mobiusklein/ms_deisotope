# cython: embedsignature=True
from collections import defaultdict
from itertools import combinations_with_replacement

cimport cython
from cpython cimport PyObject
from cpython cimport PyList_Append, PyList_Size, PyList_GetItem
from cpython cimport PyTuple_Size, PyTuple_GetItem
from cpython cimport PySet_Add
from cpython cimport PyDict_GetItem, PyDict_SetItem
from cpython cimport PyInt_AsLong

import numpy as np
cimport numpy as np

from ms_deisotope._c.peak_set cimport DeconvolutedPeak, DeconvolutedPeakSet

np.import_array()


cdef class NodeBase(object):
    def __init__(self):
        self.key = self.get_index()

    cdef double get_neutral_mass(self):
        return 0

    cdef double get_intensity(self):
        return 0

    cdef size_t get_index(self):
        return 0

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

    def __hash__(self):
        return hash(self.get_index())

    def __eq__(self, other):
        cdef:
            NodeBase other_typed
        if not isinstance(other, NodeBase):
            return False
        other_typed = <NodeBase>other
        return self._eq(other_typed)

    def __ne__(self, other):
        return not self == other

    @property
    def index(self):
        return self.get_index()

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

    cdef size_t get_index(self):
        return self.peak._index.neutral_mass

    cdef double get_neutral_mass(self):
        return self.peak.neutral_mass

    cdef double get_intensity(self):
        return self.peak.intensity

    def __iter__(self):
        yield self.peak


cdef class PeakGroupNode(NodeBase):
    def __init__(self, peaks):
        self.peaks = sorted(set(peaks), key=lambda x: x.neutral_mass)
        super(PeakGroupNode, self).__init__()

    cdef size_t get_size(self):
        return len(self.peaks)

    cdef size_t get_index(self):
        cdef:
            size_t i, n
            size_t index
            DeconvolutedPeak peak
        n = self.get_size()
        index = -1
        for i in range(n):
            peak = <DeconvolutedPeak>PyList_GetItem(self.peaks, i)
            if peak._index.neutral_mass < index:
                index = peak._index.neutral_mass
        return index

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

    def __iter__(self):
        return iter(self.peaks)

cdef class PlaceHolderNode(NodeBase):

    def __init__(self, str label, double mass, double intensity=0, size_t index=0):
        self.label = label
        self._mass = mass
        self._intensity = intensity
        self._index = index
        super(PlaceHolderNode, self).__init__()

    cdef size_t get_index(self):
        return self._index

    cdef double get_neutral_mass(self):
        return self._mass

    cdef double get_intensity(self):
        return self._intensity

    def __iter__(self):
        return iter([])

    def __repr__(self):
        return "{self.__class__.__name__}({self.label!r}, {mass}, {intensity})".format(
            self=self, mass=self.neutral_mass, intensity=self.intensity)


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
                    conv = node.start.get_index()
                    node_set.add(conv)
                    tmp = PyDict_GetItem(by_node, conv)
                    if tmp == NULL:
                        index_list = []
                        PyDict_SetItem(by_node, conv, index_list)
                    else:
                        index_list = <list>tmp
                    index_list.append(i)
                conv = node.end.get_index()
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
            conv = node.start.get_index()
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
    """An adapter class to make types whose mass calculation is a method
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
    """

    @classmethod
    def wrap(cls, components, combinations=1):
        components = list(map(cls, components))
        if combinations > 1:
            result = []
            for i in range(1, combinations + 1):
                for combos in combinations_with_replacement(components, i):
                    mass = sum(c.mass for c in combos)
                    obj = tuple(c.obj for c in combos)
                    result.append(cls(obj, mass))
            components = result
        components = sorted(components, key=lambda x: x.mass)
        return components

    def __init__(self, obj, mass=None):
        if isinstance(obj, MassWrapper):
            if mass is None:
                mass = obj.mass
            obj = obj.obj

        self.obj = obj
        if mass is not None:
            self.mass = mass
        else:
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
        self.transitions = list(edge_list)
        self.total_signal = self._total_signal()
        self.start_mass = self.get(0).start.neutral_mass
        self.end_mass = self.get(-1).end.neutral_mass
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


def create_edge_group(edges):
    if not edges:
        return None
    starts = []
    ends = []
    annotations = []
    for edge in edges:
        starts.extend(edge.start)
        ends.extend(edge.end)
        annotations.append(edge.annotation)
    start = PeakGroupNode(starts)
    end = PeakGroupNode(ends)
    annotation = annotations[0]
    return Edge(start, end, annotation)


def collect_edges(paths):
    edges = []
    n = len(paths)
    if n == 0:
        return []
    m = len(paths[0])
    for i in range(m):
        group = []
        for j in range(n):
            group.append(paths[j][i])
        edges.append(group)
    return Path([create_edge_group(g) for g in edges])


cdef double INF = float('inf')


cdef class PathFinder(object):
    cdef:
        public list components
        public double product_error_tolerance

    def __init__(self, components, product_error_tolerance=1e-5):
        self.components = sorted(map(MassWrapper, components), key=lambda x: x.mass)
        self.product_error_tolerance = product_error_tolerance

    cpdef find_edges(self, scan, max_mass=None):
        cdef:
            SpectrumGraph graph
            size_t pi, pn
            size_t ci, cn
            size_t oi, on
            DeconvolutedPeakSet deconvoluted_peak_set
            DeconvolutedPeak peak, other_peak
            MassWrapper component
            tuple complements
            double upper_limit, query_mass

        graph = SpectrumGraph()
        upper_limit = INF
        deconvoluted_peak_set = self._guess_limit_and_get_peaks(scan, &upper_limit)
        if max_mass is not None:
            upper_limit = max_mass

        cn = PyList_Size(self.components)
        pn = deconvoluted_peak_set.get_size()
        for pi in range(pn):
            peak = <DeconvolutedPeak>deconvoluted_peak_set.getitem(pi)
            if peak.neutral_mass > upper_limit:
                break
            for ci in range(cn):
                component = <MassWrapper>PyList_GetItem(self.components, ci)
                query_mass = peak.neutral_mass + component.mass
                if (query_mass - upper_limit) / upper_limit > self.product_error_tolerance:
                    break
                complements = deconvoluted_peak_set.all_peaks_for(
                    query_mass, self.product_error_tolerance)
                on = PyTuple_Size(complements)
                for oi in range(on):
                    other_peak = <DeconvolutedPeak>PyTuple_GetItem(complements, oi)
                    graph.add(peak, other_peak, component.obj)
        return graph

    cdef DeconvolutedPeakSet _guess_limit_and_get_peaks(self, scan, double* upper_limit):
        cdef:
            DeconvolutedPeakSet deconvoluted_peak_set
        if not isinstance(scan, DeconvolutedPeakSet):
            # Assume we have a Python object with Scan-like interface
            deconvoluted_peak_set = scan.deconvoluted_peak_set
            if scan.precursor_information is not None:
                upper_limit[0] = scan.precursor_information.neutral_mass + 1
            else:
                upper_limit[0] = INF
        else:
            deconvoluted_peak_set = scan
            upper_limit[0] = INF
        return deconvoluted_peak_set

    cpdef add_edges_for(self, scan, SpectrumGraph graph, MassWrapper component, max_mass=None):
        cdef:
            size_t pi, pn
            size_t oi, on
            DeconvolutedPeakSet deconvoluted_peak_set
            DeconvolutedPeak peak, other_peak
            tuple complements
            double upper_limit, query_mass

        upper_limit = INF
        deconvoluted_peak_set = self._guess_limit_and_get_peaks(scan, &upper_limit)
        if max_mass is not None:
            upper_limit = max_mass

        pn = deconvoluted_peak_set.get_size()
        for pi in range(pn):
            peak = <DeconvolutedPeak>deconvoluted_peak_set.getitem(pi)
            if peak.neutral_mass > upper_limit:
                break
            query_mass = peak.neutral_mass + component.mass
            if (query_mass - upper_limit) / upper_limit > self.product_error_tolerance:
                break
            complements = deconvoluted_peak_set.all_peaks_for(
                query_mass, self.product_error_tolerance)
            on = PyTuple_Size(complements)
            for oi in range(on):
                other_peak = <DeconvolutedPeak>PyTuple_GetItem(complements, oi)
                graph.add(peak, other_peak, component.obj)

    cpdef find_complements(self, scan, SpectrumGraph graph,  MassWrapper component, max_mass=None):
        cdef:
            size_t pi, pn
            size_t oi, on
            DeconvolutedPeakSet deconvoluted_peak_set
            DeconvolutedPeak peak, other_peak
            tuple complements
            double upper_limit, query_mass

        upper_limit = INF
        deconvoluted_peak_set = self._guess_limit_and_get_peaks(scan, &upper_limit)
        if max_mass is not None:
            upper_limit = max_mass

        pn = deconvoluted_peak_set.get_size()
        for pi in range(pn):
            peak = <DeconvolutedPeak>deconvoluted_peak_set.getitem(pi)
            if peak.neutral_mass > upper_limit:
                break
            query_mass = component.mass - peak.neutral_mass
            if (query_mass - upper_limit) / upper_limit > self.product_error_tolerance:
                break
            complements = deconvoluted_peak_set.all_peaks_for(
                query_mass, self.product_error_tolerance)
            on = PyTuple_Size(complements)
            for oi in range(on):
                other_peak = <DeconvolutedPeak>PyTuple_GetItem(complements, oi)
                graph.add(peak, other_peak, component.obj)

    cpdef init_paths(self, graph, limit=200):
        cdef:
            double min_start_mass
            list paths, edge_list
            Path path

        paths = []
        min_start_mass = (<MassWrapper>PyList_GetItem(self.components, 0)).mass + 1
        for edge_list in graph.longest_paths(limit=limit):
            path = Path(edge_list)
            if path.start_mass < min_start_mass:
                continue
            paths.append(path)
        return paths

    def collect_paths(self, paths):
        """Group together paths which share the same annotation sequence and approximate
        start and end masses.

        Parameters
        ----------
        paths: :class:`list`
            A list of :class:`Path` objects to group

        Returns
        -------
        :class:`list` of :class:`list` of :class:`Path` objects
        """
        groups = defaultdict(list)
        if not paths:
            return []
        for path in paths:
            key = tuple(e.annotation for e in path)
            groups[key].append(path)
        result = []
        for key, block_members in groups.items():
            block_members = sorted(block_members, key=lambda x: x.start_mass)
            current_path = block_members[0]
            members = [current_path]
            for path in block_members[1:]:
                if abs(current_path.start_mass - path.start_mass) / path.start_mass < self.product_error_tolerance:
                    members.append(path)
                else:
                    result.append(members)
                    current_path = path
                    members = [current_path]
            result.append(members)
        return result

    def merge_paths(self, paths):
        path_groups = self.collect_paths(paths)
        merged_paths = [collect_edges(group) for group in path_groups]
        merged_paths.sort(key=lambda x: x.total_signal, reverse=True)
        return merged_paths

    def paths(self, scan, limit=200, merge=False, max_mass=None):
        graph = self.find_edges(scan, max_mass=max_mass)
        paths = self.init_paths(graph, limit)
        if merge:
            paths = self.merge_paths(paths)
        return paths
