cimport cython

from ms_deisotope._c.peak_set cimport DeconvolutedPeak


cdef class NodeBase(object):

    cdef:
        public size_t key

    cdef double get_neutral_mass(self)
    cdef double get_intensity(self)
    cdef size_t get_index(self)
    cdef bint _eq(self, NodeBase other)


cdef class PlaceHolderNode(NodeBase):
    cdef:
        public str label
        public double _mass
        public double _intensity
        public size_t _index


cdef class PeakNode(NodeBase):
    cdef:
        public DeconvolutedPeak peak



cdef class PeakGroupNode(NodeBase):
    cdef:
        public list peaks

    cdef size_t get_size(self)


@cython.final
@cython.freelist(1000000)
cdef class Edge(object):
    cdef:
        public NodeBase start
        public NodeBase end
        public object annotation
        public (size_t, size_t) key
        public Py_hash_t _hash

    @staticmethod
    cdef Edge _create(NodeBase start, NodeBase end, object annotation)


cdef class SpectrumGraph(object):
    cdef:
        public set transitions
        public object by_first
        public object by_second

    cpdef add(self, DeconvolutedPeak p1, DeconvolutedPeak p2, object annotation)
    cpdef list topological_sort(self, adjacency_matrix=*)

    cpdef list paths_starting_at(self, ix, int limit=*)
    cpdef list paths_ending_at(self, ix, int limit=*)

    cpdef list transitive_closure(self, list paths, int limit=*)


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
    cdef:
        public object obj
        public double mass


cdef class Path(object):
    cdef:
        public list transitions
        public double total_signal
        public double start_mass
        public double end_mass
        public set _peaks_used
        public set _edges_used

    cdef double _total_signal(self)
    cdef Edge get(self, Py_ssize_t i)
    cdef size_t get_size(self)