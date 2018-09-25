from cpython.list cimport PyList_GetItem, PyList_Size
from cpython.tuple cimport PyTuple_GetItem

from ms_deisotope._c.peak_dependency_network.intervals cimport (SpanningMixin, Interval)
from ms_deisotope._c.scoring cimport IsotopicFitRecord


cdef class FitNode(SpanningMixin):

    cdef:
        public size_t index
        public IsotopicFitRecord fit
        public set edges
        public set overlap_edges
        public set peak_indices
        public long _hash
        public double score

    cpdef bint isdisjoint(self, FitNode other)
    cpdef visit(self, FitNode other)

    cdef void _init_fields(self, list experimental)

    @staticmethod
    cdef FitNode _create(IsotopicFitRecord fit, size_t index)


cdef class ConnectedSubgraph(object):
    cdef:
        public tuple nodes
        public bint maximize

    cdef void _init_nodes(self, list fits)
    cdef void populate_edges(self)
    cpdef find_heaviest_path(self, method=*)
