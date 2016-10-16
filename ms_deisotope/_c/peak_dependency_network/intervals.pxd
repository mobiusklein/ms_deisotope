cdef class Interval(SpanningMixin):
    cdef:
        public object members
        public dict data


cdef class SpanningMixin(object):
    cdef:
        public double start
        public double end

    cpdef bint overlaps(self, SpanningMixin interval)

    cpdef double overlap_size(self, SpanningMixin interval)

    cpdef bint is_contained_in_interval(self, SpanningMixin interval)

    cpdef bint contains_interval(self, SpanningMixin interval)


cdef class IntervalTreeNode(object):
    cdef:
        public double center
        public double start
        public double end
        public IntervalTreeNode left
        public IntervalTreeNode right
        public list contained
        public int level
        public IntervalTreeNode parent

    cpdef IntervalTreeNode node_contains_point(self, double x)
    cpdef list contains_point(self, double x)
    cpdef list _overlaps_interval(self, double start, double end)
    cpdef list overlaps(self, double start, double end)
    cpdef bint _eq(self, IntervalTreeNode other)
