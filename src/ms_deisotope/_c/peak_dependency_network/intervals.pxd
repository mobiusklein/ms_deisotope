cdef class SpanningMixin(object):
    cdef:
        public double start
        public double end

    cdef bint _contains(self, double i)

    cpdef bint contains(self, double i)

    cpdef bint overlaps(self, SpanningMixin interval)

    cpdef double overlap_size(self, SpanningMixin interval)

    cpdef bint is_contained_in_interval(self, SpanningMixin interval)

    cpdef bint contains_interval(self, SpanningMixin interval)


cdef class Interval(SpanningMixin):
    cdef:
        public object members
        public dict data

cdef class SimpleInterval(SpanningMixin):
    @staticmethod
    cdef SimpleInterval _create(double start, double end)

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

    cpdef list contains_point(self, double x)
    cpdef list _overlaps_interval(self, double start, double end)
    cpdef list overlaps(self, double start, double end)
    cpdef bint _eq(self, IntervalTreeNode other)


cdef class SpanningMixin2D(object):
    cdef:
        public double[2] _start
        public double[2] _end

    cdef bint _contains(self, double[2] i)
    cpdef bint contains(self, double[:] i)
    cpdef bint overlaps(self, SpanningMixin2D interval)
    cpdef bint contains_interval(self, SpanningMixin2D interval)


cdef class Interval2D(SpanningMixin2D):
    cdef:
        public list members
        public dict data


cdef class IntervalTreeNode2D(IntervalTreeNode):
    cdef:
        public object inner_organizer
        public IntervalTreeNode organized

    cpdef list _overlaps_interval_2d(self, double[::1] starts, double[::1] ends)
    cpdef list overlaps_2d(self, double[::1] starts, double[::1] ends)
