from ms_peak_picker._c.peak_set cimport FittedPeak


cdef class LCMSFeatureTreeList(object):
    cdef:
        public list roots
        public set _node_id_hash

    cdef void _invalidate(self)

    cpdef tuple find_time(self, double retention_time)
    cdef LCMSFeatureTreeNode _find_time(self, double retention_time, size_t* indexout)
    cdef LCMSFeatureTreeNode getitem(self, size_t i)


cdef class LCMSFeatureTreeNode(object):
    cdef:
        public double retention_time
        public list members
        public FittedPeak _most_abundant_member
        public double _mz
        public object node_id

    cpdef double _total_intensity_members(self)
    cpdef double max_intensity(self)
    cpdef double total_intensity(self)

    cpdef bint _eq(self, LCMSFeatureTreeNode other)
    cpdef bint _ne(self, LCMSFeatureTreeNode other)


cdef class FeatureBase(object):
    cdef:
        public LCMSFeatureTreeList nodes
        public double _mz
        public double _start_time
        public double _end_time

    cdef double get_mz(self)
    cdef double get_start_time(self)
    cdef double get_end_time(self)


cdef class LCMSFeature(FeatureBase):
    cdef:
        public double _total_intensity
        public double _last_mz
        public object _most_abundant_member
        public object _last_most_abundant_member
        public object _retention_times
        public object _peaks
        public object _scan_ids
        public object adducts
        public object used_as_adduct
        public object created_at
        public object feature_id

    cpdef bint _eq(self, other)
    cpdef bint _ne(self, other)
    cpdef bint overlaps_in_time(self, FeatureBase interval)

    cdef LCMSFeatureTreeNode getitem(self, size_t i)


cdef class EmptyFeature(FeatureBase):
    pass


cdef class FeatureSetIterator(object):
    cdef:
        public list features
        public list real_features
        public double start_time
        public double end_time
        public object time_point_queue
        public double last_time_seen
        public bint done

    @staticmethod
    cdef FeatureSetIterator _create(list features)

    cpdef double get_next_time(self)
    cpdef list get_peaks_for_time(self, double time)
    cpdef bint has_more(self)
    cpdef list get_next_value(self)
