from ms_peak_picker._c.peak_set cimport FittedPeak, PeakBase


cdef class LCMSFeatureTreeList(object):
    cdef:
        public list roots
        public set _node_id_hash

    @staticmethod
    cdef LCMSFeatureTreeList _create(list roots)

    cdef void _invalidate(self)

    cpdef tuple find_time(self, double retention_time)
    cdef inline LCMSFeatureTreeNode _find_time(self, double retention_time, size_t* indexout)
    cdef inline LCMSFeatureTreeNode getitem(self, size_t i)

    cdef inline size_t get_size(self)


cdef class LCMSFeatureTreeNode(object):
    cdef:
        public double retention_time
        public list members
        public PeakBase _most_abundant_member
        public double _mz
        public object node_id

    cpdef double _total_intensity_members(self)
    cpdef double max_intensity(self)
    cpdef double total_intensity(self)

    cdef PeakBase _find_most_abundant_member(self)
    cpdef _calculate_most_abundant_member(self)
    cpdef _recalculate(self)

    cpdef bint _eq(self, LCMSFeatureTreeNode other)
    cpdef bint _ne(self, LCMSFeatureTreeNode other)
    cdef inline FittedPeak getpeak(self, size_t i)
    cdef PeakBase getitem(self, size_t i)
    cdef inline size_t get_members_size(self)

    cdef double get_mz(self)


cdef class FeatureBase(object):
    cdef:
        public LCMSFeatureTreeList nodes
        public double _mz
        public double _start_time
        public double _end_time

    cdef double get_mz(self)
    cdef double get_start_time(self)
    cdef double get_end_time(self)

    cdef inline size_t get_size(self)
    cpdef tuple find_time(self, double retention_time)
    cdef inline LCMSFeatureTreeNode _find_time(self, double retention_time, size_t* indexout)


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
    @staticmethod
    cdef EmptyFeature _create(double mz)


cdef class FeatureSetIterator(object):
    cdef:
        public list features
        public list real_features
        public double start_time
        public double end_time
        public double last_time_seen
        size_t* index_list

        @staticmethod
        cdef FeatureSetIterator _create(list features)

        @staticmethod
        cdef FeatureSetIterator _create_with_threshold(list features, list theoretical_distribution, double detection_threshold)

        cdef void _initialize(self, list features)

        cpdef init_indices(self)
        cpdef double get_next_time(self)
        cpdef bint has_more(self)
        cpdef list get_peaks_for_time(self, double time)
        cpdef list get_next_value(self)

        cdef inline size_t get_size(self)
        cdef inline FeatureBase getitem(self, size_t i)
