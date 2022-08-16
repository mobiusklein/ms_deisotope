from ms_peak_picker._c.peak_set cimport FittedPeak, PeakBase


cdef class LCMSFeatureTreeList(object):
    cdef:
        public list roots
        public set _node_id_hash

    @staticmethod
    cdef LCMSFeatureTreeList _create(list roots)

    cdef void _invalidate(self)

    cpdef LCMSFeatureTreeNodeBase _make_node(self, double time, list peaks)

    cpdef tuple find_time(self, double time)
    cdef inline LCMSFeatureTreeNodeBase _find_time(self, double time, size_t* indexout)
    cdef inline LCMSFeatureTreeNodeBase getitem(self, size_t i)

    cdef inline size_t get_size(self)


cdef class LCMSFeatureTreeNodeBase(object):
    cdef:
        public double time
        public list members
        public PeakBase _most_abundant_member
        public object node_id

    cpdef double _total_intensity_members(self)
    cpdef double max_intensity(self)
    cpdef double total_intensity(self)

    cdef PeakBase _find_most_abundant_member(self)
    cpdef _calculate_most_abundant_member(self)
    cpdef _recalculate(self)

    cdef PeakBase getitem(self, size_t i)
    cdef inline size_t get_members_size(self)




cdef class LCMSFeatureTreeNode(LCMSFeatureTreeNodeBase):
    cdef:

        public double _mz

    cpdef bint _eq(self, LCMSFeatureTreeNode other)
    cpdef bint _ne(self, LCMSFeatureTreeNode other)
    cdef inline FittedPeak getpeak(self, size_t i)

    cdef double get_mz(self)


cdef class FeatureBase(object):
    cdef:
        public LCMSFeatureTreeList nodes
        public double _mz
        public double _start_time
        public double _end_time

    cdef double get_mz(self)
    cdef double get_neutral_mass(self)
    cdef double get_start_time(self)
    cdef double get_end_time(self)

    cdef inline size_t get_size(self)
    cpdef tuple find_time(self, double time)
    cdef inline LCMSFeatureTreeNode _find_time(self, double time, size_t* indexout)


cdef class LCMSFeature(FeatureBase):
    cdef:
        public double _total_intensity
        public double _last_mz
        public object _times
        public object _peaks
        public object adducts
        public object used_as_adduct
        public object created_at
        public object feature_id
        public RunningWeightedAverage _peak_averager

    cpdef bint _eq(self, other)
    cpdef bint _ne(self, other)
    cpdef bint overlaps_in_time(self, FeatureBase interval)
    cpdef bint spans_in_time(self, double time)

    cdef LCMSFeatureTreeNode getitem(self, size_t i)

    cpdef _initialize_averager(self)
    cpdef _update_from_averager(self)
    cpdef _reaverage_and_update(self)

    cpdef insert_node(self, LCMSFeatureTreeNode node)
    cpdef insert(self, PeakBase peak, double time)
    cpdef _invalidate(self, bint reaverage=*)
    cpdef list split_sparse(self, double delta_rt=*)
    cpdef LCMSFeature clone(self, deep=*, cls=*)
    cpdef double max_intensity(self)


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

        cdef int _initialize(self, list features) except 1

        cpdef init_indices(self)
        cpdef double get_next_time(self)
        cpdef double get_current_time(self)
        cpdef bint has_more(self)
        cpdef list get_peaks_for_time(self, double time)
        cpdef list get_next_value(self)

        cdef inline size_t get_size(self)
        cdef inline FeatureBase getitem(self, size_t i)


cdef class RunningWeightedAverage(object):
    cdef:
        public double current_mean
        public size_t current_count
        public double total_weight

    cpdef reset(self)
    cpdef _initialize(self)
    cpdef add(self, PeakBase peak)
    cpdef RunningWeightedAverage update(self, iterable)
    cpdef int feed_from_feature(self, LCMSFeature feature)

    @staticmethod
    cdef RunningWeightedAverage _create(list peaks)
    cdef void _update(self, list peaks)


cdef class RunningWeightedAverageNeutralMass(RunningWeightedAverage):

    @staticmethod
    cdef RunningWeightedAverageNeutralMass _create(list peaks)
