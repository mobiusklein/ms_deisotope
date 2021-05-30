from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin, IntervalTreeNode

from ms_deisotope._c.feature_map.lcms_feature cimport FeatureBase, EmptyFeature
from ms_deisotope._c.feature_map.feature_fit cimport LCMSFeatureSetFit, map_coord


cdef class SpanningMixin2D(object):
    cdef:
        public map_coord start
        public map_coord end

    cdef bint _contains(self, map_coord i)

    cpdef bint contains(self, map_coord i)

    cpdef bint overlaps(self, SpanningMixin2D interval)

    cpdef bint is_contained_in_interval(self, SpanningMixin2D interval)

    cpdef bint contains_interval(self, SpanningMixin2D interval)


cdef class FeatureNode(object):
    cdef:
        public FeatureBase feature
        public dict links
        Py_hash_t _hash


cdef class DependenceCluster(SpanningMixin2D):
    cdef:
        public object parent
        public list dependencies
        public LCMSFeatureSetFit best_fit
        public bint maximize

    cpdef _reset(self)
    cpdef add(self, LCMSFeatureSetFit fit)
    cpdef disjoint_subset(self, int max_size=*)
    cpdef list disjoint_best_fits(self, int max_size=*)

    cdef LCMSFeatureSetFit _best_fit(self)
    cdef map_coord _start(self)
    cdef map_coord _end(self)


cdef class FeatureSetFitNode(SpanningMixin2D):
    cdef:
        public LCMSFeatureSetFit fit
        public set edges
        public set overlap_edges
        public set feature_indices
        public Py_hash_t _hash
        public long index
        public double score

    cdef void _init_indices(self)

    cpdef visit(self, FeatureSetFitNode other)
    cpdef bint isdisjoint(self, FeatureSetFitNode other)


cdef class ConnectedSubgraph(object):
    cdef:
        public tuple nodes
        public bint maximize

    cpdef populate_edges(self)
    cpdef find_heaviest_path(self, method=*)
    cdef void _init(self)


cdef class FeatureDependenceGraphBase(object):
    cdef:
        public dict nodes
        public object feature_map
        public set dependencies
        public list clusters

        public int max_missing_features
        public bint use_monoisotopic_superceded_filtering
        public bint maximize

    cpdef list nodes_for(self, LCMSFeatureSetFit fit_record, dict cache=*)
    cpdef add_fit_dependence(self, LCMSFeatureSetFit fit_record)
    cpdef drop_superceded_fits(self)
    cpdef drop_fit_dependence(self, LCMSFeatureSetFit fit_record, dict cache=*)
    cpdef best_exact_fits(self)
    cpdef drop_gapped_fits(self, n=*)
    cpdef _gather_independent_clusters(self, dict nodes_for_cache=*)


