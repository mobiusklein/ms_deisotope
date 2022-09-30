from ms_peak_picker._c.peak_set cimport FittedPeak, PeakSet

from ms_deisotope._c.scoring cimport IsotopicFitRecord
from ms_deisotope._c.peak_dependency_network.intervals cimport SpanningMixin, IntervalTreeNode
from ms_deisotope._c.peak_dependency_network.subgraph cimport ConnectedSubgraph


cdef class PeakNode(object):
    """
    Represent a single FittedPeak and the junction of multiple
    IsotopicFitRecords which depend uponit

    Attributes
    ----------
    links : dict
        Mapping IsotopicFitRecords which depend upon this peak to those fit's scores
    peak : ms_peak_picker.FittedPeak
        The peak being depended upon
    """

    cdef:
        public FittedPeak peak
        public dict links
        long _hash

    @staticmethod
    cdef PeakNode _create(FittedPeak peak)


cdef class DependenceCluster(SpanningMixin):
    """
    Represent a cluster of isotopic fits which are overlapping

    Attributes
    ----------
    dependencies : set
        The set of all fits which are codependent
    start: float
        The smallest mz of all elements of this cluster
    end: float
        The largest mz of all elements of this cluster
    maximize: bool
        Whether the objective is to maximize or minimize the
        isotopic fit score
    """

    cdef:
        public object parent
        public list dependencies
        public bint maximize
        public IsotopicFitRecord best_fit

    cpdef _reset(self)
    cpdef add(self, IsotopicFitRecord fit)
    cpdef disjoint_subset(self)
    cpdef ConnectedSubgraph build_graph(self)
    cdef IsotopicFitRecord _best_fit(self)
    cpdef list disjoint_best_fits(self)
    cdef double _start(self)
    cdef double _end(self)


cdef class PeakDependenceGraphBase(object):
    cdef:
        public dict nodes
        public set dependencies
        public PeakSet peaklist
        public int max_missed_peaks
        public bint use_monoisotopic_superceded_filtering
        public bint maximize

        public list clusters
        public dict _solution_map
        public list _all_clusters

        public IntervalTreeNode _interval_tree

    cpdef reset(self)
    cpdef _populate_initial_graph(self)
    cpdef add_fit_dependence(self, IsotopicFitRecord fit_record)
    cpdef list nodes_for(self, IsotopicFitRecord fit_record, dict cache=*)
    cpdef drop_fit_dependence(self, IsotopicFitRecord fit_record)
    cpdef best_exact_fits(self)
    cpdef _gather_independent_clusters(self, dict nodes_for_cache=*)