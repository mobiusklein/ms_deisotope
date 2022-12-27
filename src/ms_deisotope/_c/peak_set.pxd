from ms_deisotope._c.scoring cimport IsotopicFitRecord
from ms_peak_picker._c.peak_set cimport PeakBase


cdef class _Index:
    cdef:
        public size_t neutral_mass
        public size_t mz

    @staticmethod
    cdef _Index _create(size_t neutral_mass, size_t mz)


cdef class EnvelopePair:
    cdef:
        public double mz
        public double intensity

    cpdef bint _eq(self, EnvelopePair other)

    @staticmethod
    cdef EnvelopePair _create(double mz, double intensity)


cdef class Envelope:
    cdef public tuple pairs

    cdef inline size_t get_size(self)

    cdef inline EnvelopePair getitem(self, size_t i)

    cpdef bint _eq(self, Envelope other)

    cpdef Envelope clone(self)

    @staticmethod
    cdef Envelope _create(tuple pairs)


cdef class DeconvolutedPeak(PeakBase):
    cdef:
        public double neutral_mass
        # public double intensity
        public double signal_to_noise
        public _Index _index
        public double full_width_at_half_max
        public int charge
        public double a_to_a2_ratio
        public double most_abundant_mass
        public double average_mass
        public double score
        # public double area
        public Envelope envelope
        # public double mz
        public IsotopicFitRecord fit
        public bint chosen_for_msms

    cpdef bint _eq(self, DeconvolutedPeak other)

    @staticmethod
    cdef DeconvolutedPeak _create_simple(
        double neutral_mass, double intensity, int charge, double score,
        double mz, Envelope envelope)


cdef class IonMobilityDeconvolutedPeak(DeconvolutedPeak):
    cdef:
        public double drift_time


cdef class DeconvolutedPeakSolution(DeconvolutedPeak):
    cdef:
        public object solution

cdef class DeconvolutedPeakSet(object):
    cdef:
        public tuple peaks
        public tuple _mz_ordered
        public bint indexed


    cpdef reindex(self)

    cpdef DeconvolutedPeakSet clone(self)
    cpdef DeconvolutedPeakSet copy(self)

    cdef DeconvolutedPeak _has_peak(self, double neutral_mass, double error_tolerance=*, bint use_mz=*)

    cpdef DeconvolutedPeak has_peak(self, double neutral_mass, double error_tolerance=*, bint use_mz=*)

    cpdef tuple all_peaks_for(self, double neutral_mass, double tolerance=*)

    cdef DeconvolutedPeak _get_nearest_peak(self, double neutral_mass, double* errout)

    cdef size_t get_size(self)

    cdef DeconvolutedPeak getitem(self, size_t i)
    cdef tuple getslice(self, size_t start, size_t end)


cdef int _binary_search_interval(double* array, double target, double error_tolerance, size_t n, size_t* start, size_t* end) nogil


cdef class DeconvolutedPeakSetIndexed(DeconvolutedPeakSet):
    cdef:
        double* neutral_mass_array
        double* mz_array
        index_list* interval_index

        size_t _size

    cdef void _build_index_arrays(self)
    cdef int _interval_for(self, double neutral_mass, double tolerance, size_t* start, size_t* end) nogil


cdef size_t INTERVAL_INDEX_SIZE

cdef struct index_cell:
    double center_value
    size_t start
    size_t end


cdef struct index_list:
    index_cell* index
    size_t size
    double low
    double high


cdef struct deconvoluted_peak_t:
    double neutral_mass
    double intensity
    int charge
    unsigned int index


cdef enum CPeaksFlags:
    borrowing = 0
    owns_peaks = 1
    owns_index = 2
    owns_peaks_and_index = 3

cdef struct deconvoluted_peak_set_t:
    deconvoluted_peak_t* peaks
    double* mass_index
    size_t size
    CPeaksFlags flags


cdef int create_deconvoluted_peak_set_t(DeconvolutedPeakSetIndexed peak_set, deconvoluted_peak_set_t* destination) except 1
cdef int free_deconvoluted_peak_set_t(deconvoluted_peak_set_t* destination) nogil

cdef deconvoluted_peak_set_t deconvoluted_peak_set_all_peaks_for(deconvoluted_peak_set_t* self, double neutral_mass, double error_tolerance=*) nogil
cdef deconvoluted_peak_t* deconvoluted_peak_set_has_peak(deconvoluted_peak_set_t* self, double neutral_mass, double error_tolerance=*) nogil

cdef int deconvoluted_peak_eq(deconvoluted_peak_t* self, deconvoluted_peak_t* other) nogil
cdef size_t deconvoluted_peak_hash(deconvoluted_peak_t* self) nogil


cdef class _CPeakSet:
    cdef:
        deconvoluted_peak_set_t* ptr

    cdef deconvoluted_peak_t* getitem(self, size_t i) nogil

    cdef deconvoluted_peak_t* _has_peak(self, double neutral_mass, double error_tolerance=*) nogil
    cdef deconvoluted_peak_set_t _all_peaks_for(self, double neutral_mass, double error_tolerance=*) nogil

    cpdef object has_peak(self, double neutral_mass, double error_tolerance=*)
    cpdef tuple all_peaks_for(self, double neutral_mass, double error_tolerance=*)
