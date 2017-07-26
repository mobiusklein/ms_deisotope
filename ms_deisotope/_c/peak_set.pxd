from ms_deisotope._c.scoring cimport IsotopicFitRecord
from ms_peak_picker._c.peak_set cimport PeakBase


cdef class _Index:
    cdef:
        public size_t neutral_mass
        public size_t mz


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


cdef class DeconvolutedPeakSolution(DeconvolutedPeak):
    cdef:
        public object solution

cdef class DeconvolutedPeakSet:
    cdef:
        public tuple peaks
        public tuple _mz_ordered

    cdef DeconvolutedPeak _has_peak(self, double neutral_mass, double error_tolerance=*, bint use_mz=*)

    cpdef DeconvolutedPeak has_peak(self, double neutral_mass, double error_tolerance=*, bint use_mz=*)

    cdef DeconvolutedPeak getitem(self, size_t i)
