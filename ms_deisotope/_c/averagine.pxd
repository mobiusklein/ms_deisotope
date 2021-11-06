cimport cython

from brainpy._c.isotopic_distribution cimport TheoreticalPeak

cdef double PROTON

cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=*)
cdef double neutral_mass(double mz,  int z, double charge_carrier=*)
cdef void slide(double mz, list peaklist)
cdef dict scale_dict(dict data, double factor)

cdef class Averagine(object):
    cdef:
        public double base_mass
        public dict base_composition

    cpdef dict scale(self, double mz, int charge=*, double charge_carrier=*)
    cpdef TheoreticalIsotopicPattern isotopic_cluster(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*, double ignore_below=*)
    cdef TheoreticalIsotopicPattern _isotopic_cluster(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*, double ignore_below=*)


@cython.final
cdef class TheoreticalIsotopicPattern(object):
    cdef:
        public list peaklist
        public double origin
        public double offset

    @staticmethod
    cdef TheoreticalIsotopicPattern _create(list peaklist, double origin, double offset)

    cpdef TheoreticalIsotopicPattern clone(self)

    cpdef bint _eq(self, object other)

    @cython.final
    cdef bint _eq_inst(self, TheoreticalIsotopicPattern other)

    cpdef TheoreticalIsotopicPattern ignore_below(self, double ignore_below=*)
    cpdef TheoreticalIsotopicPattern truncate_after(self, double truncate_after=*)
    cpdef TheoreticalIsotopicPattern shift(self, double mz)
    cpdef TheoreticalIsotopicPattern scale(self, list experimental_distribution, str method=*)
    cpdef TheoreticalIsotopicPattern scale_raw(self, double scale_factor)
    cpdef double drop_last_peak(self)

    @cython.final
    cpdef double total(self)

    @cython.final
    cpdef TheoreticalIsotopicPattern normalize(self)

    @cython.final
    cdef inline TheoreticalIsotopicPattern _scale(self, list experimental_distribution, str method=*)

    @cython.final
    cdef inline TheoreticalIsotopicPattern clone_shift(self, double mz)

    @cython.final
    cdef inline TheoreticalPeak get(self, ssize_t i)

    @cython.final
    cdef inline size_t get_size(self)

    @cython.final
    cdef inline double get_monoisotopic_mz(self)
    @cython.final
    cdef inline list get_processed_peaks(self)

    @cython.final
    cpdef size_t basepeak_index(self)

    @cython.final
    cpdef list incremental_truncation(TheoreticalIsotopicPattern self, double threshold)

    cdef TheoreticalIsotopicPattern clone_drop_last(self)


@cython.final
cdef class AveragineCache(object):
    cdef:
        public dict backend
        public Averagine averagine
        public double cache_truncation
        public bint enabled

    cdef double _make_cache_key(self, double mz)

    cdef TheoreticalIsotopicPattern has_mz_charge_pair(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*, double ignore_below=*)
    cpdef TheoreticalIsotopicPattern isotopic_cluster(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*, double ignore_below=*)


cpdef double isotopic_shift(int charge=*)

cpdef TheoreticalIsotopicPattern poisson_approximate(double mass, size_t n_peaks, double lambda_factor=*, int charge=*)