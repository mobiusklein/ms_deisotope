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


cdef class TheoreticalIsotopicPattern(object):
    cdef:
        public list base_tid
        public list truncated_tid

    @staticmethod
    cdef TheoreticalIsotopicPattern _create(list base_tid, list truncated_tid=*)

    cpdef TheoreticalIsotopicPattern clone(self)

    cpdef bint _eq(self, object other)

    cpdef TheoreticalIsotopicPattern ignore_below(self, double ignore_below=*)
    cpdef TheoreticalIsotopicPattern truncate_after(self, double truncate_after=*)
    cpdef TheoreticalIsotopicPattern shift(self, double mz, bint truncated=*)
    cpdef TheoreticalIsotopicPattern scale(self, list experimental_distribution, str method=*)

    cdef inline TheoreticalPeak get(self, ssize_t i)
    cdef inline TheoreticalPeak get_base(self, ssize_t i)

    cdef size_t get_size(self)
    cdef size_t get_base_size(self)

    cdef double get_monoisotopic_mz(self)
    cdef list get_processed_peaks(self)



cdef class AveragineCache(object):
    cdef:
        public dict backend
        public Averagine averagine
        public double cache_truncation
        public bint enabled
    
    cdef TheoreticalIsotopicPattern has_mz_charge_pair(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*, double ignore_below=*)
    cpdef TheoreticalIsotopicPattern isotopic_cluster(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*, double ignore_below=*)


cpdef double isotopic_shift(int charge=*)