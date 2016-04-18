
cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=*)
cdef double neutral_mass(double mz,  int z, double charge_carrier=*)
cdef void slide(double mz, list peaklist)
cdef dict scale_dict(dict data, double factor)

cdef class Averagine(object):
    cdef:
        public double base_mass
        public dict base_composition
    
    cpdef dict scale(self, double mz, int charge=*, double charge_carrier=*)
    cpdef list isotopic_cluster(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*)
    cdef list _isotopic_cluster(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*)


cdef class AveragineCache(object):
    cdef:
        public dict backend
        public Averagine averagine
    
    cdef list has_mz_charge_pair(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*)
    cpdef list isotopic_cluster(self, double mz, int charge=*, double charge_carrier=*, double truncate_after=*)


cpdef double isotopic_shift(int charge=*)