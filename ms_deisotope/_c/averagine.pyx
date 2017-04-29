# cython: embedsignature=True

cimport cython
from cpython cimport PyObject
from cpython.float cimport PyFloat_AsDouble
from cpython.list cimport PyList_New, PyList_GET_ITEM, PyList_SET_ITEM, PyList_GET_SIZE, PyList_Append
from cpython.dict cimport PyDict_Next, PyDict_SetItem, PyDict_GetItem

from libc.math cimport floor
from libc.stdlib cimport malloc, free

from brainpy import PROTON as _PROTON, isotopic_variants, calculate_mass as _py_calculate_mass
from brainpy._c.isotopic_distribution cimport _isotopic_variants
from brainpy._c.isotopic_distribution cimport TheoreticalPeak
from brainpy._speedup cimport calculate_mass


cdef double PROTON
PROTON = _PROTON


cdef inline double _round(double x):
    return floor(x + 0.5)


cdef double neutral_mass(double mz,  int z, double charge_carrier=PROTON):
    return (mz * abs(z)) - (z * charge_carrier)


@cython.cdivision
cdef double mass_charge_ratio(double neutral_mass, int z, double charge_carrier=PROTON):
    return (neutral_mass + (z * charge_carrier)) / abs(z)


cdef void slide(double mz, list peaklist):
    cdef:
        size_t n
        size_t i
        TheoreticalPeak peak, first_peak
        double delta

    n = len(peaklist) - 1

    first_peak = <TheoreticalPeak>PyList_GET_ITEM(peaklist, 0) 
    for i in range(n):
        peak = <TheoreticalPeak>PyList_GET_ITEM(peaklist, i + 1)
        delta = peak.mz - first_peak.mz
        peak.mz = mz + delta
    first_peak.mz = mz
    return


cdef dict scale_dict(dict data, double factor):
    cdef:
        dict scaled
        Py_ssize_t i
        PyObject* key
        PyObject* value

    i = 0
    scaled = {}
    while PyDict_Next(data, &i, &key, &value):
        PyDict_SetItem(scaled, <object>key, _round(PyFloat_AsDouble(<object>value) * factor))
    return scaled


cdef class Averagine(object):

    def __init__(self, object base_composition):
        self.base_composition = dict(**base_composition)
        self.base_mass = calculate_mass(self.base_composition)
        assert self.base_mass > 0

    def __reduce__(self):
        return self.__class__, (self.base_composition,), self.__getstate__()

    def __getstate__(self):
        return self.base_composition, self.base_mass

    def __setstate__(self, state):
        composition, mass = state
        self.base_composition = dict(composition)
        self.base_mass = float(mass)

    @cython.cdivision
    cpdef dict scale(self, double mz, int charge=1, double charge_carrier=PROTON):
        cdef:
            double neutral, scale, scaled_mass
            dict scaled
            Py_ssize_t i
            int delta_hydrogen
            int H
            PyObject* key
            PyObject* value
            double count

        neutral = neutral_mass(mz, charge, charge_carrier)

        scale = neutral / self.base_mass
        scaled = scale_dict(self.base_composition, scale)

        delta_hydrogen = <int>_round(calculate_mass(scaled) - neutral)
        H = <int>PyFloat_AsDouble(<object>PyDict_GetItem(scaled, "H"))
        if H > delta_hydrogen:
            PyDict_SetItem(scaled, "H", H - delta_hydrogen)
        else:
            PyDict_SetItem(scaled, "H", 0)

        return scaled

    cdef list _isotopic_cluster(self, double mz, int charge=1, double charge_carrier=PROTON, double truncate_after=0.95,
                                double ignore_below=0.0):
        cdef:
            dict composition
            list result, tid, kept_tid
            double cumsum, base_mz, total
            TheoreticalPeak peak
            size_t i, n
        composition = self.scale(mz, charge, charge_carrier)
        cumsum = 0
        result = []
        tid = _isotopic_variants(composition, npeaks=None, charge=charge, charge_carrier=PROTON)
        n = PyList_GET_SIZE(tid)
        base_mz = (<TheoreticalPeak>PyList_GET_ITEM(tid, 0)).mz
        for i in range(n):
            peak = <TheoreticalPeak>PyList_GET_ITEM(tid, i)
            cumsum += peak.intensity
            result.append(peak)
            peak.mz = (peak.mz - base_mz) + mz
            if cumsum >= truncate_after:
                break
        # Renormalize so the truncated peak list sums to 1.0
        if truncate_after < 1.0:
            n = i + 1
            for i in range(n):
                peak = <TheoreticalPeak>PyList_GET_ITEM(tid, i)
                peak.intensity *= 1. / cumsum

        if ignore_below > 0:
            total = 0.0
            kept_tid = []
            for i in range(n):
                peak = <TheoreticalPeak>PyList_GET_ITEM(result, i)
                if (peak.intensity < ignore_below) and (i > 1):
                    pass
                else:
                    PyList_Append(kept_tid, peak)
                    total += peak.intensity
            n = PyList_GET_SIZE(kept_tid)
            for i in range(n):
                peak = <TheoreticalPeak>PyList_GET_ITEM(kept_tid, i)
                peak.intensity /= total
            result = kept_tid

        return result

    cpdef list isotopic_cluster(self, double mz, int charge=1, double charge_carrier=PROTON, double truncate_after=0.95,
                                double ignore_below=0.0):
        cdef:
            list out
        out = self._isotopic_cluster(mz, charge, charge_carrier, truncate_after, ignore_below)
        return out

    def __getitem__(self, key):
        return self.base_composition[key]

    def __iter__(self):
        return iter(self.base_composition)

    def keys(self):
        return self.base_composition.keys()

    def values(self):
        return self.base_composition.values()

    def items(self):
        return self.base_composition.items()

    def __repr__(self):
        return "Averagine(%r)" % self.base_composition

    def __richcmp__(self, other, int code):
        if code == 2:
            return self.base_composition == other.base_composition
        elif code == 3:
            return self.base_composition != other.base_composition

    def __hash__(self):
        return hash(self.base_composition.items())


cdef list clone_peak_list(list peaklist):
    cdef:
        size_t i
        list result
        TheoreticalPeak peak

    result = []
    for i in range(PyList_GET_SIZE(peaklist)):
        peak = <TheoreticalPeak>PyList_GET_ITEM(peaklist, i)
        result.append(peak.clone())
    return result



cdef class AveragineCache(object):

    def __init__(self, object averagine, object backend=None, double cache_truncation=1.):
        if backend is None:
            backend = {}
        self.backend = dict(backend)
        if isinstance(averagine, AveragineCache):
            self.averagine = averagine.averagine
            self.cache_truncation = averagine.cache_truncation
        else:
            self.averagine = Averagine(averagine)
        self.cache_truncation = cache_truncation
        self.enabled = True

    def __reduce__(self):
        return self.__class__, self.__getstate__()

    def __getstate__(self):
        return self.averagine, self.backend, self.cache_truncation

    def __setstate__(self, state):
        avg, store, trunc = state
        self.averagine = Averagine(avg)
        self.store = dict(store)
        self.cache_truncation = trunc

    cdef list has_mz_charge_pair(self, double mz, int charge=1, double charge_carrier=PROTON, double truncate_after=0.95,
                                 double ignore_below=0.0):
        cdef:
            double key_mz
            tuple key_tuple
            PyObject* pvalue
            list tid
        if self.enabled:
            if self.cache_truncation == 0.0:
                key_mz = mz
            else:
                key_mz = _round(mz / self.cache_truncation) * self.cache_truncation
            key_tuple = (key_mz, charge, charge_carrier, truncate_after)
            pvalue = PyDict_GetItem(self.backend, key_tuple)
            if pvalue == NULL:
                tid = self.averagine._isotopic_cluster(mz, charge, charge_carrier, truncate_after)
                PyDict_SetItem(self.backend, key_tuple, clone_peak_list(tid))
                return tid
            else:
                tid = <list>pvalue
                tid = clone_peak_list(tid)
                slide(mz, tid)
                return tid
        else:
            tid = self.averagine._isotopic_cluster(mz, charge, charge_carrier, truncate_after)
            return tid

    cpdef list isotopic_cluster(self, double mz, int charge=1, double charge_carrier=PROTON, double truncate_after=0.95,
                                double ignore_below=0.0):
        return self.has_mz_charge_pair(mz, charge, charge_carrier, truncate_after, ignore_below)

    def __getitem__(self, key):
        return self.averagine[key]

    def __iter__(self):
        return iter(self.averagine)

    def keys(self):
        return self.averagine.keys()

    def values(self):
        return self.averagine.values()

    def items(self):
        return self.averagine.items()

    def __repr__(self):
        return "AveragineCache(%r)" % self.averagine

    def clear(self):
        self.backend.clear()


cdef double _neutron_shift
_neutron_shift = _py_calculate_mass({"C[13]": 1}) - _py_calculate_mass({"C[12]": 1})


@cython.cdivision
cpdef double isotopic_shift(int charge=1):
    return _neutron_shift / <double>(charge)

