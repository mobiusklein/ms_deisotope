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

from ms_peak_picker._c.peak_set cimport FittedPeak


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

    n = PyList_GET_SIZE(peaklist) - 1

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

    cdef TheoreticalIsotopicPattern _isotopic_cluster(self, double mz, int charge=1, double charge_carrier=PROTON, double truncate_after=0.95,
                                double ignore_below=0.0):
        cdef:
            dict composition
            list tid
            TheoreticalIsotopicPattern isotopic_pattern

        composition = self.scale(mz, charge, charge_carrier)
        tid = _isotopic_variants(composition, npeaks=None, charge=charge, charge_carrier=PROTON)

        isotopic_pattern = TheoreticalIsotopicPattern._create(
            tid, (<TheoreticalPeak>PyList_GET_ITEM(tid, 0)).mz, 0)
        isotopic_pattern.shift(mz)

        if truncate_after < 1.0:
            isotopic_pattern.truncate_after(truncate_after)

        if ignore_below > 0:
            isotopic_pattern.ignore_below(ignore_below)

        return isotopic_pattern

    cpdef TheoreticalIsotopicPattern isotopic_cluster(self, double mz, int charge=1, double charge_carrier=PROTON, double truncate_after=0.95,
                                double ignore_below=0.0):
        cdef:
            TheoreticalIsotopicPattern out
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
        return hash(frozenset(self.base_composition.items()))


cdef list clone_peak_list(list peaklist):
    cdef:
        size_t i
        list result
        TheoreticalPeak peak

    result = []
    for i in range(PyList_GET_SIZE(peaklist)):
        result.append((<TheoreticalPeak>PyList_GET_ITEM(peaklist, i)).clone())
    return result


cdef double sum_intensity(list peaklist):
    cdef:
        size_t i
        double total
        FittedPeak peak
    total = 0
    for i in range(PyList_GET_SIZE(peaklist)):
        total += (<FittedPeak>PyList_GET_ITEM(peaklist, i)).intensity
    return total


@cython.freelist(1000000)
cdef class TheoreticalIsotopicPattern(object):

    @staticmethod
    cdef TheoreticalIsotopicPattern _create(list peaklist, double origin, double offset):
        cdef:
            TheoreticalIsotopicPattern self

        self = TheoreticalIsotopicPattern.__new__(TheoreticalIsotopicPattern)
        self.peaklist = peaklist
        self.origin = origin
        self.offset = offset
        return self

    def __init__(self, peaklist, origin, offset=None):
        self.peaklist = list(peaklist)
        self.origin = origin
        if offset is None:
            offset = self.peaklist[0].mz - origin
        self.offset = offset


    def __getitem__(self, i):
        return self.peaklist[i]

    def __iter__(self):
        return iter(self.peaklist)

    def __len__(self):
        return self.get_size()

    @cython.final
    cdef inline TheoreticalPeak get(self, ssize_t i):
        return <TheoreticalPeak>PyList_GET_ITEM(self.peaklist, i)

    @cython.final
    cdef inline size_t get_size(self):
        return PyList_GET_SIZE(self.peaklist)

    cpdef TheoreticalIsotopicPattern clone(self):
        cdef:
            size_t i, n
            TheoreticalPeak p
            list peaklist
        peaklist = []
        n = self.get_size()
        for i in range(n):
            p = self.get(i)
            peaklist.append(TheoreticalPeak._create(p.mz, p.intensity, p.charge))
        return TheoreticalIsotopicPattern._create(peaklist, self.origin, self.offset)

    def __reduce__(self):
        return self.__class__, (self.peaklist, self.origin, self.offset)

    @cython.final
    cdef inline list get_processed_peaks(self):
        return self.peaklist

    cdef inline double get_monoisotopic_mz(self):
        return self.origin

    @property
    def monoisotopic_mz(self):
        return self.get_monoisotopic_mz()

    @cython.cdivision
    cpdef TheoreticalIsotopicPattern ignore_below(self, double ignore_below=0.0):
        cdef:
            double total
            list kept_tid
            size_t i, n
            TheoreticalPeak p
        total = 0
        kept_tid = []
        n = self.get_size()
        for i in range(n):
            p = self.get(i)
            if (p.intensity < ignore_below) and (i > 1):
                continue
            else:
                total += p.intensity
                PyList_Append(kept_tid, p)
        self.peaklist = kept_tid
        self.offset = self.origin - self.get(0).mz
        n = self.get_size()
        for i in range(n):
            p = self.get(i)
            p.intensity /= total
        return self

    cpdef TheoreticalIsotopicPattern shift(self, double mz):
        cdef:
            TheoreticalPeak peak
            size_t i, n
            double delta

        delta = mz - self.origin
        self.origin = mz

        n = self.get_size()
        for i in range(n):
            peak = self.get(i)
            peak.mz += delta
        return self

    cdef TheoreticalIsotopicPattern clone_shift(self, double mz):
        cdef:
            size_t i, n
            TheoreticalPeak p
            list peaklist
            double delta

        delta = mz - self.origin
        peaklist = []
        n = self.get_size()
        for i in range(n):
            p = self.get(i)
            peaklist.append(TheoreticalPeak._create(p.mz + delta, p.intensity, p.charge))
        return TheoreticalIsotopicPattern._create(peaklist, mz, self.offset)


    @cython.cdivision
    cpdef TheoreticalIsotopicPattern truncate_after(self, double truncate_after=0.95):
        cdef:
            double cumsum, normalizer
            TheoreticalPeak peak
            list result
            size_t i, n
        cumsum = 0
        result = []
        n = self.get_size()
        for i in range(n):
            peak = self.get(i)
            cumsum += peak.intensity
            PyList_Append(result, peak)
            if cumsum >= truncate_after:
                break
        self.peaklist = result
        n = self.get_size()
        normalizer = 1. / cumsum
        for i in range(n):
            peak = self.get(i)
            peak.intensity *= normalizer
        return self

    cpdef TheoreticalIsotopicPattern scale(self, list experimental_distribution, str method="sum"):
        return self._scale(experimental_distribution, method)

    @cython.final
    @cython.cdivision
    cdef inline TheoreticalIsotopicPattern _scale(self, list experimental_distribution, str method="sum"):
        cdef:
            size_t i, j, n
            TheoreticalPeak peak
            FittedPeak expeak
            double total_abundance, maximum, scale_factor
            double weights, scales, w

        n = self.get_size()
        if n == 0:
            raise ValueError("Isotopic Pattern has length 0 (%f, %r)" % (self.origin, self.peaklist))
        if method == "sum":
            total_abundance = sum_intensity(experimental_distribution)
            for i in range(n):
                peak = self.get(i)
                peak.intensity *= total_abundance
        elif method == "max":
            i = 0
            maximum = 0
            for j in range(n):
                peak = self.get(j)
                if peak.intensity > maximum:
                    maximum = peak.intensity
                    j = i
            scale_factor = (<FittedPeak>PyList_GET_ITEM(
                experimental_distribution, i)).intensity / maximum
            for j in range(n):
                peak = self.get(j)
                peak.intensity *= scale_factor
        elif method == "meanscale":
            scales = 0
            weights = 0
            total_abundance = 0
            for i in range(n):
                expeak = <FittedPeak>PyList_GET_ITEM(experimental_distribution, i)
                peak = self.get(i)
                total_abundance += expeak.intensity
                w = (peak.intensity * expeak.intensity ** 2)
                scales += expeak.intensity / peak.intensity * w
                weights += w
            scale_factor = scales / weights
            for i in range(n):
                peak = self.get(i)
                peak.intensity *= scale_factor
        return self

    def _scale_raw(self, double scale_factor):
        for peak in self:
            peak.intensity *= scale_factor

    def __repr__(self):
        return "TheoreticalIsotopicPattern(%0.4f, charge=%d, (%s))" % (
            self.monoisotopic_mz,
            self.peaklist[0].charge,
            ', '.join("%0.3f" % p.intensity for p in self.peaklist))

    cpdef bint _eq(self, object other):
        cdef:
            list peaklist
            TheoreticalIsotopicPattern other_typed
        if isinstance(other, list):
            peaklist = other
        elif isinstance(other, TheoreticalIsotopicPattern):
            peaklist = other.peaklist
        else:
            raise TypeError(type(other))
        return self.get_processed_peaks() == peaklist

    def __richcmp__(self, object other, int code):
        if other is None:
            if code == 3:
                return True
            else:
                return False

        if code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))
        else:
            return NotImplemented


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

    @cython.cdivision
    cdef TheoreticalIsotopicPattern has_mz_charge_pair(self, double mz, int charge=1, double charge_carrier=PROTON, double truncate_after=0.95,
                                 double ignore_below=0.0):
        cdef:
            double key_mz
            tuple cache_key
            PyObject* pvalue
            TheoreticalIsotopicPattern tid
        if self.enabled:
            if self.cache_truncation == 0.0:
                key_mz = mz
            else:
                key_mz = _round(mz / self.cache_truncation) * self.cache_truncation

            # Attempting to replace this tuple construction (which in turn necessitates packing each
            # numeric argument as a Python object) with a hand-written extension class that can compute
            # its own hash value without invoking any Python operations turns out to be just a bit slower
            # than the bare tuple itself.
            cache_key = (key_mz, charge, charge_carrier, truncate_after)
            pvalue = PyDict_GetItem(self.backend, cache_key)
            if pvalue == NULL:
                tid = self.averagine._isotopic_cluster(mz, charge, charge_carrier, truncate_after)
                PyDict_SetItem(self.backend, cache_key, tid.clone())
                return tid
            else:
                tid = <TheoreticalIsotopicPattern>pvalue
                tid = tid.clone_shift(mz)
                return tid
        else:
            tid = self.averagine._isotopic_cluster(mz, charge, charge_carrier, truncate_after)
            return tid

    cpdef TheoreticalIsotopicPattern isotopic_cluster(
                                self, double mz, int charge=1, double charge_carrier=PROTON,
                                double truncate_after=0.95, double ignore_below=0.0):
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

