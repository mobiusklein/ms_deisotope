# cython: embedsignature =True

cimport cython

from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_GetSlice


cdef class LCMSFeatureMap(object):
    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.mz)

    def __len__(self):
        return len(self.features)

    def __iter__(self):
        return iter(self.features)

    def __getitem__(self, i):
        if isinstance(i, (int, slice)):
            return self.features[i]
        else:
            return [self.features[j] for j in i]

    cpdef LCMSFeature _search(self, double mz, double error_tolerance):
        cdef:
            long i
            LCMSFeature out
        i = binary_search(self.features, mz, error_tolerance)
        if i == -1:
            return None
        else:
            out = <LCMSFeature>PyList_GET_ITEM(self.features, i)
            return out

    def search(self, mz, error_tolerance=2e-5):
        return self._search(mz, error_tolerance)

    cpdef list _find_all(self, double mz, double error_tolerance):
        cdef:
            long loout, hiout
            list result
        search_sweep(self.features, mz, error_tolerance, &loout, &hiout)
        result = <list>PyList_GetSlice(self.features, loout, hiout)
        return result

    def find_all(self, mz, error_tolerance=2e-5):
        return self._find_all(mz, error_tolerance)

    def index_range(self, lo, hi, error_tolerance=2e-5):
        return (
            binary_search_with_flag(
                self.features, lo, error_tolerance)[0][0],
            binary_search_with_flag(
                self.features, hi, error_tolerance)[0][0])

    def between(self, lo, hi, error_tolerance=2e-5):
        lo_ix = binary_search_with_flag(
            self.features, lo, error_tolerance)[0][0]
        hi_ix = binary_search_with_flag(
            self.features, hi, error_tolerance)[0][0]
        if self[lo_ix].mz < lo:
            lo_ix += 1
        if self[hi_ix] > hi:
            hi_ix -= 1
        if hi_ix < 0:
            hi_ix = 0
        if lo_ix > len(self):
            lo_ix = len(self) - 1
        return self[lo_ix:hi_ix]

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))


cpdef tuple binary_search_with_flag(list array, double mz, double error_tolerance):
    cdef:
        int n, lo, hi, mid, low_end, high_end
        double err
        LCMSFeature x
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) / 2
        x = <LCMSFeature>PyList_GET_ITEM(array, mid)
        err = (x.get_mz() - mz) / mz
        if abs(err) <= error_tolerance:
            i = mid - 1
            # Begin Sweep forward
            while i > 0:
                x = array[i]
                err = (x.get_mz() - mz) / mz
                if abs(err) <= error_tolerance:
                    i -= 1
                    continue
                else:
                    break
            low_end = i
            i = mid + 1

            # Begin Sweep backward
            while i < n:
                x = array[i]
                err = (x.get_mz() - mz) / mz
                if abs(err) <= error_tolerance:
                    i += 1
                    continue
                else:
                    break
            high_end = i
            return list(range(low_end, high_end)), True
        elif (hi - lo) == 1:
            return [mid], False
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0, False


@cython.cdivision(True)
cdef long binary_search(list array, double mz, double error_tolerance):
    cdef:
        int lo, hi, n, mid, i
        double err
        LCMSFeature x
    lo = 0
    n = hi = PyList_GET_SIZE(array)
    while hi != lo:
        mid = (hi + lo) / 2
        x = <LCMSFeature>PyList_GET_ITEM(array, mid)
        err = (x.get_mz() - mz) / mz
        if abs(err) <= error_tolerance:
            best_index = mid
            best_error = err
            i = mid - 1
            while i > 0:
                x = <LCMSFeature>PyList_GET_ITEM(array, i)
                err = abs((x.get_mz() - mz) / mz)
                if err < best_error:
                    best_error = err
                    best_index = i
                i -= 1

            i = mid + 1
            while i < n:
                x = <LCMSFeature>PyList_GET_ITEM(array, i)
                err = abs((x.get_mz() - mz) / mz)
                if err < best_error:
                    best_error = err
                    best_index = i
                i += 1
            return best_index
        elif (hi - lo) == 1:
            return -1
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


@cython.cdivision(True)
cdef void search_sweep(list array, double mz, double error_tolerance, long* loout, long* hiout):
    cdef:
        int i, n, lo, hi, mid, start, end
        double err
        LCMSFeature x
    lo = 0
    n = hi = PyList_GET_SIZE(array)
    loout[0] = 0
    hiout[0] = 0
    while hi != lo:
        mid = (hi + lo) / 2
        x = <LCMSFeature>PyList_GET_ITEM(array, mid)
        err = (x.get_mz() - mz) / mz
        if abs(err) <= error_tolerance:
            i = mid - 1
            while i > 0:
                x = <LCMSFeature>PyList_GET_ITEM(array, i)
                err = abs((x.get_mz() - mz) / mz)
                if err > error_tolerance:
                    break
                i -= 1
            start = i + 1
            i = mid + 1
            while i < n:
                x = <LCMSFeature>PyList_GET_ITEM(array, i)
                err = abs((x.get_mz() - mz) / mz)
                if err > error_tolerance:
                    break
                i += 1
            end = i
            loout[0] = start
            hiout[0] = end
            return
        elif (hi - lo) == 1:
            return
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return


