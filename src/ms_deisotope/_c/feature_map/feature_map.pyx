# cython: embedsignature =True

cimport cython

from ms_peak_picker._c.peak_set cimport PeakBase
from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature, LCMSFeatureTreeNode
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_GetSlice, PyList_GetItem

cimport numpy as np
import numpy as np

np.import_array()

from cpython cimport array
import array



cdef class LCMSFeatureMap(object):
    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.mz)

    def __len__(self):
        return PyList_GET_SIZE(self.features)

    def __iter__(self):
        return iter(self.features)

    def __getitem__(self, i):
        if isinstance(i, (int, slice)):
            return self.features[i]
        else:
            return [self.features[j] for j in i]

    def __eq__(self, other):
        if other is None:
            return False
        return self.features == other.features

    def __ne__(self, other):
        return not self == other

    cdef Py_ssize_t get_size(self):
        return PyList_GET_SIZE(self.features)

    cdef LCMSFeature get(self, size_t i):
        return <LCMSFeature>PyList_GetItem(self.features, i)

    cpdef LCMSFeature search(self, double mz, double error_tolerance=2e-5):
        """Search for a single feature within `error_tolerance` of `mz`.

        Parameters
        ----------
        mz : float
            The m/z to search for
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM

        Returns
        -------
        :class:`~.LCMSFeature`
            The a feature matching the query m/z, or :const:`None`.

        See Also
        --------
        find_all
        between
        """

        cdef:
            long i
            LCMSFeature out
        i = binary_search(self.features, mz, error_tolerance)
        if i == -1:
            return None
        else:
            out = <LCMSFeature>PyList_GET_ITEM(self.features, i)
            return out

    cpdef list find_all(self, double mz, double error_tolerance=2e-5):
        """Search for all features within `error_tolerance` of `mz`.

        Parameters
        ----------
        mz : float
            The m/z to search for
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM

        Returns
        -------
        list of :class:`~.LCMSFeature`
            A list of all features matching the query m/z, or the empty list.

        See Also
        --------
        between
        """
        cdef:
            long loout, hiout
            list result
        search_sweep(self.features, mz, error_tolerance, &loout, &hiout)
        result = <list>PyList_GetSlice(self.features, loout, hiout)
        return result

    def index_range(self, lo, hi, error_tolerance=2e-5):
        return (
            binary_search_with_flag(
                self.features, lo, error_tolerance)[0][0],
            binary_search_with_flag(
                self.features, hi, error_tolerance)[0][0])

    cpdef list spanning_time(self, double time_point):
        cdef:
            size_t i, n
            LCMSFeature feature
            list result
        result = []
        n = self.get_size()
        for i in range(n):
            feature = self.get(i)
            if feature.spans_in_time(time_point):
                result.append(feature)
        return result

    cpdef list between(self, double lo, double hi, double error_tolerance=2e-5):
        """Search for all features between `lo` and `hi`, allowing `error_tolerance`
        around the edges.

        Parameters
        ----------
        lo : float
            The m/z to search for the lower bound
        hi : float
            The m/z to search for the upper bound
        error_tolerance : float, optional
            The error tolerance to search with, defaults to 2e-5, 20 PPM

        Returns
        -------
        list of :class:`~.LCMSFeature`
            The features between the query bounds with the permitted error, or an empty list.

        See Also
        --------
        find_all
        """
        cdef:
            Py_ssize_t n, i, lo_ix, hi_ix
            LCMSFeature f

        n = self.get_size()
        if n == 0:
            return []
        lo_ix = binary_search_with_flag(
            self.features, lo, error_tolerance)[0][0]
        if self.get(lo_ix).get_mz() < lo:
            lo_ix += 1
        if lo_ix > n:
            lo_ix = n - 1
        i = lo_ix
        while i < n:
            f = self.get(i)
            if f.get_mz() > hi:
                break
            i += 1
        hi_ix = i
        return self[lo_ix:hi_ix]

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))

    cpdef LCMSFeatureMap clone(self, bint deep=True):
        cdef:
            size_t i, n
            list result
        result = []
        n = self.get_size()
        for i in range(n):
            feature = self.get(i).clone(deep)
            result.append(feature)
        return self.__class__(result)

    def as_arrays(self):
        cdef:
            size_t i, n, j, m, k, q
            double time
            LCMSFeature feature
            LCMSFeatureTreeNode node
            PeakBase peak

        mz_array = array.array('d')
        intensity_array = array.array('d')
        ion_mobility_array = array.array('d')
        feature_id_array = array.array('L')
        n = self.get_size()
        for i in range(n):
            feature = <LCMSFeature>self.get(i)
            m = feature.get_size()
            for j in range(m):
                node = <LCMSFeatureTreeNode>feature.getitem(j)
                time = node.time
                q = node.get_members_size()
                for k in range(q):
                    peak = <PeakBase>node.getitem(k)
                    ion_mobility_array.append(time)
                    mz_array.append(peak.mz)
                    intensity_array.append(peak.intensity)
                    feature_id_array.append(i)

        mz_array = np.array(mz_array, copy=False)
        intensity_array = np.array(intensity_array, copy=False)
        ion_mobility_array = np.array(ion_mobility_array, copy=False)
        feature_id_array = np.array(feature_id_array, copy=False)
        mask = np.lexsort(np.stack((ion_mobility_array, mz_array)))
        return (mz_array[mask], intensity_array[mask], ion_mobility_array[mask], feature_id_array[mask])


@cython.binding(True)
cpdef list split_sparse(LCMSFeatureMap self, double delta_rt=1.0, size_t min_size=2):
    cdef:
        list result, chunks
        size_t i, n, j, m
        LCMSFeature feature, chunk
    result = []
    n = self.get_size()
    for i in range(n):
        feature = self.get(i)
        chunks = feature.split_sparse(delta_rt)
        m = PyList_GET_SIZE(chunks)
        for j in range(m):
            chunk = <LCMSFeature>PyList_GET_ITEM(chunks, j)
            if chunk.get_size() >= min_size:
                result.append(chunk)
    return result


cpdef tuple binary_search_with_flag(list array, double mz, double error_tolerance):
    cdef:
        int i, n, lo, hi, mid, low_end, high_end
        double err
        LCMSFeature x
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
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
    return [0], False


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
            best_error = abs(err)
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


