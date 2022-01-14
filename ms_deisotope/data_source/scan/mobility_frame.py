import array
import warnings

from abc import abstractmethod
from weakref import WeakValueDictionary
from itertools import chain

from collections import namedtuple, defaultdict

import numpy as np

from ms_peak_picker import average_signal
from ms_deisotope.data_source.scan.scan import WrappedScan

from ms_deisotope.peak_set import IonMobilityDeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.peak_dependency_network import IntervalTreeNode, Interval

from .base import RawDataArrays, ScanBunch

class IonMobilitySource(object):
    @abstractmethod
    def _frame_id(self, data):
        raise NotImplementedError()

    @abstractmethod
    def _frame_index(self, data):
        raise NotImplementedError()

    def _frame_time(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.scan_time

    def _frame_ms_level(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.ms_level

    @abstractmethod
    def _frame_start_scan_index(self, data):
        raise NotImplementedError()

    @abstractmethod
    def _frame_end_scan_index(self, data):
        raise NotImplementedError()

    def _frame_precursor_information(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.precursor_information

    def _frame_activation(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.activation

    def _frame_isolation_window(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.isolation_window

    def _frame_polarity(self, data):
        scan = self.get_scan_by_index(
            self._frame_start_scan_index(data))
        return scan.polarity

    def _frame_drift_times(self, data):
        scans = []
        for i in range(self._frame_start_scan_index(data), self._frame_end_scan_index(data)):
            scan = self.get_scan_by_index(i)
            scans.append(scan.drift_time)
        return scans

    def _frame_scans(self, data):
        scans = []
        for i in range(self._frame_start_scan_index(data), self._frame_end_scan_index(data)):
            scan = self.get_scan_by_index(i)
            scans.append(scan)
        return scans

    def _frame_annotations(self, data):
        return {}

    def _frame_arrays(self, data):
        scans = self._frame_scans(data)
        return RawDataArrays3D.stack(scans)

    def _frame_ion_mobilities(self, data):
        scans = self._frame_scans(data)
        return np.array([s.drift_time for s in scans])

    def _frame_acquisition_information(self, data):
        scan = self.get_scan_by_index(self._frame_start_scan_index(data))
        acq = scan.acquisition_information
        if acq is None:
            return acq
        acq = acq.copy()
        for event in acq:
            im = event._ion_mobility
            mt = im.ion_mobility_type()
            if mt is not None:
                im.remove_ion_mobility_type(mt)
        return acq


class IonMobilitySourceRandomAccessFrameSource(IonMobilitySource):
    @abstractmethod
    def get_frame_by_index(self, index):
        raise NotImplementedError()

    @abstractmethod
    def get_frame_by_time(self, time):
        raise NotImplementedError()

    def _validate_frame(self, data):
        return True

    def _make_frame(self, data):
        return IonMobilityFrame(data, self)

    def _cache_frame(self, frame):
        pass

    @abstractmethod
    def _default_frame_iterator(self, start_index=None):
        raise NotImplementedError()

    @abstractmethod
    def make_frame_iterator(self, iterator=None, grouped=False):
        raise NotImplementedError()

    @abstractmethod
    def start_from_frame(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        raise NotImplementedError()

    def initialize_frame_cache(self):
        '''Initialize a cache which keeps track of which :class:`~.IonMobilityFrame`
        objects are still in memory using a :class:`weakref.WeakValueDictionary`.

        When a frame is requested, if the frame object is found in the cache, the
        existing object is returned rather than re-read from disk.
        '''
        self._frame_cache = WeakValueDictionary()

    @property
    def frame_cache(self):
        '''A :class:`weakref.WeakValueDictionary` mapping used to retrieve
        frames from memory if available before re-reading them from disk.
        '''
        return self._frame_cache

    @frame_cache.setter
    def frame_cache(self, value):
        self._frame_cache = value


class _upto_index(object):
    def __init__(self, values, last_value):
        self.values = values
        self.last_value = last_value
        self.index = np.searchsorted(self.values, last_value)

    def values_up_to(self, value):
        i = np.searchsorted(self.values, value)
        intermediate = self.values[self.index:i]
        self.index = i
        self.last_value = self.values[i - 1]
        return intermediate

    def __call__(self, value):
        return self.values_up_to(value)


class RawDataArrays3D(namedtuple("RawDataArrays3D", ['mz', 'intensity', 'ion_mobility'])):
    """Represent the m/z, intensity, and ion mobility arrays associated with a raw
    ion mobility frame of mass spectra.

    Thin wrapper around a ``namedtuple``, so this object supports
    the same interfaces as a tuple.

    Attributes
    ----------
    mz : :class:`np.ndarray`
        The m/z axis of a mass spectrum
    intensity : :class:`np.ndarray`
        The intensity measured at the corresponding m/z of a mass spectrum
    ion_mobility : :class:`np.ndarray`
        The ion mobility of each data point
    distinct_ion_mobility : :class:`np.ndarray`
        The distinct ion mobility values in the multidimensional space, independent
        of whether there are data points at that point
    ion_mobility_array_type : :class:`~.Term`
        The type of ion mobility array
    data_arrays : dict
        Any other data arrays
    """
    def __new__(cls, mz, intensity, ion_mobility, distinct_ion_mobility, ion_mobility_array_type=None, data_arrays=None):
        inst = super(RawDataArrays3D, cls).__new__(
            cls, mz, intensity, ion_mobility)
        inst.distinct_ion_mobility = distinct_ion_mobility
        inst.ion_mobility_array_type = ion_mobility_array_type
        inst.data_arrays = dict()
        if data_arrays:
            inst.data_arrays.update(data_arrays)
        return inst

    @property
    def drift_time(self):
        return self.ion_mobility

    def has_array(self, array_type):
        '''Check if this array set contains an array of the
        requested type.

        This method uses the semantic lookup mechanism to test
        "is-a" relationships so if a more abstract term is used,
        a wider range of terms may be matched.

        Parameters
        ----------
        array_type : str or :class:`~.Term`
            The array type name to test.

        Returns
        -------
        bool
        '''
        from ms_deisotope.data_source.metadata.scan_traits import binary_data_arrays
        try:
            term = binary_data_arrays[array_type]
        except KeyError:
            warnings.warn(
                "Array type %r could not be resolved, treating as a plain string" % (array_type, ))
            return array_type in self.binary_data_arrays
        if self.mz is not None and len(self.mz):
            k = binary_data_arrays['m/z array']
            if term.is_a(k):
                return k
        if self.intensity is not None and len(self.intensity):
            k = binary_data_arrays['intensity array']
            if term.is_a(k):
                return k
        if self.ion_mobility_array_type is not None:
            k = binary_data_arrays[self.ion_mobility_array_type]
            if term.is_a(k):
                return k
        for k in self.data_arrays:
            try:
                k = binary_data_arrays[k]
                if term.is_a(k):
                    return k
            except KeyError:
                if term == k:
                    return k
        return False

    def __copy__(self):
        inst = self.__class__(self.mz.copy(), self.intensity.copy(), self.ion_mobility.copy(),
                              self.distinct_ion_mobility.copy(), self.ion_mobility_array_type, {
            k: v.copy() for k, v in self.data_arrays.items()
        })
        return inst

    def copy(self):
        """Make a deep copy of this object.

        Returns
        -------
        :class:`RawDataArrays3D`
        """
        return self.__copy__()

    def __eq__(self, other):
        try:
            return np.allclose(
                self[0], other[0]) and np.allclose(
                    self[1], other[1]) and np.allclose(self[2], other[2])
        except ValueError:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __mul__(self, i):
        return self.__class__(self.mz, self.intensity * i, self.ion_mobility,
                              self.distinct_ion_mobility, self.ion_mobility_array_type, self.data_arrays)

    def __div__(self, d):
        return self.__class__(self.mz, self.intensity / d, self.ion_mobility,
                              self.distinct_ion_mobility, self.ion_mobility_array_type, self.data_arrays)
    def __getitem__(self, i):
        if isinstance(i, int):
            return super(RawDataArrays3D, self).__getitem__(i)
        elif isinstance(i, (slice, list, tuple, np.ndarray)):
            return self._slice(i)
        else:
            return self.data_arrays[i]

    @classmethod
    def stack(cls, scans, ion_mobility_array_type=None):
        '''Combine multiple :class:`~.Scan` objects or (ion mobility, :class:`~.RawDataArrays`)
        pairs into a single :class:`~.RawDataArrays3D`

        Parameters
        ----------
        scans : list of :class:`~.Scan`
            The single ion mobility point arrays to combine
        ion_mobility_array_type : :class:`~.Term` or :class:`str`
            The type of ion mobility array
        '''
        mz = array.array('d')
        intensity = array.array('d')
        ion_mobility = array.array('d')
        data_arrays = defaultdict(lambda: array.array('d'))

        distinct_ion_mobility = array.array('d')
        for scan in scans:
            try:
                arrays = scan.arrays
                im = scan.drift_time
            except AttributeError:
                im = scan[0]
                arrays = scan[1]
            distinct_ion_mobility.append(im)
            if arrays.mz.size > 0:
                mz.extend(arrays.mz)
                intensity.extend(arrays.intensity)
                ion_mobility.extend([im] * arrays.mz.size)
                for key, val in arrays.data_arrays.items():
                    data_arrays[key].extend(val)
        mz = np.array(mz)
        ion_mobility = np.array(ion_mobility)
        mask = np.lexsort(np.stack((ion_mobility, mz,)))
        return cls(mz[mask], np.array(intensity)[mask], ion_mobility[mask],
                   np.sort(distinct_ion_mobility), ion_mobility_array_type,
                   {k: np.array(v)[mask] for k, v in data_arrays.items()})

    @classmethod
    def from_arrays(cls, mz_array, intensity_array, ion_mobility_array, ion_mobility_array_type=None, data_arrays=None):
        '''Build a new :class:`~.RawDataArrays3D` from parallel arrays.

        This will sort all arrays w.r.t. m/z and ion mobility.

        Parameters
        ----------
        mz_array : np.ndarray
            The m/z array
        intensity_array : np.ndarray
            The intensity array
        ion_mobility_array : np.ndarray
            The ion mobility array
        data_arrays : dict
            Any extra data arrays

        Returns
        -------
        :class:`~.RawDataArrays3D`
        '''
        mz_array = np.asarray(mz_array)
        ion_mobility_array = np.asarray(ion_mobility_array)
        distinct_ion_mobility = np.sort(np.unique(ion_mobility_array))
        mask = np.lexsort(np.stack((ion_mobility_array, mz_array)))
        if data_arrays:
            data_arrays = {
                k: np.asarray(v)[mask] for k, v in data_arrays.items()
            }
        return cls(
            mz_array[mask], np.asarray(intensity_array)[mask], ion_mobility_array[mask],
            distinct_ion_mobility, ion_mobility_array_type, data_arrays=data_arrays)


    @classmethod
    def empty(cls):
        return cls(np.array([]), np.array([]), np.array([]), np.array([]))

    def _slice(self, i, include_ion_mobility=True):
        mz = self.mz[i]
        intensity = self.intensity[i]
        data_arrays = {k: v[i] for k, v in self.data_arrays.items()}
        if not include_ion_mobility:
            return RawDataArrays(mz, intensity, data_arrays)
        return RawDataArrays3D(
            mz, intensity, self.ion_mobility[i], self.distinct_ion_mobility,
            self.ion_mobility_array_type, data_arrays)

    def unstack(self, include_empty=True):
        '''Convert this 3D array into a list of (ion mobility, :class:`~.RawDataArrays`) pairs

        Parameters
        ----------
        include_empty : bool
            Whether to include ion mobility values with no signal. Defaults to :const:`True`.

        Returns
        -------
        arrays : list[tuple[float, :class:`~.RawDataArrays`]]
            The split arrays arranged by drift time
        '''
        # Sort by m/z first, then by ion mobility
        mask = np.lexsort(np.stack((self.mz, self.ion_mobility)))
        sorted_dim = self.ion_mobility[mask]
        acc = []
        # re-sort all arrays along the mask, and then slice out in chunks.
        boundaries = np.where(np.diff(sorted_dim) != 0)[0]
        if len(boundaries) == 0:
            return [RawDataArrays(self.mz, self.intensity, self.data_arrays)]
        if boundaries[-1] < len(sorted_dim) - 1:
            boundaries = np.append(boundaries, [len(sorted_dim) - 1])
        last_boundary = 0
        index = _upto_index(self.distinct_ion_mobility, -1)
        for boundary in boundaries:
            drift_time = self.ion_mobility[mask[last_boundary]]
            for i in index(drift_time):
                acc.append((i, self._slice(slice(0, 0), False)))
            chunk = self._slice(mask[last_boundary:boundary + 1], False)
            acc.append((drift_time, chunk))
            last_boundary = boundary + 1
        return acc

    def grid(self):
        '''Convert this 3D array into an intensity grid with
        m/z along the rows and ion mobility along the columns.

        Returns
        -------
        mz_axis : np.ndarray
            The m/z value for each row
        im_axis : np.ndarray
            The ion mobility value for each column
        intensity : np.ndarray
            A 2D array of intensity values at each axis position
        '''
        mz_axis = np.unique(self.mz)
        im_axis = (self.distinct_ion_mobility)
        intensity = np.zeros((mz_axis.shape[0], im_axis.shape[0]))
        i1 = np.searchsorted(mz_axis, self.mz)
        i2 = np.searchsorted(im_axis, self.ion_mobility)
        intensity[i1, i2] += self.intensity
        return mz_axis, im_axis, intensity

    def find_mz(self, mz):
        """Find the nearest index to the query ``mz``

        Parameters
        ----------
        mz : float
            The m/z value to search for

        Returns
        -------
        int
            The index nearest to the query m/z
        """
        n = len(self.mz)
        lo = 0
        hi = n

        while hi != lo:
            mid = int((hi + lo) // 2)
            y = self.mz[mid]
            err = y - mz
            if abs(err) < 0.1:
                best_index = mid
                best_err = abs(err)
                i = mid
                while i >= 0:
                    y = self.mz[i]
                    err = y - mz
                    if err <= -0.1:
                        break
                    abs_err = abs(err)
                    if abs_err < best_err:
                        best_err = abs_err
                        best_index = i
                    i -= 1
                i = mid
                while i < n:
                    y = self.mz[i]
                    err = y - mz
                    if err >= 0.1:
                        break
                    abs_err = abs(err)
                    if abs_err < best_err:
                        best_err = abs_err
                        best_index = i
                    i += 1
                return best_index
            elif hi - lo == 1:
                return mid
            elif err > 0:
                hi = mid
            else:
                lo = mid
        return 0

    def between_mz(self, low, high):
        """Returns a slice of the arrays between ``low`` and ``high``
        m/z

        Parameters
        ----------
        low : float
            The lower bound m/z
        high : float
            The upper bound m/z

        Returns
        -------
        :class:`.RawDataArrays`
        """
        i = self.find_mz(low)
        mz_i = self.mz[i]
        while self.mz[i] == mz_i:
            i -= 1
        i += 1
        j = self.find_mz(high)
        mz_j = self.mz[j]
        while self.mz[j] == mz_j:
            j += 1
        if not (low <= self.mz[i] <= high):
            i += 1
        return self._slice(slice(i, j), True)


class FrameBase(object):

    @property
    def scan_time(self):
        return self.time

    @property
    def ion_mobilities(self):
        return self.source._frame_ion_mobilities(self._data)

    def scans(self):
        return self.source._frame_scans(self._data)

    def bind(self, source):
        self.source = source
        if self.precursor_information is not None:
            self.precursor_information.bind(source)
        return self

    def unbind(self):
        self.source = None
        if self.precursor_information:
            self.precursor_information.unbind()
        return self

    def __repr__(self):
        template = ("{self.__class__.__name__}({self.id}, index={self.index},"
                    " time={self.time}, ms_level={self.ms_level})")
        return template.format(self=self)

    def __eq__(self, other):
        if other is None:
            return False
        if self.id != other.id:
            return False
        elif self.ms_level != other.ms_level:
            return False
        elif self.time != other.time:
            return False
        elif self.deconvoluted_features != other.deconvoluted_features:
            return False
        elif self.features != other.features:
            return False
        return True

    def __ne__(self, other):
        return not self == other

    __hash__ = None


class IonMobilityFrame(FrameBase):
    '''A :class:`IonMobilityFrame` represents an single time point acquisition of
    multiple mass spectra across multiple ion mobility drift time points. Because
    of its drift time by m/z structure, it does not have 1-D peak sets, but pseudo-
    :class:`~.LCMSFeatureMap`-like data structures which conserve the over-time
    property.
    '''
    def __init__(self, data, source):
        self._data = data
        self.source = source
        self.features = None
        self.deconvoluted_features = None

        self._unload()

    def _unload(self):
        self._id = None
        self._ms_level = None
        self._time = None
        self._polarity = None
        self._index = None
        self._start_scan_index = None
        self._end_scan_index = None

        self._arrays = None
        self._annotations = None

        self._isolation_window = None
        self._precursor_information = None
        self._activation = None
        self._acquisition_information = None

    @property
    def id(self):
        if self._id is None:
            self._id = self.source._frame_id(self._data)
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    @property
    def index(self):
        if self._index is None:
            self._index = self.source._frame_index(self._data)
        return self._index

    @index.setter
    def index(self, value):
        self._index = value

    @property
    def time(self):
        if self._time is None:
            self._time = self.source._frame_time(self._data)
        return self._time

    @time.setter
    def time(self, value):
        self._time = value

    @property
    def ms_level(self):
        if self._ms_level is None:
            self._ms_level = self.source._frame_ms_level(self._data)
        return self._ms_level

    @ms_level.setter
    def ms_level(self, value):
        self._ms_level = value

    @property
    def polarity(self):
        if self._polarity is None:
            self._polarity = self.source._frame_polarity(self._data)
        return self._polarity

    @polarity.setter
    def polarity(self, value):
        self._polarity = value

    @property
    def start_scan_index(self):
        if self._start_scan_index is None:
            self._start_scan_index = self.source._frame_start_scan_index(self._data)
        return self._start_scan_index

    @start_scan_index.setter
    def start_scan_index(self, value):
        self._start_scan_index = value

    @property
    def end_scan_index(self):
        if self._end_scan_index is None:
            self._end_scan_index = self.source._frame_end_scan_index(self._data)
        return self._end_scan_index

    @end_scan_index.setter
    def end_scan_index(self, value):
        self._end_scan_index = value

    @property
    def precursor_information(self):
        if self._precursor_information is None:
            self._precursor_information = self.source._frame_precursor_information(self._data)
        return self._precursor_information

    @precursor_information.setter
    def precursor_information(self, value):
        self._precursor_information = value

    @property
    def activation(self):
        if self._activation is None:
            self._activation = self.source._frame_activation(self._data)
        return self._activation

    @activation.setter
    def activation(self, value):
        self._activation = value

    @property
    def isolation_window(self):
        if self._isolation_window is None:
            self._isolation_window = self.source._frame_isolation_window(self._data)
        return self._isolation_window

    @isolation_window.setter
    def isolation_window(self, value):
        self._isolation_window = value

    @property
    def acquisition_information(self):
        if self._acquisition_information is None:
            self._acquisition_information = self.source._frame_acquisition_information(self._data)
        return self._acquisition_information

    @acquisition_information.setter
    def acquisition_information(self, value):
        self._acquisition_information = value

    @property
    def arrays(self):
        if self._arrays is None:
            self._arrays = self.source._frame_arrays(self._data)
        return self._arrays

    @arrays.setter
    def arrays(self, value):
        self._arrays = value

    @property
    def annotations(self):
        if self._annotations is None:
            self._annotations = self.source._frame_annotations(self._data)
        return self._annotations

    @annotations.setter
    def annotations(self, value):
        self._annotations = dict(value)

    def flatten_to_scan(self):
        scans = self.scans()
        mz, intensity = average_signal([scan.arrays for scan in scans], dx=0.001)
        intensity *= len(scans)
        arrays = RawDataArrays(mz, intensity)
        return WrappedScan(scans[0]._data, self.source, arrays, id="merged=%d" % self.index)

    def get_scan_by_drift_time(self, drift_time):
        dt_axis = self.ion_mobilities
        lo = 0
        hi = n = len(dt_axis)
        best_match = None
        best_error = float('inf')

        # Handles infinities
        if drift_time > dt_axis[-1]:
            return self.source.get_scan_by_index(self.end_scan_index - 1)

        while hi != lo:
            mid = (hi + lo) // 2
            dt = dt_axis[mid]
            err = abs(dt - drift_time)
            if err < best_error:
                best_error = err
                best_match = mid
            if (dt == drift_time) or (hi - 1 == lo):
                i = mid - 1
                while i >= 0:
                    err = abs(dt_axis[i] - drift_time)
                    if err < best_error:
                        best_error = err
                        best_match = i
                        i -= 1
                    else:
                        break
                i = mid + 1
                while i < n:
                    err = abs(dt_axis[i] - drift_time)
                    if err < best_error:
                        best_error = err
                        best_match = i
                        i += 1
                    else:
                        break
                return self.source.get_scan_by_index(self.start_scan_index + best_match)
            elif dt > drift_time:
                hi = mid
            else:
                lo = mid

    def extract_features(self, error_tolerance=1.5e-5, max_gap_size=0.25, min_size=2, average=0, dx=0.001, **kwargs):
        from ms_deisotope.feature_map import feature_map
        scans = self.scans()
        lff = feature_map.LCMSFeatureForest(error_tolerance=error_tolerance)
        n = len(scans)
        for i, scan in enumerate(scans):
            if average:
                acc = []
                if i > 0:
                    for j in range(max((i - average, 0)), i):
                        acc.append(scans[j])
                    for j in range(i + 1, min(i + average + 1, n)):
                        acc.append(scans[j])
                scan = scan.average_with(acc, dx=dx)
            scan.pick_peaks(**kwargs)
            for peak in scan:
                lff.handle_peak(peak, scan.drift_time)
        lff.split_sparse(max_gap_size, min_size).smooth_overlaps(
            error_tolerance)
        self.features = lff
        return self

    def deconvolute_features(self, averagine=None, scorer=None, truncate_after=0.95,
                             minimum_intensity=5, min_size=2, max_gap_size=0.25, **kwargs):
        from ms_deisotope.feature_map import feature_processor
        if averagine is None:
            from ms_deisotope.averagine import peptide
            averagine = peptide
        if scorer is None:
            from ms_deisotope.scoring import MSDeconVFitter
            scorer = MSDeconVFitter(1)
        if self.features is None:
            raise ValueError(
                "IM-MS Features must be extracted before they can be charge state deconvoluted")
        decon = feature_processor.LCMSFeatureProcessor(
            self.features, averagine, scorer, minimum_size=min_size,
            maximum_time_gap=max_gap_size)
        self.deconvoluted_features = decon.deconvolute(
            minimum_intensity=minimum_intensity, truncate_after=truncate_after, **kwargs)
        return self

    def to_deconvoluted_peak_set(self, time_bound=None):
        if self.deconvoluted_features is None:
            raise ValueError("Must first deconvolute IMS features before converting to a time-binned peak set!")
        if time_bound is not None:
            cft = IntervalTreeNode.build(
                [Interval(f.start_time, f.end_time, [f]) for f in self.deconvoluted_features])
            features = chain.from_iterable(cft.overlaps(*time_bound))
        else:
            features = self.deconvoluted_features
        return features_to_peak_set(features)

    def pack(self, bind=False):
        packed = ProcessedIonMobilityFrame(
            self.id, self.precursor_information, self.ms_level,
            self.time, self.index, self.features, self.deconvoluted_features,
            self.start_scan_index, self.end_scan_index,
            self.polarity, self.activation, self.acquisition_information,
            self.isolation_window, self.annotations,
            self.source if bind else None)
        return packed


class ProcessedIonMobilityFrame(FrameBase):
    def __init__(self, id, precursor_information, ms_level, time, index,
                 features=None, deconvoluted_features=None,
                 start_scan_index=None, end_scan_index=None,
                 polarity=None, activation=None, acquisition_information=None,
                 isolation_window=None, annotations=None, source=None):
        if annotations is None:
            annotations = dict()
        self.source = source

        self.id = id
        self.precursor_information = precursor_information
        self.ms_level = ms_level
        self.time = time
        self.index = index
        self.features = features
        self.deconvoluted_features = deconvoluted_features

        self.start_scan_index = start_scan_index
        self.end_scan_index = end_scan_index
        self.polarity = polarity
        self.activation = activation
        self.acquisition_information = acquisition_information
        self.isolation_window = isolation_window
        self.annotations = annotations


def weighted_centroid(feature):
    total = 0
    normalizer = 0
    for node in feature:
        weight = node.total_intensity()
        total += node.time * weight
        normalizer += weight
    return total / normalizer


def merge_envelopes(envelopes):
    base = envelopes[0].clone()
    for env in envelopes[1:]:
        for i, p in enumerate(env):
            base[i].intensity += p.intensity

    return base


def feature_to_peak(feature):
    peak_cluster = feature.peaks
    peak_cluster = [pi for p in peak_cluster for pi in p]
    total_intensity = sum(p.intensity for p in peak_cluster)
    mz = sum(p.mz * p.intensity for p in peak_cluster) / total_intensity
    neutral_mass = sum(
        p.neutral_mass * p.intensity for p in peak_cluster) / total_intensity
    most_abundant_mass = sum(
        p.most_abundant_mass * p.intensity for p in peak_cluster) / total_intensity
    a_to_a2_ratio = sum(
        p.a_to_a2_ratio * p.intensity for p in peak_cluster) / total_intensity
    average_mass = sum(
        p.average_mass * p.intensity for p in peak_cluster) / total_intensity
    signal_to_noise = sum(p.signal_to_noise *
                          p.intensity for p in peak_cluster) / total_intensity
    fwhm = sum(p.full_width_at_half_max *
               p.intensity for p in peak_cluster) / total_intensity
    area = sum(p.area * p.intensity for p in peak_cluster) / total_intensity
    score = sum(p.score * p.intensity for p in peak_cluster) / total_intensity
    charge = peak_cluster[0].charge
    envelope = merge_envelopes([p.envelope for p in peak_cluster])
    return IonMobilityDeconvolutedPeak(
        neutral_mass=neutral_mass, intensity=total_intensity, charge=charge, signal_to_noise=signal_to_noise,
        full_width_at_half_max=fwhm, index=-1, a_to_a2_ratio=a_to_a2_ratio, most_abundant_mass=most_abundant_mass,
        average_mass=average_mass, score=score, envelope=envelope, mz=mz, fit=None, chosen_for_msms=False,
        area=area, drift_time=weighted_centroid(feature))


def features_to_peak_set(features):
    peak_set = DeconvolutedPeakSet(tuple(map(feature_to_peak, features)))
    peak_set.reindex()
    return peak_set


class Generic3DIonMobilityFrameSource(IonMobilitySourceRandomAccessFrameSource):
    def __init__(self, loader, **kwargs):
        self.loader = loader
        self.initialize_frame_cache()
        self._producer = self._wrap_iterator(self.loader)

    def _frame_acquisition_information(self, data):
        return self.loader._acquisition_information(data)

    def _frame_activation(self, data):
        return self.loader._activation(data)

    def _frame_id(self, data):
        return self.loader._scan_id(data)

    def _frame_index(self, data):
        return self.loader._scan_index(data)

    def _frame_time(self, data):
        return self.loader._scan_time(data)

    def _frame_precursor_information(self, data):
        return self.loader._precursor_information(data)

    def _frame_ms_level(self, data):
        return self.loader._ms_level(data)

    def _frame_polarity(self, data):
        return self.loader._polarity(data)

    def _frame_isolation_window(self, data):
        return self.loader._isolation_window(data)

    def _frame_start_scan_index(self, data):
        return self._frame_index(data)

    def _frame_end_scan_index(self, data):
        return self._frame_index(data) + 1

    def _frame_arrays(self, data):
        from ms_deisotope.data_source.metadata.scan_traits import binary_data_arrays

        mz, intensity, data_arrays = self.loader._scan_arrays(data)
        ion_mobility = None
        ion_mobility_type = None
        for key, value in data_arrays.items():
            try:
                acc = key.accession
                term = binary_data_arrays[acc]
            except (KeyError, AttributeError):
                continue
            if term.is_a('ion mobility array'):
                ion_mobility = value
                ion_mobility_type = key
                break
        if ion_mobility_type is None:
            raise ValueError(
                "Expected an ion mobility array type, but none were found.")
        d = data_arrays.copy()
        d.pop(ion_mobility_type)
        array3d = RawDataArrays3D.from_arrays(
            mz, intensity, ion_mobility, ion_mobility_type, d)
        return array3d

    def _frame_scans(self, data):
        from ms_deisotope.data_source.metadata.scan_traits import ion_mobility_attribute

        arrays = self._frame_arrays(data)
        scan_id_base = self._frame_id(data)
        i = 1
        scans = []
        for drift_time, arrays_i in arrays.unstack():
            subdat = data.copy()
            subdat['m/z array'] = arrays_i.mz
            subdat['intensity array'] = arrays_i.intensity
            subdat['id'] = scan_id_base + (' scan=%d' % i)
            i += 1
            scan = self.loader._make_scan(subdat)
            acq = scan.acquisition_information
            event = acq[0]
            tp = event.ion_mobility_type
            if tp:
                event._ion_mobility.remove_ion_mobility_type(tp)
            event._ion_mobility.add_ion_mobility(
                ion_mobility_attribute, drift_time)
            scans.append(scan)
        return scans

    def get_frame_by_index(self, index):
        scan = self.loader.get_scan_by_index(index)
        frame = self._make_frame(scan._data)
        self._cache_frame(frame)
        return frame

    def get_frame_by_id(self, id):
        scan = self.loader.get_scan_by_id(id)
        frame = self._make_frame(scan._data)
        self._cache_frame(frame)
        return frame

    def get_frame_by_time(self, time):
        scan = self.loader.get_scan_by_time(time)
        frame = self._make_frame(scan._data)
        self._cache_frame(frame)
        return frame

    def __len__(self):
        return len(self.loader)

    def __getitem__(self, i):
        result = self.loader[i]
        if isinstance(result, list):
            out = []
            for scan in result:
                frame = self._make_frame(scan._data)
                self._cache_frame(frame)
                out.append(frame)
            return out
        else:
            scan = result
            frame = self._make_frame(scan._data)
            self._cache_frame(frame)
            return frame

    def _wrap_iterator(self, iterator):
        for val in iterator:
            if isinstance(val, ScanBunch):
                precursor = self._make_frame(val.precursor._data)
                self._cache_frame(precursor)
                products = []
                for product in val.products:
                    product = self._make_frame(product._data)
                    self._cache_frame(product)
                    products.append(product)

                yield ScanBunch(precursor, products)
            else:
                frame = self._make_frame(val._data)
                self._cache_frame(frame)
                yield frame

    def _default_frame_iterator(self, start_index=None):
        self.loader.start_from_scan(index=start_index, grouped=False)
        return self._wrap_iterator(self.loader)

    def make_frame_iterator(self, iterator=None, grouped=False):
        self.loader.make_iterator(iterator, grouped=grouped)
        return self._wrap_iterator(self.loader)

    def start_from_frame(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        self.loader.start_from_scan(
            scan_id=scan_id, rt=rt, index=index, require_ms1=require_ms1, grouped=grouped)
        self._producer = self._wrap_iterator(self.loader)
        return self

    def __next__(self):
        return next(self._producer)

    def next(self):
        return next(self._producer)

    def reset(self):
        self.loader.reset()
        self.initialize_frame_cache()
        return self


class FramedIonMobilityFrameSource(IonMobilitySourceRandomAccessFrameSource):
    def __init__(self, loader, **kwargs):
        self.loader = loader
        self.initialize_frame_cache()
        self.scan_index_range_to_frame = dict()

    def _find_cycle_start(self, scan):
        dt = scan.drift_time
        if dt is None:
            return scan
        index = scan.index
        best_scan = scan
        while index >= 0:
            prev = self.loader[index]
            p_dt = prev.drift_time
            if p_dt is None:
                return best_scan
            if p_dt < dt:
                best_scan = scan
            else:
                return best_scan
            index -= 1
        return best_scan

    def _collect_cycle(self, scan):
        acc = []
        start_scan = self._find_cycle_start(scan)
        dt = start_scan.drift_time
        if dt is None:
            raise ValueError("Could not find a start scan with drift time from %r" % (scan, ))
        acc.append(start_scan)
        index = start_scan.index + 1
        while True:
            next_scan = self.loader.get_scan_by_index(index)
            next_dt = next_scan.drift_time
            if next_dt is None:
                break
            if next_dt < dt:
                break
            dt = next_dt
            acc.append(next_scan)
            index += 1
        return acc

    def _frame_id(self, frame):
        scans = self._frame_scans(frame)
        return ','.join([s.id for s in scans])

    def _frame_acquisition_information(self, data):
        return self.loader._acquisition_information(data)

    def _frame_activation(self, data):
        scans = self._frame_scans(data)
        return scans[0].activation

    def _frame_index(self, data):
        scans = self._frame_scans(data)
        return scans[0].index

    def _frame_time(self, data):
        scans = self._frame_scans(data)
        return scans[0].scan_time

    def _frame_ms_level(self, data):
        scans = self._frame_scans(data)
        return scans[0].ms_level

    def _frame_polarity(self, data):
        scans = self._frame_scans(data)
        return scans[0].polarity

    def _frame_isolation_window(self, data):
        scans = self._frame_scans(data)
        return scans[0].isolation_window

    def _frame_start_scan_index(self, data):
        scans = self._frame_scans(data)
        return scans[0].index

    def _frame_end_scan_index(self, data):
        scans = self._frame_scans(data)
        return scans[-1].index + 1

    def _frame_scans(self, data):
        try:
            return data['scans']
        except KeyError:
            pass
        source_id = data['source_scan_id']
        scan = self.loader.get_scan_by_id(source_id)
        scans = self._collect_cycle(scan)
        data['scans'] = scans
        data['scan_index_range'] = (scans[0].index, scans[-1].index + 1)
        return data['scans']

    def _make_frame_from_scan(self, scan):
        return {
            "source_scan_id": scan.id,
            "source_scan_index": scan.index,
        }

    def get_scan_by_index(self, index):
        return self.loader.get_scan_by_index(index)
