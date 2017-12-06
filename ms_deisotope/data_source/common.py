import abc
import warnings

from collections import namedtuple

import numpy as np

from ms_peak_picker import (
    pick_peaks, reprofile, average_signal, scan_filter)
from ms_peak_picker import PeakIndex, PeakSet
from ..averagine import neutral_mass, mass_charge_ratio
from ..utils import Constant, add_metaclass
from ..deconvolution import deconvolute_peaks

from .metadata.instrument_components import Component, component, all_components
from .metadata.file_information import FileInformation, SourceFile

try:
    from ..utils import draw_raw, draw_peaklist, annotate_scan as _annotate_precursors
    has_plot = True
except Exception:
    has_plot = False


class ScanBunch(namedtuple("ScanBunch", ["precursor", "products"])):

    def __new__(cls, *args, **kwargs):
        inst = super(ScanBunch, cls).__new__(cls, *args, **kwargs)
        inst._id_map = {}
        if inst.precursor is not None:
            inst._id_map[inst.precursor.id] = inst.precursor
        for scan in inst.products:
            inst._id_map[scan.id] = scan
        return inst

    def precursor_for(self, scan):
        scan_id = scan.precursor_information.precursor_scan_id
        return self.get_scan_by_id(scan_id)

    def get_scan_by_id(self, scan_id):
        return self._id_map[scan_id]

    def annotate_precursors(self, nperrow=4, ax=None):
        return _annotate_precursors(
            self.precursor, self.products, nperrow=nperrow, ax=ax)

    def _repr_pretty_(self, p, cycle):  # pragma: no cover
        if cycle:
            p.text("ScanBunch(...)")
            return
        p.text("ScanBunch(\n")
        with p.group(2):
            with p.group(4, "precursor=\n"):
                p.pretty(self.precursor)
            with p.group(4, ",\nproducts=\n"):
                p.pretty(self.products)
        p.text(")")


class RawDataArrays(namedtuple("RawDataArrays", ['mz', 'intensity'])):

    def plot(self, *args, **kwargs):
        ax = draw_raw(self, *args, **kwargs)
        return ax


DEFAULT_CHARGE_WHEN_NOT_RESOLVED = 1
ChargeNotProvided = Constant("ChargeNotProvided")


@add_metaclass(abc.ABCMeta)
class ScanDataSource(object):
    """An Abstract Base Class describing an object
    which can provide a consistent set of accessors
    for a particular format of mass spectrometry data.

    Data files come in many shapes and sizes, with different
    underlying structures. This class provides an API that
    should make features as consistent as possible to clients
    of Scan objects, making the format those Scan objects
    were read from unimportant.
    """
    @abc.abstractmethod
    def _scan_arrays(self, scan):
        """Returns raw data arrays for m/z and intensity

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        mz: np.array
            An array of m/z values for this scan
        intensity: np.array
            An array of intensity values for this scan
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _precursor_information(self, scan):
        """Returns information about the precursor ion,
        if any, that this scan was derived form.

        Returns `None` if this scan has no precursor ion

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        PrecursorInformation
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _scan_title(self, scan):
        """Returns a verbose name for this scan, if one
        were stored in the file. Usually includes both the
        scan's id string, as well as information about the
        original file and format.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        str
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _scan_id(self, scan):
        """Returns the scan's id string, a unique
        identifier for this scan in the context of
        the data file it is recordered in

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        str
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _scan_index(self, scan):
        """Returns the base 0 offset from the start
        of the data file in number of scans to reach
        this scan.

        If the original format does not natively include
        an index value, this value may be computed from
        the byte offset index.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        int
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _ms_level(self, scan):
        """Returns the degree of exponential fragmentation
        used to produce this scan. 1 refers to a survey scan
        of unfragmented ions, 2 refers to a tandem scan derived
        from an ms level 1 ion, and so on.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        int
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _scan_time(self, scan):
        """Returns the time in minutes from the start of data
        acquisition to when this scan was acquired.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        float
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _is_profile(self, scan):
        """Returns whether the scan contains profile data (`True`)
        or centroided data (`False`).

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        bool
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _polarity(self, scan):
        """Returns whether this scan was acquired in positive mode (+1)
        or negative mode (-1).

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        int
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _activation(self, scan):
        """Returns information about the activation method used to
        produce this scan, if any.

        Returns `None` for MS1 scans

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        ActivationInformation
        """
        raise NotImplementedError()

    def _acquisition_information(self, scan):
        return None

    def _isolation_window(self, scan):
        return None

    def _instrument_configuration(self, scan):
        return None

    def _annotations(self, scan):
        return dict()


@add_metaclass(abc.ABCMeta)
class ScanIterator(ScanDataSource):
    """An Abstract Base Class that extends ScanDataSource
    with additional requirements that enable clients of the
    class to treat the object as an iterator over the underlying
    data file.
    """

    iteration_mode = 'single'

    @abc.abstractmethod
    def next(self):
        raise NotImplementedError()

    def _make_scan(self, data):
        return Scan(data, self)

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self

    def reset(self):
        raise NotImplementedError()

    def _make_default_iterator(self):
        raise NotImplementedError()

    def make_iterator(self, iterator=None, grouped=True):
        if grouped:
            self._producer = self._scan_group_iterator(iterator)
            self.iteration_mode = 'group'
        else:
            self._producer = self._single_scan_iterator(iterator)
            self.iteration_mode = 'single'

    def _single_scan_iterator(self, iterator=None):
        if iterator is None:
            iterator = self._make_default_iterator()

        _make_scan = self._make_scan

        for scan in iterator:
            packed = _make_scan(scan)
            if not self._validate(packed):
                continue
            self._scan_cache[packed.id] = packed
            yield packed

    def _scan_group_iterator(self, iterator=None):
        if iterator is None:
            iterator = self._make_default_iterator()
        precursor_scan = None
        product_scans = []

        current_level = 1

        _make_scan = self._make_scan

        for scan in iterator:
            packed = _make_scan(scan)
            if not self._validate(packed):
                continue
            self._scan_cache[packed.id] = packed
            if packed.ms_level > 1:
                # inceasing ms level
                if current_level < packed.ms_level:
                    current_level = packed.ms_level
                # decreasing ms level
                elif current_level > packed.ms_level:
                    current_level = packed.ms_level.ms_level
                product_scans.append(packed)
            elif packed.ms_level == 1:
                if current_level > 1:
                    precursor_scan.product_scans = list(product_scans)
                    yield ScanBunch(precursor_scan, product_scans)
                else:
                    if precursor_scan is not None:
                        precursor_scan.product_scans = list(product_scans)
                        yield ScanBunch(precursor_scan, product_scans)
                precursor_scan = packed
                product_scans = []
            else:
                raise ValueError("Could not interpret MS Level %r" % (packed.ms_level,))
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)


@add_metaclass(abc.ABCMeta)
class RandomAccessScanSource(ScanDataSource):

    @abc.abstractmethod
    def get_scan_by_id(self, scan_id):
        raise NotImplementedError()

    @abc.abstractmethod
    def get_scan_by_time(self, time):
        raise NotImplementedError()

    @abc.abstractmethod
    def get_scan_by_index(self, index):
        raise NotImplementedError()

    @abc.abstractmethod
    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        raise NotImplementedError()

    def _scan_cleared(self, scan):
        pass

    def _locate_ms1_scan(self, scan):
        while scan.ms_level != 1:
            if scan.index <= 0:
                break
            scan = self.get_scan_by_index(scan.index - 1)
        if scan.ms_level == 1:
            return scan
        while scan.ms_level != 1:
            try:
                scan = self.get_scan_by_index(scan.index + 1)
            except IndexError:
                raise IndexError("Cannot locate MS1 Scan")
        return scan

    def find_previous_ms1(self, start_index):
        index = start_index - 1
        while index >= 0:
            try:
                scan = self.get_scan_by_index(index)
                if scan.ms_level == 1:
                    return scan
                index -= 1
            except (IndexError, KeyError):
                return None
        return None

    def find_next_ms1(self, start_index):
        index = start_index + 1
        n = len(self.index)
        while index < n:
            try:
                scan = self.get_scan_by_index(index)
                if scan.ms_level == 1:
                    return scan
                index += 1
            except (IndexError, KeyError):
                return None
        return None


class DetachedAccessError(Exception):
    pass


class DataAccessProxy(object):
    def __init__(self, source):
        self.source = None
        self.attach(source)

    def attach(self, source):
        self.source = source

    def detach(self):
        self.source = None

    def __getstate__(self):
        return ()

    def __setstate__(self, state):
        self.source = None

    def raise_if_detached(self):
        if self.source is None:
            raise DetachedAccessError("Cannot perform operation. Instance is detached.")

    def get_scan_by_id(self, scan_id):
        self.raise_if_detached()
        return self.source.get_scan_by_id(scan_id)


class ScanBase(object):

    def has_ion_mobility(self):
        acq = self.acquisition_information
        if acq is None:
            return False
        scan_event = acq[0]
        return scan_event.has_ion_mobility()

    @property
    def drift_time(self):
        acq = self.acquisition_information
        if acq is None:
            return None
        scan_event = acq[0]
        return scan_event.drift_time

    def copy(self):
        return self.clone()

    def __copy__(self):
        return self.clone()


class Scan(ScanBase):
    """Container for mass spectral data and associated descriptive information.

    A :class:`Scan` object is a generic object intended to be created by any `ScanDataSource` and describes
    a mass spectrum at each level of processing (Profile -> Peak Fitted -> Deconvoluted). The raw object
    provided by the source is wrapped and queried lazily when an attribute is requested, and retrieved
    using methods specific to that source type.

    Attributes
    ----------
    deconvoluted_peak_set : ms_deisotope.peak_set.DeconvolutedPeakSet or None
        Deconvoluted peaks resulting from charge state deconvolution and deisotoping. Will
        be `None` if deconvolution has not been done.
    peak_set : ms_peak_picker.peak_index.PeakIndex or None
        Picked peaks and (possibly) associated raw data points as produced by :meth:`pick_peaks`.
        Will be `None` if peak picking has not been done.
    product_scans : list
        A list of :class:`Scan` instances which were produced by fragmenting ions from this one.
    source : ScanDataSource
        The object which produced this scan and which defines the methods for retrieving common
        attributes from the underlying data structures.
    precursor_information: PrecursorInformation or None
        Descriptive metadata for the ion which was chosen for fragmentation, and a reference to
        the precursor scan
    arrays: list of numpy.ndarray
        A pair of `numpy.ndarray` objects corresponding to the raw m/z and intensity data points
    id: str
        The unique identifier for this scan as given by the source
    title: str
        The human-readable display string for this scan as shown in some external software
    ms_level: int
        The degree of fragmentation performed. 1 corresponds to a MS1 or "Survey" scan, 2 corresponds
        to MS/MS, and so on. If `ms_level` > 1, the scan is considered a "tandem scan" or "MS^n" scan
    scan_time: float
        The time the scan was acquired during data acquisition. The unit of time depends upon the source.
        May be minutes or seconds.
    index: int
        The integer number indicating how many scans were acquired prior to this scan.
    is_profile: bool
        Whether this scan's raw data points corresponds to a profile scan or whether the raw data was
        pre-centroided.
    polarity: int
        If the scan was acquired in positive mode, the value `+1`.  If the scan was acquired in negative
        mode, the value `-1`. May be used to indicating how to calibrate charge state determination methods.
    activation: ActivationInformation or None
        If this scan is an MS^n scan, this attribute will contain information about the process
        used to produce it from its parent ion.
    acquisition_information: ScanAcquisitionInformation or None
        Describes the type of event that produced this scan, as well as the scanning method
        used.
    isolation_window: IsolationWindow or None:
        Describes the range of m/z that were isolated from a parent scan to create this scan
    """
    def __init__(self, data, source, peak_set=None, deconvoluted_peak_set=None, product_scans=None, annotations=None):
        if product_scans is None:
            product_scans = []
        if annotations is None:
            annotations = dict()
        self.source = source
        self.peak_set = peak_set
        self.deconvoluted_peak_set = deconvoluted_peak_set

        self._data = data

        self._arrays = None
        self._id = None
        self._title = None
        self._ms_level = None
        self._scan_time = None
        self._precursor_information = None
        self._index = None
        self._is_profile = None
        self._polarity = None
        self._activation = None
        self._acquisition_information = None
        self._isolation_window = None
        self._instrument_configuration = None

        self._annotations = None
        self._external_annotations = annotations

        self.product_scans = product_scans

    def clone(self):
        dup = self.__class__(
            self._data, self.source,
            self.peak_set.clone() if self.peak_set is not None else None,
            self.deconvoluted_peak_set.clone() if self.deconvoluted_peak_set is not None else None,
            [s.clone() for s in self.product_scans], self._external_annotations.copy())
        return dup

    def _load(self):
        self.arrays
        self.id
        self.title
        self.ms_level
        self.scan_time
        self.index
        self.polarity
        self.precursor_information
        self.activation
        self.acquisition_information

    def _unload(self):
        self._arrays = None
        self._id = None
        self._title = None
        self._ms_level = None
        self._scan_time = None
        self._precursor_information = None
        self._index = None
        self._is_profile = None
        self._polarity = None
        self._activation = None
        self._acquisition_information = None
        self._isolation_window = None
        self._instrument_configuration = None

    def clear(self):
        if self.source is not None:
            self.source._scan_cleared(self)
        self._unload()

    @property
    def ms_level(self):
        if self._ms_level is None:
            self._ms_level = self.source._ms_level(self._data)
        return self._ms_level

    @property
    def is_profile(self):
        if self._is_profile is None:
            self._is_profile = self.source._is_profile(self._data)
        return self._is_profile

    @property
    def polarity(self):
        if self._polarity is None:
            self._polarity = self.source._polarity(self._data)
        return self._polarity

    @property
    def scan_time(self):
        if self._scan_time is None:
            self._scan_time = self.source._scan_time(self._data)
        return self._scan_time

    @property
    def arrays(self):
        if self._arrays is None:
            self._arrays = RawDataArrays(*self.source._scan_arrays(self._data))
        return self._arrays

    @property
    def title(self):
        if self._title is None:
            self._title = self.source._scan_title(self._data)
        return self._title

    @property
    def id(self):
        if self._id is None:
            self._id = self.source._scan_id(self._data)
        return self._id

    @property
    def scan_id(self):
        return self.id

    @property
    def index(self):
        if self._index is None:
            self._index = self.source._scan_index(self._data)
        return self._index

    @index.setter
    def index(self, value):
        self._index = value

    @property
    def precursor_information(self):
        if self.ms_level < 2:
            return None
        if self._precursor_information is None:
            self._precursor_information = self.source._precursor_information(self._data)
        return self._precursor_information

    @property
    def activation(self):
        if self.ms_level < 2:
            return None
        if self._activation is None:
            self._activation = self.source._activation(self._data)
        return self._activation

    @property
    def isolation_window(self):
        if self.ms_level < 2:
            return None
        if self._isolation_window is None:
            self._isolation_window = self.source._isolation_window(self._data)
        return self._isolation_window

    @property
    def acquisition_information(self):
        if self._acquisition_information is None:
            self._acquisition_information = self.source._acquisition_information(self._data)
        return self._acquisition_information

    @property
    def instrument_configuration(self):
        if self._instrument_configuration is None:
            self._instrument_configuration = self.source._instrument_configuration(
                self._data)
        return self._instrument_configuration

    @property
    def annotations(self):
        if self._annotations is None:
            self._annotations = self.source._annotations(self._data)
            self._annotations.update(self._external_annotations)
        return self._annotations

    def __repr__(self):
        return "Scan(%r, index=%d, time=%0.4f, ms_level=%r%s)" % (
            self.id, self.index, self.scan_time, self.ms_level,
            ", " + repr(self.precursor_information) if self.precursor_information else '')

    # peak manipulation

    def __iter__(self):
        if self.peak_set is None:
            raise ValueError("Cannot iterate over peaks in a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        return iter(self.peak_set)

    def has_peak(self, *args, **kwargs):
        """A wrapper around :meth:`ms_peak_picker.PeakIndex.has_peak` to query the
        :class:`ms_peak_picker.FittedPeak` objects picked for this scan.

        Parameters
        ----------
        mz: float
            The m/z to search for
        error_tolerance: float
            The parts per million mass error tolerance to use

        Returns
        -------
        ms_peak_picker.FittedPeak or None
            The peak closest to the query m/z within the error tolerance window, or None if not found
            or if peaks have not yet been picked

        Raises
        ------
        ValueError:
            If the scan has not yet had peaks picked yet

        See Also
        --------
        :meth:`.Scan.pick_peaks`
        """
        if self.peak_set is None:
            raise ValueError("Cannot search for peaks in a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        return self.peak_set.has_peak(*args, **kwargs)

    def pick_peaks(self, *args, **kwargs):
        """A wrapper around :func:`ms_peak_picker.pick_peaks` which will populate the
        :attr:`peak_set` attribute of this scan.

        Parameters
        ----------
        *args :
            Passed along to :func:`ms_peak_picker.pick_peaks`
        **kwargs :
            Passed along to :func:`ms_peak_picker.pick_peaks`

        Returns
        -------
        self
        """
        mzs, intensities = self.arrays
        if len(mzs) == 0:
            self.peak_set = PeakIndex(mzs, intensities, PeakSet([]))
            return self
        if self.is_profile:
            peak_mode = 'profile'
        else:
            peak_mode = 'centroid'

        kwargs.setdefault('peak_mode', peak_mode)

        self.peak_set = pick_peaks(mzs, intensities, *args, **kwargs)
        return self

    def deconvolute(self, *args, **kwargs):
        if self.peak_set is None:
            raise ValueError("Cannot deconvolute a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        decon_results = deconvolute_peaks(self.peak_set, *args, **kwargs)
        self.deconvoluted_peak_set = decon_results.peak_set
        return self

    def pack(self):
        precursor_info = self.precursor_information
        return ProcessedScan(
            self.id, self.title, precursor_info,
            self.ms_level, self.scan_time, self.index,
            self.peak_set.pack() if self.peak_set is not None else None,
            self.deconvoluted_peak_set,
            self.polarity,
            self.activation,
            self.acquisition_information,
            self.isolation_window,
            self.instrument_configuration,
            self.product_scans,
            self.annotations)

    # signal transformation

    def reprofile(self, max_fwhm=0.2, dx=0.01, model_cls=None):
        if self.peak_set is None and self.is_profile:
            raise ValueError("Cannot reprofile a scan that has not been centroided")
        elif self.peak_set is None and not self.is_profile:
            self.pick_peaks()
        arrays = reprofile(self.peak_set, max_fwhm, dx, model_cls)
        scan = WrappedScan(
            self._data, self.source, arrays,
            list(self.product_scans), is_profile=True,
            annotations=self._external_annotations)
        return scan

    def denoise(self, scale=5.0, window_length=2.0, region_width=10):
        mzs, intensities = self.arrays
        mzs = mzs.astype(float)
        intensities = intensities.astype(float)
        transform = scan_filter.FTICRBaselineRemoval(
            window_length=window_length, scale=scale, region_width=region_width)
        mzs, intensities = transform(mzs, intensities)
        return WrappedScan(self._data, self.source,
                           (mzs, intensities), list(self.product_scans),
                           annotations=self._external_annotations)

    def transform(self, filters=None):
        mzs, intensities = self.arrays
        mzs = mzs.astype(float)
        intensities = intensities.astype(float)
        mzs, intensities = scan_filter.transform(mzs, intensities, filters=filters)
        return WrappedScan(self._data, self.source,
                           (mzs, intensities), list(self.product_scans),
                           is_profile=True,
                           annotations=self._external_annotations)

    def _average_with(self, scans, dx=0.01, weight_sigma=None):
        scans = [self] + list(scans)
        arrays = []
        for scan in scans:
            if scan.is_profile:
                arrays.append(scan.arrays)
            else:
                arrays.append(scan.reprofile(dx=dx).arrays)
        if weight_sigma:
            if weight_sigma == 1:
                weight_sigma = 0.025
            weights = self._compute_smoothing_weights(
                scans, mean=self.scan_time, sigma=weight_sigma)
        else:
            weights = None
        new_arrays = average_signal(arrays, dx=dx, weights=weights)
        indices = [scan.index for scan in scans]
        return AveragedScan(
            self._data, self.source, new_arrays,
            indices, list(self.product_scans), is_profile=True,
            annotations=self._external_annotations)

    def _get_adjacent_scans(self, index_interval=None, rt_interval=None):
        if index_interval is None and rt_interval is None:
            raise ValueError("One of `index_interval` or `rt_interval` must be provided")
        if self.ms_level > 1:
            raise ValueError("Cannot average MSn scans at this time")
        if not self.source:
            raise ValueError("Can't average an unbound scan")
        before = []
        after = []
        if index_interval is not None:
            before = []
            current_index = self.index
            for i in range(index_interval):
                next_scan = self.source.find_previous_ms1(current_index)
                if next_scan is None:
                    break
                before.append(next_scan)
                current_index = next_scan.index
            before = before[::-1]
            after = []
            current_index = self.index
            for i in range(index_interval):
                try:
                    next_scan = self.source.find_next_ms1(current_index)
                except ValueError:
                    break
                if next_scan is None:
                    break
                after.append(next_scan)
                current_index = next_scan.index
        elif rt_interval is not None:
            reference_time = self.scan_time
            before = []
            current_index = self.index
            current_time = self.scan_time
            while abs(reference_time - current_time) < rt_interval and current_index > 0:
                next_scan = self.source.find_previous_ms1(current_index)
                if next_scan is None:
                    break
                before.append(next_scan)
                current_index = next_scan.index
                current_time = next_scan.scan_time

            before = before[::-1]

            after = []
            current_index = self.index
            current_time = self.scan_time
            while abs(reference_time - current_time) < rt_interval and current_index > 0:
                try:
                    next_scan = self.source.find_next_ms1(current_index)
                except ValueError:
                    break
                if next_scan is None:
                    break
                after.append(next_scan)
                current_index = next_scan.index
                current_time = next_scan.scan_time
        else:
            raise ValueError("One of `index_interval` or `rt_interval` must be provided")
        return before, after

    def _compute_smoothing_weights(self, scans, mean, sigma=0.025):
        sigma_sqrd_2 = (2 * sigma ** 2)
        time_array = np.array([s.scan_time for s in scans])
        weights = np.exp((-(time_array - mean) ** 2) / sigma_sqrd_2)
        return weights

    def average(self, index_interval=None, rt_interval=None, dx=0.01, weight_sigma=None):
        before, after = self._get_adjacent_scans(index_interval, rt_interval)
        scans = before + [self] + after
        arrays = []
        for scan in scans:
            if scan.is_profile:
                arrays.append(scan.arrays)
            else:
                arrays.append(scan.reprofile(dx=dx).arrays)
        if weight_sigma:
            if weight_sigma == 1:
                weight_sigma = 0.025
            weights = self._compute_smoothing_weights(
                scans, mean=self.scan_time, sigma=weight_sigma)
        else:
            weights = None
        new_arrays = average_signal(arrays, dx=dx, weights=weights)
        indices = [scan.index for scan in scans]
        return AveragedScan(
            self._data, self.source, new_arrays,
            indices, list(self.product_scans), is_profile=True,
            annotations=self._external_annotations)


class WrappedScan(Scan):
    overridable_keys = [
        "_arrays",
        "_id",
        "_title",
        "_ms_level",
        "_scan_time",
        "_precursor_information",
        "_index",
        "_is_profile",
        "_polarity",
        "_activation",
        "_acquisition_information",
        "_isolation_window",
        "_instrument_configuration"
    ]

    def __init__(self, data, source, array_data, product_scans=None, annotations=None, **overrides):
        super(WrappedScan, self).__init__(
            data, source, peak_set=None,
            deconvoluted_peak_set=None,
            annotations=annotations,
            product_scans=product_scans)
        self._arrays = RawDataArrays(*array_data)
        self._overrides = overrides
        for key, value in overrides.items():
            if not key.startswith("_"):
                key = "_" + key
            if key in self.overridable_keys:
                setattr(self, key, value)

    def clone(self):
        dup = self.__class__(
            self._data, self.source, self.arrays,
            [s.clone() for s in self.product_scans],
            annotations=self._external_annotations,
            **self._overrides)
        dup.peak_set = self.peak_set.clone()
        dup.deconvoluted_peak_set = self.deconvoluted_peak_set.clone()
        return dup


class AveragedScan(WrappedScan):
    def __init__(self, data, source, array_data, scan_indices, product_scans=None, annotations=None, **overrides):
        super(AveragedScan, self).__init__(
            data, source, array_data,
            product_scans=product_scans, annotations=annotations, **overrides)
        self.scan_indices = scan_indices

    def clone(self):
        dup = self.__class__(
            self._data, self.source, self.arrays,
            self.scan_indices,
            [s.clone() for s in self.product_scans],
            annotations=self._external_annotations,
            **self._overrides)
        dup.peak_set = self.peak_set.clone()
        dup.deconvoluted_peak_set = self.deconvoluted_peak_set.clone()
        return dup


class PrecursorInformation(object):
    """Store information relating a tandem MS scan to its precursor MS scan.

    .. note:
        The attributes prefixed with `extracted_` refer to the quantities estimated
        from the data, while those unprefixed are the values read directly from the
        data source. These values regularly do not agree. When available, the extracted
        values should be more accurate.

    Attributes
    ----------
    charge : int
        The charge reported in the source metadata
    defaulted : bool
        Whether the information in the extracted fields reflects empirical
        information or fell back on the vendor-reported values.
    extracted_charge : int
        The charge estimated from the source data
    extracted_intensity : float
        The sum of the peak heights of the extracted isotopic pattern
    extracted_neutral_mass : float
        The monoisotopic neutral mass estimated from the source data
    extracted_peak : DeconvolutedPeak
        The deconvoluted peak summarizing the precursor ion
    intensity : float
        The abundance reported in the source metadata
    mz : float
        The m/z reported in the source metadata
    orphan : bool
        Whether there was an isotopic pattern to extract in the precursor scan. Usually
        paired with `defaulted`
    peak : FittedPeak
        The peak nearest `mz`, and the starting point for estimating information
        about the precursor ion
    precursor_scan_id : str
        The id string for the precursor scan
    source : ScanIteratorBase
        Any object implementing the `ScanIteratorBase` interface to be used to look up
        the precursor scan with `precursor_scan_id`
    """
    def __init__(self, mz, intensity, charge, precursor_scan_id=None, source=None,
                 extracted_neutral_mass=0, extracted_charge=0, extracted_intensity=0,
                 peak=None, extracted_peak=None, defaulted=False, orphan=False,
                 product_scan_id=None):
        try:
            charge = int(charge)
        except Exception:
            pass
        try:
            extracted_charge = int(extracted_charge)
        except Exception:
            pass

        self.mz = mz
        self.intensity = intensity
        self.charge = charge

        self.precursor_scan_id = precursor_scan_id
        self.source = source

        self.extracted_neutral_mass = extracted_neutral_mass
        self.extracted_charge = extracted_charge
        self.extracted_intensity = extracted_intensity

        self.peak = peak
        self.extracted_peak = extracted_peak
        self.defaulted = defaulted
        self.orphan = orphan
        self.product_scan_id = product_scan_id

    def __repr__(self):
        return "PrecursorInformation(mz=%0.4f/%0.4f, intensity=%0.4f/%0.4f, charge=%r/%r, scan_id=%r)" % (
            self.mz,
            self.extracted_mz if self.extracted_neutral_mass != 0. else 0.,
            self.intensity, self.extracted_intensity or 0., self.charge,
            self.extracted_charge or 0., self.precursor_scan_id)

    def __reduce__(self):
        return self.__class__, self.__getstate__()

    def __getstate__(self):
        return (self.mz, self.intensity, self.charge, self.precursor_scan_id, None, self.extracted_neutral_mass,
                self.extracted_charge, self.extracted_intensity, self.peak, self.extracted_peak,
                self.defaulted, self.orphan, self.product_scan_id)

    def __setstate__(self, state):
        (self.mz, self.intensity, self.charge, self.precursor_scan_id, self.source, self.extracted_neutral_mass,
         self.extracted_charge, self.extracted_intensity, self.peak, self.extracted_peak,
         self.defaulted, self.orphan, self.product_scan_id) = state

    def extract(self, peak, override_charge=None):
        self.extracted_neutral_mass = peak.neutral_mass
        self.extracted_charge = int(peak.charge) if override_charge is None else override_charge
        self.extracted_intensity = peak.intensity
        self.extracted_peak = peak

    def default(self, orphan=False):
        if self.charge == ChargeNotProvided:
            warnings.warn("A precursor has been defaulted with an unknown charge state.")
            self.extracted_charge = ChargeNotProvided
            self.extracted_neutral_mass = neutral_mass(self.mz, DEFAULT_CHARGE_WHEN_NOT_RESOLVED)
            self.extracted_intensity = self.intensity
            self.defaulted = True
        else:
            self.extracted_charge = int(self.charge)
            self.extracted_neutral_mass = self.neutral_mass
            self.extracted_intensity = self.intensity
            self.defaulted = True
        if orphan:
            self.orphan = True

    @property
    def neutral_mass(self):
        if self.charge == ChargeNotProvided:
            warnings.warn("A precursor with an unknown charge state was used to compute a neutral mass.")
            return neutral_mass(self.mz, DEFAULT_CHARGE_WHEN_NOT_RESOLVED)
        return neutral_mass(self.mz, self.charge)

    @property
    def extracted_mz(self):
        if self.extracted_charge == ChargeNotProvided:
            warnings.warn("A precursor with an unknown charge state was used to compute a m/z.")
            return mass_charge_ratio(self.mz, DEFAULT_CHARGE_WHEN_NOT_RESOLVED)
        return mass_charge_ratio(self.extracted_neutral_mass, self.extracted_charge)

    @property
    def precursor(self):
        return self.source.get_scan_by_id(self.precursor_scan_id)

    @property
    def product(self):
        return self.source.get_scan_by_id(self.product_scan_id)


class ProcessedScan(ScanBase):
    def __init__(self, id, title, precursor_information,
                 ms_level, scan_time, index, peak_set,
                 deconvoluted_peak_set, polarity=None, activation=None,
                 acquisition_information=None, isolation_window=None,
                 instrument_configuration=None, product_scans=None,
                 annotations=None):
        if product_scans is None:
            product_scans = []
        if annotations is None:
            annotations = {}
        self.id = id
        self.title = title
        self.precursor_information = precursor_information
        self.ms_level = ms_level
        self.scan_time = scan_time
        self.index = index
        self.peak_set = peak_set
        self.deconvoluted_peak_set = deconvoluted_peak_set
        self.polarity = polarity
        self.activation = activation
        self.acquisition_information = acquisition_information
        self.isolation_window = isolation_window
        self.instrument_configuration = instrument_configuration
        self.product_scans = product_scans
        self.annotations = annotations

    def clear(self):
        self.peak_set = None
        self.deconvoluted_peak_set = None
        self.activation = None
        self.acquisition_information = None
        self.isolation_window = None
        self.instrument_configuration = None
        self.product_scans = None

    @property
    def scan_id(self):
        return self.id

    def __iter__(self):
        return iter(self.deconvoluted_peak_set)

    def __getitem__(self, index):
        return self.deconvoluted_peak_set[index]

    def __repr__(self):
        if self.deconvoluted_peak_set is not None:
            peaks = self.deconvoluted_peak_set
        elif self.peak_set is not None:
            peaks = self.peak_set
        else:
            peaks = []

        pinfo = self.precursor_information
        if pinfo:
            pinfo_string = ", %s" % pinfo
        else:
            pinfo_string = ""

        return "ProcessedScan(id=%s, ms_level=%d, %d peaks%s)" % (
            self.id, self.ms_level, len(peaks), pinfo_string)

    def clone(self):
        dup = self.__class__(
            self.id, self.title, self.precursor_information, self.ms_level,
            self.scan_time, self.index, self.peak_set, self.deconvoluted_peak_set,
            self.polarity, self.activation, self.acquisition_information,
            self.isolation_window, self.instrument_configuration,
            list(self.product_scans), self.annotations.copy())
        return dup


class ActivationInformation(object):
    def __init__(self, method, energy, data=None):
        if data is None:
            data = dict()
        self.method = dissociation_methods.get(str(method).lower(), method)
        self.energy = energy
        self.data = data

    def __repr__(self):
        return "ActivationInformation(%r, %r%s)" % (
            str(self.method), self.energy,
            "" if not self.data else ", %r" % self.data)

    def __str__(self):
        return str(self.method)


CID = Constant("collision-induced dissociation")
HCD = Constant("beam-type collision-induced dissociation")
ETD = Constant("electron transfer dissociation")
ECD = Constant("electron capture dissociation")
UnkownDissociation = Constant("unknown dissociation")


dissociation_methods = {
    "cid": CID,
    'cad': CID,
    CID.name.lower(): CID,
    'hcd': HCD,
    HCD.name.lower(): HCD,
    'etd': ETD,
    ETD.name.lower(): ETD,
    "ecd": ECD,
    ECD.name.lower(): ECD,
    "": UnkownDissociation,
    None: UnkownDissociation
}


ActivationInformation.dissociation_methods = dissociation_methods


class IsolationWindow(namedtuple("IsolationWindow", ['lower', 'target', 'upper'])):

    @property
    def lower_bound(self):
        return self.target - self.lower

    @property
    def upper_bound(self):
        return self.target + self.upper

    def __contains__(self, x):
        return self.lower_bound <= x <= self.upper_bound

    def is_empty(self):
        if self.lower is None:
            return self.upper is None
        return self.lower == self.upper == 0.0


class ScanAcquisitionInformation(object):
    def __init__(self, combination, scan_list):
        self.combination = combination
        self.scan_list = scan_list

    def __getitem__(self, i):
        return self.scan_list[i]

    def __iter__(self):
        return iter(self.scan_list)

    def __len__(self):
        return len(self.scan_list)

    def __repr__(self):
        return "ScanAcquisitionInformation(combination=%r, scan_list=%r)" % (
            self.combination, self.scan_list)


class ScanEventInformation(object):
    def __init__(self, start_time, window_list, drift_time=None):
        self.start_time = start_time
        self.window_list = window_list or []
        self.drift_time = drift_time

    def has_ion_mobility(self):
        return self.drift_time is not None and self.drift_time > 0

    def __getitem__(self, i):
        return self.window_list[i]

    def __iter__(self):
        return iter(self.window_list)

    def __len__(self):
        return len(self.window_list)

    def __repr__(self):
        template = "ScanEventInformation(start_time={}, window_list={}{})"
        if self.has_ion_mobility():
            tail = ", drift_time={}".format(self.drift_time)
        else:
            tail = ''
        form = template.format(self.start_time, self.window_list, tail)
        return form

    def total_scan_window(self):
        low = float('inf')
        high = 0
        for window in self:
            low = min(low, window.lower)
            high = max(high, window.upper)
        return ScanWindow(low, high)


class ScanWindow(namedtuple("ScanWindow", ['lower', 'upper'])):
    def __contains__(self, i):
        return self.lower <= i <= self.upper


class IteratorFacadeBase(DataAccessProxy, ScanIterator):
    def __init__(self, source, **kwargs):
        DataAccessProxy.__init__(self, source)
        self._producer = None

    def make_iterator(self, iterator=None, grouped=True):
        if grouped:
            self._producer = self.source._scan_group_iterator(iterator)
        else:
            self._producer = self.source._single_scan_iterator(iterator)

    def _transform(self, scan_bunch):
        return scan_bunch

    def next(self):
        return self._transform(next(self._producer))


class ComponentGroup(object):
    def __init__(self, type, members, order):
        self.type = type
        self.members = list(members)
        self.order = int(order)

    def __repr__(self):
        t = "{s.__class__.__name__}({s.type!r}, {s.members}, order={s.order})"
        return t.format(s=self)

    def __getitem__(self, i):
        return self.members[i]

    def __setitem__(self, i, v):
        self.members[i] = v

    def add(self, v):
        self.members.append(v)

    def __len__(self):
        return len(self.members)


class InstrumentInformation(object):
    def __init__(self, id, groups):
        self.id = id
        self.groups = sorted(groups, key=lambda x: x.order)
        self.analyzers = []

        for group in self.groups:
            if group.type == 'analyzer':
                self.analyzers.extend(group)

    def __getitem__(self, i):
        return self.groups[i]

    def __len__(self):
        return len(self.grou)

    def __repr__(self):
        return "{self.__class__.__name__}({self.id!r}, {self.groups})".format(
            self=self)
