"""A collection of common base classes for types that
load data for :class:`~.Scan` objects.
"""
import abc
import logging
import os
from typing import Any, Dict, Hashable, Iterator, List, Optional, Tuple, Union, Generic, TypeVar

from weakref import WeakValueDictionary

from six import string_types as basestring

import numpy as np

from ms_deisotope.utils import add_metaclass

from ms_deisotope.utils import Constant
from ms_deisotope.data_source.metadata.file_information import FileInformation
from ms_deisotope.data_source.metadata.activation import ActivationInformation
from ms_deisotope.data_source.metadata.instrument_components import InstrumentInformation
from ms_deisotope.data_source.metadata.data_transformation import DataProcessingInformation
from ms_deisotope.data_source.metadata.software import Software
from ms_deisotope.data_source.metadata.scan_traits import IsolationWindow, ScanAcquisitionInformation
from ms_deisotope.data_source._compression import MaybeFastRandomAccess

from .base import PrecursorInformation, ScanBase, ScanBunch
from .scan import Scan
from .scan_iterator import (
    _ScanIteratorImplBase,
    _SingleScanIteratorImpl,
    _GroupedScanIteratorImpl,
    _FakeGroupedScanIteratorImpl,
    _InterleavedGroupedScanIteratorImpl,
    ITERATION_MODE_GROUPED,
    ITERATION_MODE_SINGLE,
    MSEIterator)


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


DataPtrType = TypeVar("DataPtrType")
ScanType = TypeVar("ScanType", bound=ScanBase)


@add_metaclass(abc.ABCMeta)
class ScanDataSource(Generic[DataPtrType, ScanType]):
    """An Abstract Base Class describing an object
    which can provide a consistent set of accessors
    for a particular format of mass spectrometry data.

    Data files come in many shapes and sizes, with different
    underlying structures. This class provides an API that
    should make features as consistent as possible to clients
    of :class:`Scan` objects.
    """

    def _make_scan(self, data: DataPtrType) -> ScanType:
        return Scan(data, self)

    def _pick_peaks_vendor(self, scan: DataPtrType, *args, **kwargs):
        """Invoke the underlying data access library's peak picking procedure.

        Not available for open format readers, where behavior will default to
        the :mod:`ms_peak_picker` algorithm.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        peak_set: :class:`ms_peak_picker.PeakIndex`

        Raises
        ------
        NotImplementedError:
            When there is no method available for the given scan and/or data source
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def _scan_arrays(self, scan: DataPtrType) -> Union[Tuple[np.ndarray, np.ndarray], Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]]:
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
    def _precursor_information(self, scan: DataPtrType) -> Optional[PrecursorInformation]:
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
    def _scan_title(self, scan: DataPtrType) -> str:
        """Returns a verbose name for this scan, if one
        was stored in the file. Usually includes both the
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
    def _scan_id(self, scan: DataPtrType) -> str:
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
    def _scan_index(self, scan: DataPtrType) -> int:
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
    def _ms_level(self, scan: DataPtrType) -> int:
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
    def _scan_time(self, scan: DataPtrType) -> float:
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
    def _is_profile(self, scan: DataPtrType) -> bool:
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
    def _polarity(self, scan: DataPtrType) -> int:
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
    def _activation(self, scan: DataPtrType) -> Optional[ActivationInformation]:
        """Returns information about the activation method used to
        produce this scan, if any.

        Returns :const:`None` for MS1 scans

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a :class:`dict`

        Returns
        -------
        ActivationInformation
        """
        raise NotImplementedError()

    def _acquisition_information(self, scan: DataPtrType) -> Optional[ScanAcquisitionInformation]:
        return None

    def _isolation_window(self, scan: DataPtrType) -> Optional[IsolationWindow]:
        return None

    def _instrument_configuration(self, scan: DataPtrType) -> Optional[InstrumentInformation]:
        return None

    def _annotations(self, scan: DataPtrType) -> Dict[str, Any]:
        return dict()

    @property
    def source_file_name(self) -> Optional[str]:
        """Return the name of the file that backs this data source, if available.

        Returns
        -------
        :class:`str` or :const:`None`
        """
        try:
            file_ = self.source_file
        except AttributeError:
            return None
        if isinstance(file_, (basestring, os.PathLike)):
            return file_
        try:
            name = file_.name
        except AttributeError:
            return None
        return name

    def close(self):
        """Close the underlying scan data stream, which may be a file or other
        system resource.

        A closed data source may not be able to serve data requests, but not all
        :class:`ScanDataSource` implementations require the data stream be open
        for all operations.
        """

    def _dispose(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, type, exc_value, traceback):
        self.close()


@add_metaclass(abc.ABCMeta)
class ScanIterator(ScanDataSource[DataPtrType, ScanType]):
    """An Abstract Base Class that extends ScanDataSource
    with additional requirements that enable clients of the
    class to treat the object as an iterator over the underlying
    data file.

    Attributes
    ----------
    iteration_mode: str
        A string denoting :const:`~.ITERATION_MODE_GROUPED` or :const:`~.ITERATION_MODE_SINGLE`
        that controls whether :class:`~.ScanBunch` or :class:`~.Scan` are produced
        by iteration.
    """

    iteration_mode = ITERATION_MODE_GROUPED

    def has_ms1_scans(self) -> bool:
        """Checks if this :class:`ScanDataSource` contains MS1 spectra.

        Returns
        -------
        :class:`bool` or :const:`None`
            Returns a boolean value if the presence of MS1 scans is known for certain, or :const:`None`
            if it cannot be determined in the case of missing metadata.
        """
        return True

    def has_msn_scans(self) -> bool:
        """Checks if this :class:`ScanDataSource` contains MSn spectra.

        Returns
        -------
        :class:`bool` or :const:`None`
            Returns a boolean value if the presence of MSn scans is known for certain, or :const:`None`
            if it cannot be determined in the case of missing metadata.
        """
        return True

    @abc.abstractmethod
    def next(self) -> Union[ScanType, ScanBunch]:
        """Advance the iterator, fetching the next :class:`~.ScanBunch` or :class:`~.ScanBase`
        depending upon iteration strategy.

        Returns
        -------
        :class:`~.ScanBunch` or :class:`~.ScanBase`
        """
        raise NotImplementedError()

    def __next__(self) -> Union[ScanType, ScanBunch]:
        """Advance the iterator, fetching the next :class:`~.ScanBunch` or :class:`~.ScanBase`
        depending upon iteration strategy.

        Returns
        -------
        :class:`~.ScanBunch` or :class:`~.ScanBase`
        """
        return self.next()

    def __iter__(self) -> Iterator[Union[ScanType, ScanBunch]]:
        return self

    def reset(self):
        """Reset the iterator, if possible, and clear any caches.
        """
        raise NotImplementedError()

    def _dispose(self):
        extant = len(self.scan_cache)
        if extant > 0:
            logger.info("Disposing of %s with %d extant scans attached to it.", self, extant)
            for _key, value in list(self.scan_cache.items()):
                value.clear()
        self.scan_cache.clear()
        # Break reference cycle between the iterator and the iteration strategy
        self._producer = iter([])

    @abc.abstractmethod
    def _make_default_iterator(self):
        """Set up the default iterator for the :class:`ScanIterator`.
        """
        raise NotImplementedError()

    def make_iterator(self, iterator=None, grouped=None, **kwargs) -> 'ScanIterator':
        """Configure the :class:`ScanIterator`'s behavior, selecting it's iteration strategy over
        either its default iterator or the provided ``iterator`` argument.

        Parameters
        ----------
        iterator : Iterator, optional
            The iterator to manipulate. If missing, the default
            iterator will be used.
        grouped : bool, optional
            Whether the iterator should be grouped and produce :class:`.ScanBunch` objects
            or single :class:`.Scan`. If :const:`None` is passed, :meth:`has_ms1_scans` will be
            be used instead. Defaults to :const:`None`.
        """
        if grouped is None:
            grouped = self.has_ms1_scans()

        if grouped:
            self._producer = self._scan_group_iterator(iterator, grouped, **kwargs)
            self.iteration_mode = ITERATION_MODE_GROUPED
        else:
            self._producer = self._single_scan_iterator(iterator, grouped, **kwargs)
            self.iteration_mode = ITERATION_MODE_SINGLE
        return self

    def _make_cache_key(self, scan) -> Hashable:
        return scan.id

    def _cache_scan(self, scan):
        key = self._make_cache_key(scan)
        self._scan_cache[key] = scan
        return key

    def _validate(self, scan) -> bool:
        return True

    def _single_scan_iterator(self, iterator: Iterator=None, mode=None, **kwargs) -> _SingleScanIteratorImpl:
        if iterator is None:
            iterator = self._make_default_iterator()

        impl = _SingleScanIteratorImpl.from_scan_source(iterator, self, **kwargs)
        return impl

    def _scan_group_iterator(self, iterator: Iterator=None, mode=None, **kwargs) -> Union[_InterleavedGroupedScanIteratorImpl, _FakeGroupedScanIteratorImpl]:
        if iterator is None:
            iterator = self._make_default_iterator()
        if isinstance(mode, _ScanIteratorImplBase):
            impl = mode.from_scan_source(iterator, self, **kwargs)
            return impl
        elif callable(mode):
            impl = mode(iterator, self, **kwargs)
            return impl
        elif mode == "mse":
            impl = MSEIterator.from_scan_source(iterator, self, **kwargs)
            return impl
        elif self.has_ms1_scans():
                impl = _InterleavedGroupedScanIteratorImpl.from_scan_source(iterator, self, **kwargs)
        else:
            impl = _FakeGroupedScanIteratorImpl.from_scan_source(iterator, self, **kwargs)
        return impl

    def _scan_cleared(self, scan):
        self.scan_cache.pop(self._make_cache_key(scan), None)

    def initialize_scan_cache(self):
        """Initialize a cache which keeps track of which :class:`~.Scan`
        objects are still in memory using a :class:`weakref.WeakValueDictionary`.

        When a scan is requested, if the scan object is found in the cache, the
        existing object is returned rather than re-read from disk.
        """
        self._scan_cache = WeakValueDictionary()

    @property
    def scan_cache(self):
        """A :class:`weakref.WeakValueDictionary` mapping used to retrieve
        scans from memory if available before re-reading them from disk.
        """
        return self._scan_cache

    @scan_cache.setter
    def scan_cache(self, value):
        self._scan_cache = value


@add_metaclass(abc.ABCMeta)
class RandomAccessScanSource(ScanIterator[DataPtrType, ScanType]):
    """An Abstract Base Class that extends ScanIterator
    with additional requirements that the implementation support
    random access to individual scans. This should be doable by unique
    identifier, sequential index, or by scan time.
    """

    @property
    def has_fast_random_access(self) -> Constant:
        """Check whether the underlying data stream supports fast random access
        or not.

        Even if the file format supports random access, it may be impractical due
        to overhead in parsing the underlying data stream, e.g. calling :meth:`gzip.GzipFile.seek`
        can force the file to be decompressed from the *beginning of the file* on each call. This
        property can be used to signal to the caller whether or not it should use a different
        strategy.

        Returns
        -------
        :class:`Constant`:
            One of :data:`~.DefinitelyNotFastRandomAccess`, :data:`~.MaybeFastRandomAccess`, or
            :data:`~.DefinitelyFastRandomAccess`. The first is a False-y value, the latter two
            will evaluate to :const:`True`
        """
        return MaybeFastRandomAccess

    @abc.abstractmethod
    def get_scan_by_id(self, scan_id: str) -> ScanType:
        """Retrieve the scan object for the specified scan id.

        If the scan object is still bound and in memory somewhere,
        a reference to that same object will be returned. Otherwise,
        a new object will be created.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        Scan
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_scan_by_time(self, time: float) -> ScanType:
        """Retrieve the scan object for the specified scan time.

        This internally calls :meth:`get_scan_by_id` which will
        use its cache.

        Parameters
        ----------
        time : float
            The time to get the nearest scan from

        Returns
        -------
        Scan
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def get_scan_by_index(self, index: int) -> ScanType:
        """Retrieve the scan object for the specified scan index.

        This internally calls :meth:`get_scan_by_id` which will
        use its cache.

        Parameters
        ----------
        index: int
            The index to get the scan for

        Returns
        -------
        Scan
        """
        raise NotImplementedError()

    @abc.abstractmethod
    def start_from_scan(self, scan_id: Optional[str]=None, rt: Optional[float]=None, index: Optional[int]=None,
                        require_ms1: bool=True, grouped=True, **kwargs) -> 'RandomAccessScanSource':
        """Reconstruct an iterator which will start from the scan matching one of ``scan_id``,
        ``rt``, or ``index``. Only one may be provided.

        After invoking this method, the iterator this object wraps will be changed to begin
        yielding scan bunchs (or single scans if ``grouped`` is ``False``).

        This method will trigger several random-access operations, making it prohibitively
        expensive for normally compressed files.

        Arguments
        ---------
        scan_id: str, optional
            Start from the scan with the specified id.
        rt: float, optional
            Start from the scan nearest to specified time (in minutes) in the run. If no
            exact match is found, the nearest scan time will be found, rounded up.
        index: int, optional
            Start from the scan with the specified index.
        require_ms1: bool, optional
            Whether the iterator must start from an MS1 scan. True by default.
        grouped: bool, optional
            whether the iterator should yield scan bunches or single scans. True by default.
        """
        raise NotImplementedError()

    def _locate_ms1_scan(self, scan: ScanType, search_range: int=150) -> Optional[ScanType]:
        i = 0
        initial_scan = scan
        if (self.has_ms1_scans() is False):
            raise IndexError('Cannot locate MS1 Scan')
        while scan.ms_level != 1 and i < search_range:
            i += 1
            if scan.index <= 0:
                break
            scan = self.get_scan_by_index(scan.index - 1)
        if scan.ms_level == 1:
            return scan
        scan = initial_scan
        i = 0
        while scan.ms_level != 1 and i < search_range:
            i += 1
            try:
                scan = self.get_scan_by_index(scan.index + 1)
            except IndexError:
                raise IndexError("Cannot locate MS1 Scan")
        return scan

    def find_previous_ms1(self, start_index: int) -> Optional[ScanType]:
        """Locate the MS1 scan preceding ``start_index``, iterating backwards through
        scans until either the first scan is reached or an MS1 scan is found.

        Returns
        -------
        :class:`~.ScanBase` or :const:`None` if not found
        """
        if self.has_ms1_scans() is False:
            return None
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

    def find_next_ms1(self, start_index: int) -> Optional[ScanType]:
        """Locate the MS1 scan following ``start_index``, iterating forwards through
        scans until either the last scan is reached or an MS1 scan is found.

        Returns
        -------
        :class:`~.ScanBase` or :const:`None` if not found
        """
        if self.has_ms1_scans() is False:
            return None
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

    @abc.abstractmethod
    def __len__(self) -> int:
        raise NotImplementedError()

    def __getitem__(self, i: int) -> Union[ScanType, List[ScanType]]:
        """Retrieve the scan object for the specified scan index.

        This internally calls :meth:`get_scan_by_index` but supports
        slicing and negative indexing.

        Parameters
        ----------
        index: int
            The index to get the scan for

        Returns
        -------
        Scan
        """
        if isinstance(i, slice):
            n = len(self)
            scans = []
            start, stop, step = i.indices(n)
            for i in range(start, stop, step):
                scans.append(self[i])
            return scans
        elif i < 0:
            n = len(self)
            i = n + i
        return self.get_scan_by_index(i)

    @property
    def time(self):
        """A indexer facade that lets you index and slice by scan time.

        Returns
        -------
        TimeIndex
        """
        return TimeIndex(self)


class TimeIndex(object):
    """A facade that translates ``[x]`` into
    scan time access, and supports slicing over
    a time range.

    """
    def __init__(self, scan_loader):
        self.scan_loader = scan_loader

    def is_sorted_by_time(self):
        lo = 0
        hi = len(self.scan_loader)

        last = 0.0
        for i in range(lo, hi):
            current = self.scan_loader[i].scan_time
            if last <= current:
                last = current
            else:
                return False
        return True

    def __len__(self):
        return len(self.scan_loader)

    def __getitem__(self, time):
        if isinstance(time, slice):
            start_scan = self.scan_loader.get_scan_by_time(time.start)
            end_scan = self.scan_loader.get_scan_by_time(time.stop)
            return self.scan_loader[start_scan.index:end_scan.index]
        else:
            return self.scan_loader.get_scan_by_time(time)


@add_metaclass(abc.ABCMeta)
class ScanFileMetadataBase(object):
    """Objects implementing this interface can describe the original source
    files, instrument configuration, and data processing parameters used to
    create the current spectral data file.

    Patterned after the provenance features of mzML that could also be mapped
    onto mzXML and other complete vendor readers.
    """

    @abc.abstractmethod
    def file_description(self):
        """Describe the file and its components, as well
        as any content types it has.

        Returns
        -------
        :class:`~.FileInformation`
        """
        return FileInformation()

    @property
    def id_format(self):
        file_desc = self.file_description()
        return file_desc.id_format

    @abc.abstractmethod
    def instrument_configuration(self) -> List[InstrumentInformation]:
        """Describe the different instrument components and configurations used
        to acquire scans in this run.

        Returns
        -------
        :class:`list` of :class:`~.InstrumentInformation`
        """
        return []

    @abc.abstractmethod
    def data_processing(self) -> List[DataProcessingInformation]:
        """Describe any preprocessing steps applied to the data described by this
        instance.

        Returns
        -------
        :class:`list` of :class:`~.DataProcessingInformation`
        """
        return []

    def software_list(self) -> List[Software]:
        """Describe any software used on the data described by this instance.

        Returns
        -------
        :class:`list` of :class:`~.Software`
        """
        return []
