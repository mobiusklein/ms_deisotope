'''A collection of common base classes for types that
load data for :class:`~.Scan` objects.
'''
import abc
import weakref
from weakref import WeakValueDictionary

from collections import OrderedDict

from ms_deisotope.utils import add_metaclass

from ms_deisotope.data_source.metadata.file_information import FileInformation


from .scan import Scan, ScanBase
from .scan_iterator import (
    _SingleScanIteratorImpl,
    _GroupedScanIteratorImpl)


@add_metaclass(abc.ABCMeta)
class ScanDataSource(object):
    """An Abstract Base Class describing an object
    which can provide a consistent set of accessors
    for a particular format of mass spectrometry data.

    Data files come in many shapes and sizes, with different
    underlying structures. This class provides an API that
    should make features as consistent as possible to clients
    of :class:`Scan` objects.
    """

    def _make_scan(self, data):
        return Scan(data, self)

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

    iteration_mode = 'group'

    @abc.abstractmethod
    def next(self):
        raise NotImplementedError()

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self

    def reset(self):
        '''Reset the iterator, if possible, and clear any caches.
        '''
        raise NotImplementedError()

    @abc.abstractmethod
    def _make_default_iterator(self):
        raise NotImplementedError()

    def make_iterator(self, iterator=None, grouped=True):
        """Configure the iterator's behavior.

        Parameters
        ----------
        iterator : Iterator, optional
            The iterator to manipulate. If missing, the default
            iterator will be used.
        grouped : bool, optional
            Whether the iterator should be grouped and produce :class:`.ScanBunch` objects
            or single :class:`.Scan`. Defaults to True
        """
        if grouped:
            self._producer = self._scan_group_iterator(iterator)
            self.iteration_mode = 'group'
        else:
            self._producer = self._single_scan_iterator(iterator)
            self.iteration_mode = 'single'

    def _make_cache_key(self, scan):
        return scan.id

    def _cache_scan(self, scan):
        key = self._make_cache_key(scan)
        self._scan_cache[key] = scan
        return key

    def _validate(self, scan):
        return True

    def _single_scan_iterator(self, iterator=None):
        if iterator is None:
            iterator = self._make_default_iterator()

        impl = _SingleScanIteratorImpl(
            iterator, self._make_scan, self._validate, self._cache_scan)
        return impl

    def _scan_group_iterator(self, iterator=None):
        if iterator is None:
            iterator = self._make_default_iterator()

        impl = _GroupedScanIteratorImpl(
            iterator, self._make_scan, self._validate, self._cache_scan)
        return impl

    def _scan_cleared(self, scan):
        self.scan_cache.pop(self._make_cache_key(scan), None)

    def initialize_scan_cache(self):
        '''Initialize a cache which keeps track of which :class:`~.Scan`
        objects are still in memory using a :class:`weakref.WeakValueDictionary`.

        When a scan is requested, if the scan object is found in the cahce, the
        existing object is returned rather than re-read from disk.
        '''
        self._scan_cache = WeakValueDictionary()

    @property
    def scan_cache(self):
        '''A :class:`weakref.WeakValueDictionary` mapping used to retrieve
        scans from memory if available before re-reading them from disk.
        '''
        return self._scan_cache

    @scan_cache.setter
    def scan_cache(self, value):
        self._scan_cache = value


@add_metaclass(abc.ABCMeta)
class RandomAccessScanSource(ScanIterator):
    """An Abstract Base Class that extends ScanIterator
    with additional requirements that the implementation support
    random access to individual scans. This should be doable by unique
    identifier, sequential index, or by scan time.
    """

    @abc.abstractmethod
    def get_scan_by_id(self, scan_id):
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
    def get_scan_by_time(self, time):
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
    def get_scan_by_index(self, index):
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
    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        '''Reconstruct an iterator which will start from the scan matching one of ``scan_id``,
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
        '''
        raise NotImplementedError()

    def has_ms1_scans(self):
        return True

    def has_msn_scans(self):
        return True

    def _locate_ms1_scan(self, scan, search_range=150):
        i = 0
        initial_scan = scan
        if not self.has_ms1_scans():
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

    def find_previous_ms1(self, start_index):
        if not self.has_ms1_scans():
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

    def find_next_ms1(self, start_index):
        if not self.has_ms1_scans():
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
    def __len__(self):
        raise NotImplementedError()

    def __getitem__(self, i):
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


class DetachedAccessError(Exception):
    pass


class DataAccessProxy(object):  # pragma: no cover
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
        '''Describe the file and its components, as well
        as any content types it has.

        Returns
        -------
        :class:`~.FileInformation`
        '''
        return FileInformation()

    @abc.abstractmethod
    def instrument_configuration(self):
        return []

    @abc.abstractmethod
    def data_processing(self):
        return []


class IteratorFacadeBase(DataAccessProxy, ScanIterator):  # pragma: no cover
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


class ScanProxyContext(object):
    """A memory-conserving wrapper around an existing :class:`RandomAccessScanSource`
    object for serving :class:`ScanProxy` objects.

    Attributes
    ----------
    cache : :class:`collections.OrderedDict`
        A strong-reference maintaining cache from scan ID to :class:`ScanBase`
    cache_size : int
        The number of scans to keep strong references to.
    source : :class:`RandomAccessScanSource`
        The source to load scans from.
    """

    def __init__(self, source, cache_size=2 ** 10):
        self.source = source
        self.cache_size = cache_size
        self.cache = OrderedDict()

    def get_scan_by_id(self, scan_id):
        '''Retrieve a real scan by its identifier.

        This method first checks the in-memory cache for the scan identifier
        and returns it if already loaded. If not, the scan is retrieved from
        :attr:`source`, usually from disk and saved to the cache using
        :meth:`_save_scan`. This may cause the cache eviction.

        Parameters
        ----------
        scan_id: :class:`str`
            The scan to retrieve

        Returns
        -------
        :class:`~.Scan`
        '''
        if scan_id in self.cache:
            return self.cache[scan_id]
        scan = self.source.get_scan_by_id(scan_id)
        self._save_scan(scan_id, scan)
        return scan

    def _save_scan(self, scan_id, scan):
        if len(self.cache) > self.cache_size:
            _, evicted_scan = self.cache.popitem()
            self.source._scan_cleared(evicted_scan)

        self.cache[scan_id] = scan

    def __call__(self, scan_id):
        """Create a proxy for the scan referenced by scan_id

        Parameters
        ----------
        scan_id: :class:`str`

        Returns
        -------
        :class:`ScanProxy`
        """
        return ScanProxy(scan_id, self)

    def clear(self):
        '''Clear the reference cache.

        If the scans that were previously held in the cache were not
        strong-referenced anywhere else, any extant proxies which depend
        upon this context will have to reload their spectra from disk.
        '''
        self.cache.clear()


class proxyproperty(object):
    '''An descriptor that integrates with :class:`ScanProxy` to retrieve attributes
    from a wrapped :class:`~.ScanBase` object referenced by :class:`ScanProxy`,
    potentially loading the scan from disk if it has been purged from the memory
    cache.

    Attributes
    ----------
    name: :class:`str`
        The name of the attribute to retrieve from the wrapped scan
    '''

    def __init__(self, name, caching=False):
        self.name = name
        self.caching = caching

    def __get__(self, instance, owner):
        if self.caching:
            is_null_slot = "_%s_null" % self.name
            try:
                is_null = getattr(instance, is_null_slot)
            except AttributeError:
                setattr(instance, is_null_slot, None)
                is_null = None
            try:
                if is_null:
                    return None
                value = getattr(instance, "_" + self.name)
                return value
            except AttributeError:
                value = getattr(instance.scan, self.name)
                if is_null is None:
                    if value is None:
                        setattr(instance, is_null_slot, True)
                    else:
                        setattr(instance, is_null_slot, False)
                setattr(instance, "_" + self.name, value)
                return value
        return getattr(instance.scan, self.name)

    def __set__(self, instance, value):
        if self.caching:
            setattr(instance, "_" + self.name, value)
        else:
            raise TypeError("Cannot set attribute \"%s\"" % (self.name, ))


class ScanProxy(ScanBase):
    """A proxy for a :class:`ScanBase` object, allowing transparent access to the
    scan, potentially loading it from disk through the attached :class:`ScanProxyContext`
    :attrib:`context`.

    Attributes
    ----------
    context : :class:`ScanProxyContext`
        The context to use to retrieve scans from, whose cache layer will keep the
        scan objects alive.
    scan : :class:`weakref.ProxyType`
        A weakref proxy to the :class:`~.ScanBase` object.
    """

    _names = [
        "arrays",
        "id",
        "title",
        "ms_level",
        "scan_time",
        "index",
        "is_profile",
        "polarity",
        "activation",
        "acquisition_information",
        "isolation_window",
        "instrument_configuration",
        "annotations",
    ]

    def __init__(self, scan_id, context):
        self._target_scan_id = scan_id
        self.context = context
        self._scan = None
        self._precursor_information = None
        self._peak_set = None

    @property
    def source(self):
        '''The source of the scan data this proxy is bound to.

        Returns
        -------
        :class:`RandomAccessScanSource`
        '''
        return self.context.source

    def _clear_scan(self, *args, **kwargs):
        self._scan = None
        self._peak_set = None
        self._deconvoluted_peak_set = None

    def _require_scan(self):
        if self._scan is None:

            self._scan = weakref.proxy(
                self.context.get_scan_by_id(self._target_scan_id),
                self._clear_scan)

    @property
    def scan(self):
        '''The proxied scan.

        Accessing this property may require loading the scan's data
        from disk and/or triggering a cache eviction.
        '''
        self._require_scan()
        return self._scan

    def pick_peaks(self, *args, **kwargs):
        '''Request the proxied scan picked peaks if they are not
        already available and then cache them in memory.
        '''
        self._require_scan()
        peaks = self.peak_set
        if peaks is None:
            peaks = self.scan.pick_peaks(*args, **kwargs).peak_set
        self.peak_set = peaks
        return self

    def _resolve_sequence(self):
        deconvoluted_peak_set = self.deconvoluted_peak_set
        if deconvoluted_peak_set is not None:
            return deconvoluted_peak_set
        peak_set = self.peak_set
        if peak_set is not None:
            return peak_set
        self.pick_peaks()
        return self.peak_set

    def __getitem__(self, i):
        return self._resolve_sequence()[i]

    def __iter__(self):
        return iter(self._resolve_sequence())

    def __len__(self):
        return len(self._resolve_sequence())

    def has_peak(self, *args, **kwargs):
        '''Query the wrapped scan's peaks using it's :meth:`Scan.has_peak`
        method. If no peaks are available, this will call :meth:`pick_peaks`
        first to resolve the peak set and then query its peaks instead.

        Parameters
        ----------
        mz: float
            The m/z to search for
        error_tolerance: float
            The parts per million mass error tolerance to use

        Returns
        -------
        ms_peak_picker.FittedPeak or None
            The peak closest to the query m/z within the error tolerance window or None
            if there are no peaks satisfying the requirements

        See Also
        --------
        :meth:`Scan.has_peak`
        '''
        try:
            return self.scan.has_peak(*args, **kwargs)
        except ValueError:
            return self.pick_peaks().has_peak(*args, **kwargs)

    peak_set = proxyproperty('peak_set', True)
    deconvoluted_peak_set = proxyproperty('deconvoluted_peak_set', True)
    precursor_information = proxyproperty("precursor_information", True)

    def __repr__(self):
        template = "{self.__class__.__name__}({self._target_scan_id!r})"
        return template.format(self=self)

    @classmethod
    def _configure_proxy_attributes(cls):
        for name in cls._names:
            setattr(cls, name, proxyproperty(name))

ScanProxy._configure_proxy_attributes()
