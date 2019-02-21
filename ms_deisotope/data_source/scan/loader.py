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
        self._scan_cache = WeakValueDictionary()

    @property
    def scan_cache(self):
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

    def __init__(self, source, cache_size=24):
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
            self.cache.popitem()
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
        instance._require_scan()
        if self.caching:
            try:
                value = getattr(instance, "_" + self.name)
            except AttributeError:
                value = getattr(instance.scan, self.name)
                setattr(instance, "_" + self.name, value)
                return value
        return getattr(instance.scan, self.name)


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
        "peak_set",
        "deconvoluted_peak_set",
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
    ]

    def __init__(self, scan_id, context):
        self._target_scan_id = scan_id
        self.context = context
        self.scan = None
    
    @property
    def source(self):
        return self.context.source

    def _clear_scan(self, *args, **kwargs):
        self.scan = None

    def _require_scan(self):
        if self.scan is None:
            self.scan = weakref.proxy(
                self.context.get_scan_by_id(self._target_scan_id),
                self._clear_scan)

    precursor_information = proxyproperty("precursor_information", True)


for name in ScanProxy._names:
    setattr(ScanProxy, name, proxyproperty(name))
