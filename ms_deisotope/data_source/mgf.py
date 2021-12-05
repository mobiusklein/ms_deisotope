'''MGF is a simple human-readable format for MS/MS data. It
allows storing MS/MS peak lists and exprimental parameters.

This module provides :class:`MGFLoader`, a :class:`~.RandomAccessScanSource`
implementation.

The parser is based on :mod:`pyteomics.mgf`.
'''

from pyteomics import mgf
from pyteomics.auxiliary import OffsetIndex
import numpy as np

from six import string_types as basestring

from .scan import (
    ScanFileMetadataBase, RandomAccessScanSource, ScanDataSource,
    PrecursorInformation, _FakeGroupedScanIteratorImpl,
    ChargeNotProvided)

from .metadata.file_information import (
    FileInformation, MS_MSn_Spectrum)

from ._compression import test_if_file_has_fast_random_access


class _MGFParser(mgf.IndexedMGF):

    def parse_charge(self, charge_text, list_only=False):
        '''Pyteomics _parse_charge is very general-purpose, and
        can't be sped up, so we specialize it here.'''
        try:
            if not list_only:
                return int(charge_text.replace('+', ''))
            return list(map(self.parse_charge, charge_text.split(" ")))
        except Exception:
            if '-' in charge_text:
                return int(charge_text.replace("-", '')) * -1
            raise


class _MGFMetadata(ScanFileMetadataBase):
    """Objects implementing this interface can describe the original source
    files, instrument configuration, and data processing parameters used to
    create the current spectral data file.

    Patterned after the provenance features of mzML that could also be mapped
    onto mzXML and other complete vendor readers.
    """

    def file_description(self):
        '''Describe the file and its components, as well
        as any content types it has.

        Returns
        -------
        :class:`~.FileInformation`
        '''
        finfo = FileInformation()
        finfo.add_content("centroid spectrum")
        finfo.add_content(MS_MSn_Spectrum)
        if isinstance(self.source_file, basestring):
            finfo.add_file(self.source_file)
        elif hasattr(self.source_file, 'name'):
            finfo.add_file(self.source_file.name)
        return finfo

    def instrument_configuration(self):
        '''Describe the different instrument components and configurations used
        to acquire scans in this run.

        Returns
        -------
        :class:`list` of :class:`~.InstrumentInformation`
        '''
        return super(_MGFMetadata, self).data_processing()

    def data_processing(self):
        '''Describe any preprocessing steps applied to the data described by this
        instance.

        Returns
        -------
        :class:`list` of :class:`~.DataProcessingInformation`
        '''
        return super(_MGFMetadata, self).data_processing()


class MGFInterface(ScanDataSource):
    '''Provides a basic set of widely used MASCOT Generic File (MGF)
    data accessor mechanisms. Because MGF files lack any form of standardization,
    no strong guarantees of correctness can be made.

    This dialect does not know how to use the charge column of the peak data
    section, see :class:`~.ProcessedMGFLoader`.
    '''

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
        try:
            return scan['m/z array'], scan["intensity array"]
        except KeyError:
            return np.array([]), np.array([])

    def _ms_level(self, scan):
        return 2

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
        return scan['params']["title"].strip('.')

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
        return scan['params']["title"].strip('.')

    def _scan_time(self, scan):
        try:
            return float(scan['params']['rtinseconds']) / 60.0
        except KeyError:
            return -1

    def _is_profile(self, scan):
        return False

    def _precursor_information(self, scan):
        mz, intensity = scan['params']['pepmass']
        charge = scan['params'].get('charge', [ChargeNotProvided])[0]
        pinfo = PrecursorInformation(
            mz, intensity, charge, source=self,
            product_scan_id=self._scan_id(scan),
            defaulted=True, orphan=True)
        return pinfo

    def _polarity(self, scan):
        pinfo = self._precursor_information(scan)
        if pinfo is not None:
            if pinfo.charge:
                if pinfo.charge > 0:
                    return 1
                return -1
            return 1
        return 1

    def _activation(self, scan):
        return None

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
        try:
            return self._title_to_index[self._scan_title(scan)]
        except KeyError:
            try:
                return self._title_to_index[self._scan_title(scan) + '.']
            except KeyError:
                return -1
        return -1

    def _annotations(self, scan):
        annots = dict()
        params = scan['params']
        for key, value in params.items():
            if key in ("pepmass", "charge", "title", "rtinseconds"):
                continue
            else:
                try:
                    value = float(value)
                except ValueError:
                    if value == 'None':
                        value = None
            annots[key] = value
        return annots


class MGFLoader(MGFInterface, RandomAccessScanSource, _MGFMetadata):
    """Reads scans from MASCOT Generic File (MGF) Format files. Provides both iterative
    and random access.

    .. note::
        If the file is not sorted by retention time, :meth:`get_scan_by_time` and any
        other time-based accessors will fail.

    Attributes
    ----------
    source_file: str
        Path to file to read from.
    source: pyteomics.mgf.MGFBase
        Underlying scan data source
    header: dict
        Any top-of-the-file parameters
    """
    def __init__(self, source_file, encoding='utf-8', use_index=True, **kwargs):
        self.source_file = source_file
        self.encoding = encoding
        self._use_index = use_index
        self._source = self._create_parser()
        self.initialize_scan_cache()
        self.make_iterator()
        self._title_to_index = self._prepare_index_lookup()

    @property
    def has_fast_random_access(self):
        return test_if_file_has_fast_random_access(self.source.file)

    def _prepare_index_lookup(self):
        title_to_index = dict()
        for i, key in enumerate(self.index):
            title_to_index[key] = i
        return title_to_index

    @property
    def header(self):
        '''Any top-of-the-file parameters

        Returns
        -------
        dict
        '''
        return self._source.header

    def __reduce__(self):
        return self.__class__, (self.source_file, self.encoding, self._use_index, )

    def has_msn_scans(self):
        return True

    def has_ms1_scans(self):
        return False

    def _create_parser(self):
        if self._use_index:
            return _MGFParser(self.source_file, read_charges=False,
                              convert_arrays=1, encoding=self.encoding)
        simple_reader = mgf.MGF(
            self.source_file, read_charges=False,
            convert_arrays=1, encoding=self.encoding)
        simple_reader.index = OffsetIndex()
        return simple_reader

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
        try:
            return self.scan_cache[scan_id]
        except KeyError:
            pass
        try:
            scan = self.source.get_spectrum(scan_id)
        except KeyError:
            scan = self.source.get_spectrum(scan_id + '.')
        scan = self._make_scan(scan)
        self.scan_cache[scan_id] = scan
        return scan

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
        if not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")
        id_str = self.index.from_index(index)
        return self.get_scan_by_id(id_str)

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
        if not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")

        scan_ids = tuple(self.index)
        lo = 0
        hi = len(scan_ids)

        best_match = None
        best_error = float('inf')

        if time == float('inf'):
            return self.get_scan_by_id(scan_ids[-1])

        while hi != lo:
            mid = (hi + lo) // 2
            sid = scan_ids[mid]
            scan = self.get_scan_by_id(sid)
            scan_time = scan.scan_time
            err = abs(scan_time - time)
            if err < best_error:
                best_error = err
                best_match = scan
            if scan_time == time:
                return scan
            elif (hi - lo) == 1:
                return best_match
            elif scan_time > time:
                hi = mid
            else:
                lo = mid
        if hi == 0 and not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")

    @property
    def source(self):
        '''The file parser that this reader consumes.
        '''
        return self._source

    @property
    def index(self):
        '''The byte offset index used to achieve fast random access.

        Maps :class:`~.ScanBase` IDs to the byte offsets, implying
        the order the scans reside in the file.

        Returns
        -------
        :class:`pyteomics.xml.ByteEncodingOrderedDict`
        '''
        return self.source.index

    def __len__(self):
        return len(self.index)

    def close(self):
        '''Close the underlying reader.
        '''
        self._source.close()
        self._dispose()

    def reset(self):
        """Reset the object, clearing out any existing
        state.

        This resets the underlying file iterator, then
        calls :meth:`make_iterator`, and clears the scan
        cache.
        """
        self._source.reset()
        try:
            self.source.seek(0)
        except (IOError, AttributeError):
            pass
        self.make_iterator(None)
        self.initialize_scan_cache()

    def _make_default_iterator(self):
        return iter(self._source)

    def make_iterator(self, iterator=None, grouped=False):
        """Configure the iterator's behavior.

        Parameters
        ----------
        iterator : Iterator, optional
            The iterator to manipulate. If missing, the default
            iterator will be used.
        grouped : bool, optional
            Whether the iterator should be grouped and produce :class:`.ScanBunch` objects
            or single :class:`.Scan`. Defaults to False
        """
        return super(MGFLoader, self).make_iterator(iterator, grouped)

    def _yield_from_index(self, scan_source, start):
        offset_provider = self.index
        keys = list(offset_provider.keys())
        if start is not None:
            if isinstance(start, basestring):
                try:
                    start = keys.index(start)
                except ValueError:
                    start = keys.index(start + '.')
            elif isinstance(start, int):
                start = start
            else:
                raise TypeError("Cannot start from object %r" % start)
        else:
            start = 0
        for key in keys[start:]:
            yield scan_source.get_by_id(key)

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
        if scan_id is None:
            if rt is not None:
                scan = self.get_scan_by_time(rt)
            elif index is not None:
                try:
                    scan = self.get_scan_by_index(index)
                except IndexError:
                    if index > len(self.index):
                        index = len(self.index) - 1
                    else:
                        index = 0
                    scan = self.get_scan_by_index(index)

            else:
                raise ValueError("Must provide a scan locator, one of (scan_id, rt, index)")

            scan_id = scan.id
        else:
            scan = self.get_scan_by_id(scan_id)

        # MGF files do not contain MS1 scans
        if require_ms1:
            pass

        iterator = self._yield_from_index(self._source, scan_id)
        self.make_iterator(iterator, grouped=grouped)
        return self

    def _scan_group_iterator(self, iterator=None, mode=None):
        if iterator is None:
            iterator = self._make_default_iterator()

        impl = _FakeGroupedScanIteratorImpl(
            iterator, self._make_scan, self._validate, self._cache_scan)
        return impl

    def next(self):
        return next(self._producer)

    def _validate(self, scan):
        return True
