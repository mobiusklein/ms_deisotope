import json
import os

from weakref import WeakValueDictionary
from .common import (
    ScanIterator, RandomAccessScanSource)
from lxml.etree import XMLSyntaxError
from pyteomics import xml


class XMLReaderBase(RandomAccessScanSource, ScanIterator):
    @property
    def index(self):
        return self._source._offset_index

    @property
    def source(self):
        return self._source

    def close(self):
        self._source.close()

    def reset(self):
        """Reset the object, clearing out any existing
        state.

        This resets the underlying file iterator, then
        calls :meth:`make_iterator`, and clears the scan
        cache.
        """
        self._source.reset()
        self.make_iterator(None)
        self._scan_cache = WeakValueDictionary()

    def make_iterator(self, iterator=None, grouped=True):
        if grouped:
            self._producer = self._scan_group_iterator(iterator)
        else:
            self._producer = self._single_scan_iterator(iterator)

    def _validate(self, scan):
        raise NotImplementedError()

    def _make_default_iterator(self):
        return iter(self._source)

    def next(self):
        try:
            return next(self._producer)
        except XMLSyntaxError:
            raise StopIteration(
                "This iterator may need to be reset by calling `reset` to continue using it after"
                " using a random-access function like `get_by_id`")

    def __next__(self):
        return self.next()

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
            return self._scan_cache[scan_id]
        except KeyError:
            packed = self._make_scan(self._get_scan_by_id_raw(scan_id))
            self._scan_cache[packed.id] = packed
            return packed

    def _get_scan_by_id_raw(self, scan_id):
        return self._source.get_by_id(scan_id)

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
        scan_ids = tuple(self.index)
        lo = 0
        hi = len(scan_ids)
        while hi != lo:
            mid = (hi + lo) // 2
            sid = scan_ids[mid]
            sid = sid.decode('utf-8')
            scan = self.get_scan_by_id(sid)
            if not self._validate(scan):
                sid = scan_ids[mid - 1]
                scan = self.get_scan_by_id(sid)
                if not self._validate(scan):
                    sid = scan_ids[mid - 2]
                    scan = self.get_scan_by_id(sid)

            scan_time = scan.scan_time
            if scan_time == time:
                return scan
            elif (hi - lo) == 1:
                return scan
            elif scan_time > time:
                hi = mid
            else:
                lo = mid
        if hi == 0 and not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")

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
        index_keys = tuple(self.index)
        id_bytes = index_keys[index]
        id_str = id_bytes.decode("utf-8")
        return self.get_scan_by_id(id_str)

    def _yield_from_index(self, scan_source, start):
        raise NotImplementedError()

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
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

        # We must start at an MS1 scan, so backtrack until we reach one
        if require_ms1:
            scan = self._locate_ms1_scan(scan)
            scan_id = scan.id

        iterator = self._yield_from_index(self._source, scan_id)
        self.make_iterator(iterator, grouped=grouped)
        return self

    def __repr__(self):
        return "{self.__class__.__name__}({self.source_file!r})".format(self=self)

    def __reduce__(self):
        return self.__class__, (self.source_file, self._use_index)


def save_byte_index(index, fp):
    """Write the byte offset index to the provided
    file

    Parameters
    ----------
    index : ByteEncodingOrderedDict
        The byte offset index to be saved
    fp : file
        The file to write the index to

    Returns
    -------
    file
    """
    encoded_index = dict()
    for key, offset in index.items():
        encoded_index[key.decode("utf8")] = offset
    json.dump(encoded_index, fp)
    return fp


def load_byte_index(fp):
    """Read a byte offset index from a file

    Parameters
    ----------
    fp : file
        The file to read the index from

    Returns
    -------
    ByteEncodingOrderedDict
    """
    data = json.load(fp)
    index = xml.ByteEncodingOrderedDict()
    for key, value in sorted(data.items(), key=lambda x: x[1]):
        index[key] = value
    return index


class PrebuiltOffsetIndex(xml.FlatTagSpecificXMLByteIndex):
    """An Offset Index class which just holds offsets
    and performs no extra scanning effort.

    Attributes
    ----------
    offsets : ByteEncodingOrderedDict
    """

    def __init__(self, offsets):
        self.offsets = offsets


class IndexSavingXML(xml.IndexedXML):
    """An extension to the IndexedXML type which
    adds facilities to read and write the byte offset
    index externally.
    """

    _save_byte_index_to_file = staticmethod(save_byte_index)
    _load_byte_index_from_file = staticmethod(load_byte_index)

    @property
    def _byte_offset_filename(self):
        path = self._source.name
        byte_offset_filename = os.path.splitext(path)[0] + '-byte-offsets.json'
        return byte_offset_filename

    def _check_has_byte_offset_file(self):
        path = self._byte_offset_filename
        return os.path.exists(path)

    def _read_byte_offsets(self):
        with open(self._byte_offset_filename, 'r') as f:
            index = PrebuiltOffsetIndex(self._load_byte_index_from_file(f))
            self._offset_index = index

    def _write_byte_offsets(self):
        with open(self._byte_offset_filename, 'w') as f:
            self._save_byte_index_to_file(self._offset_index, f)

    @xml._keepstate
    def _build_index(self):
        try:
            self._read_byte_offsets()
        except IOError:
            super(IndexSavingXML, self)._build_index()

    @classmethod
    def prebuild_byte_offset_file(cls, path):
        inst = cls(path, use_index=True)
        inst._write_byte_offsets()
