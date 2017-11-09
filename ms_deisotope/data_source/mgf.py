import re
from weakref import WeakValueDictionary
from collections import OrderedDict
import codecs

from pyteomics import mgf
import numpy as np

from .common import (
    Scan, RandomAccessScanSource, PrecursorInformation,
    ScanDataSource, ScanIterator, ChargeNotProvided)


class MGFInterface(ScanDataSource):
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
        return scan['params']["title"]

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
        return scan['params']["title"]

    def _scan_time(self, scan):
        return float(scan['params']['rtinseconds']) / 60.0

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
                else:
                    return -1
            else:
                return 1
        else:
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
            return self._index[self._scan_title(scan)]
        except KeyError:
            return -1


def _remove_bom(bstr):
    return bstr.replace(codecs.BOM_LE, b'').lstrip(b"\x00")


def chunk_mgf(path, encoding='latin-1', read_size=1000000):
    with open(path, 'rb') as fh:
        delim = _remove_bom(u"BEGIN IONS".encode(encoding))
        pattern = re.compile(delim)
        buff = fh.read(read_size)
        parts = pattern.split(buff)
        started_with_delim = buff.startswith(delim)
        tail = parts[-1]
        front = parts[:-1]
        i = 0
        for part in front:
            i += 1
            if part == b"":
                continue
            if i == 1:
                if started_with_delim:
                    yield delim + part
                else:
                    yield part
            else:
                yield delim + part
        running = True
        while running:
            buff = fh.read(read_size)
            if len(buff) == 0:
                running = False
                buff = tail
            else:
                buff = tail + buff
            parts = pattern.split(buff)
            tail = parts[-1]
            front = parts[:-1]
            for part in front:
                yield delim + part
        yield delim + tail


def index_mgf(path, encoding='latin-1', read_size=1000000):
    gen = chunk_mgf(path, encoding, read_size)
    i = 0
    index = OrderedDict()
    pattern = re.compile(_remove_bom(u"TITLE=".encode(encoding)) + b"([^\n]+)\n")
    for chunk in gen:
        match = pattern.search(chunk)
        if match:
            title = match.group(1)
            index[title.decode(encoding)] = i
        i += len(chunk)
    return index


class MGFLoader(MGFInterface, ScanIterator):

    def __init__(self, source_file, encoding='latin-1'):
        self.source_file = source_file
        self._index = index_mgf(source_file, encoding=encoding)
        self._source = mgf.read(source_file, read_charges=False, convert_arrays=1, encoding=encoding)
        self._scan_cache = WeakValueDictionary()
        self.make_iterator()

    @property
    def source(self):
        return self._source

    @property
    def index(self):
        return self._index

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
        try:
            self.source.seek(0)
        except (IOError, AttributeError):
            pass
        self.make_iterator(None)
        self._scan_cache = WeakValueDictionary()

    def _make_default_iterator(self):
        return iter(self._source)

    def make_iterator(self, iterator=None, grouped=True):
        self.iteration_mode = 'single'
        self._producer = self._single_scan_iterator(iterator)

    def next(self):
        return next(self._producer)

    def _validate(self, scan):
        return True
