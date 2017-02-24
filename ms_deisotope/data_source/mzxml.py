import numpy as np
from pyteomics import mzxml
from .common import (
    PrecursorInformation, ScanDataSource, ChargeNotProvided,
    ActivationInformation)
from .xml_reader import XMLReaderBase
from weakref import WeakValueDictionary


class MzXMLDataInterface(ScanDataSource):
    """Provides implementations of all of the methods needed to implement the
    :class:`ScanDataSource` for mzXML files. Not intended for direct instantiation.
    """
    def _scan_arrays(self, scan):
        try:
            return scan['m/z array'], scan["intensity array"]
        except KeyError:
            return np.array([]), np.array([])

    def _precursor_information(self, scan):
        pinfo_dict = scan['precursorMz'][0]
        precursor_scan_id = pinfo_dict['precursorScanNum']
        pinfo = PrecursorInformation(
            mz=float(pinfo_dict['precursorMz']),
            intensity=float(pinfo_dict.get('precursorIntensity', 0.0)),
            charge=int(pinfo_dict.get('precursorCharge')) if pinfo_dict.get('precursorCharge') else ChargeNotProvided,
            precursor_scan_id=precursor_scan_id,
            source=self,
            product_scan_id=self._scan_id(scan))
        return pinfo

    def _scan_title(self, scan):
        return self._scan_id(scan)

    def _scan_id(self, scan):
        return scan["num"]

    def _scan_index(self, scan):
        try:
            if self._scan_index_lookup is None:
                raise ValueError("Index Not Built")
            scan_index = self._scan_index_lookup[self._scan_id(scan)]
            return scan_index
        except KeyError:
            return -1
        except ValueError:
            return -1

    def _ms_level(self, scan):
        return int(scan['msLevel'])

    def _scan_time(self, scan):
        try:
            return scan['retentionTime']
        except KeyError:
            return None

    def _is_profile(self, scan):
        return not bool(int(scan['centroided']))

    def _polarity(self, scan):
        if scan['polarity'] == '+':
            return 1
        else:
            return -1

    def _activation(self, scan):
        try:
            return ActivationInformation(
                scan['precursorMz'][0]['activationMethod'], scan['collisionEnergy'])
        except KeyError:
            return None


class MzXMLLoader(MzXMLDataInterface, XMLReaderBase):
    """Reads scans from mzXML files. Provides both iterative and
    random access.

    Attributes
    ----------
    source_file: str
        Path to file to read from.
    source: pyteomics.mzxml.MzXML
        Underlying scan data source
    """
    __data_interface__ = MzXMLDataInterface

    def __init__(self, source_file, use_index=True):
        self.source_file = source_file
        self._source = mzxml.MzXML(source_file, read_schema=True, iterative=True, use_index=use_index)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._use_index = use_index
        self._scan_index_lookup = None
        if self._use_index:
            self._build_scan_index_lookup()

    def __reduce__(self):
        return self.__class__, (self.source_file, self._use_index)

    def _get_scan_by_id_raw(self, scan_id):
        return self._source.get_by_id(scan_id, "num")

    def _build_scan_index_lookup(self):
        if not self._use_index:
            raise ValueError("Must index the entire file before sequential indices may computed.")
        index = dict()
        i = 0
        for scan, offset in self.index.items():
            index[scan] = i
            i += 1
        self._scan_index_lookup = index

    def _validate(self, scan):
        return "m/z array" in scan._data

    def _yield_from_index(self, scan_source, start=None):
        offset_provider = scan_source._offset_index.offsets
        keys = offset_provider.keys()
        if start is not None:
            if isinstance(start, basestring):
                start = keys.index(start)
            elif isinstance(start, int):
                start = start
            else:
                raise TypeError("Cannot start from object %r" % start)
        else:
            start = 0
        for key in keys[start:]:
            scan = scan_source.get_by_id(key, "num")
            yield scan
