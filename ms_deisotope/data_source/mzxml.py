import numpy as np
from pyteomics import mzxml
from .common import (
    PrecursorInformation, ScanDataSource, ChargeNotProvided,
    ActivationInformation)
from .xml_reader import XMLReaderBase, IndexSavingXML
from weakref import WeakValueDictionary


class _MzXMLParser(mzxml.MzXML, IndexSavingXML):
    pass


class MzXMLDataInterface(ScanDataSource):
    """Provides implementations of all of the methods needed to implement the
    :class:`ScanDataSource` for mzXML files. Not intended for direct instantiation.
    """
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
        try:
            pinfo_dict = scan['precursorMz'][0]
            precursor_scan_id = pinfo_dict['precursorScanNum']
            pinfo = PrecursorInformation(
                mz=float(pinfo_dict['precursorMz']),
                intensity=float(pinfo_dict.get('precursorIntensity', 0.0)),
                charge=int(pinfo_dict.get('precursorCharge')) if pinfo_dict.get(
                    'precursorCharge') else ChargeNotProvided,
                precursor_scan_id=precursor_scan_id,
                source=self,
                product_scan_id=self._scan_id(scan))
            return pinfo
        except KeyError:
            return None

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
        return self._scan_id(scan)

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
        return scan["num"]

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
            if self._scan_index_lookup is None:
                raise ValueError("Index Not Built")
            scan_index = self._scan_index_lookup[self._scan_id(scan)]
            return scan_index
        except KeyError:
            return -1
        except ValueError:
            return -1

    def _ms_level(self, scan):
        """Returns the degree of exponential fragmentation
        used to produce this scan.

        1 refers to a survey scan of unfragmented ions, 2
        refers to a tandem scan derived from an ms level 1
        ion, and so on.

        Parameters
        ----------
        scan : Mapping
            The underlying scan information storage,
            usually a `dict`

        Returns
        -------
        int
        """
        return int(scan['msLevel'])

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
        try:
            return scan['retentionTime']
        except KeyError:
            return None

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
        return not bool(int(scan['centroided']))

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
        if scan['polarity'] == '+':
            return 1
        else:
            return -1

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

    @staticmethod
    def prebuild_byte_offset_file(path):
        return _MzXMLParser.prebuild_byte_offset_file(path)

    def __init__(self, source_file, use_index=True):
        self.source_file = source_file
        self._source = _MzXMLParser(source_file, read_schema=True, iterative=True, use_index=use_index)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._use_index = use_index
        self._scan_index_lookup = None
        if self._use_index:
            self._build_scan_index_lookup()

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
