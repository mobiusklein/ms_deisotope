import numpy as np
from pyteomics import mzml
from .common import (
    PrecursorInformation, ScanDataSource,
    ChargeNotProvided, ActivationInformation,
    ScanAcquisitionInformation, ScanEventInformation,
    ScanWindow, IsolationWindow)
from weakref import WeakValueDictionary
from .xml_reader import XMLReaderBase, IndexSavingXML


class _MzMLParser(mzml.MzML, IndexSavingXML):
    pass


class MzMLDataInterface(ScanDataSource):
    """Provides implementations of all of the methods needed to implement the
    :class:`ScanDataSource` for mzML files. Not intended for direct instantiation.
    """

    def _stray_cvs(self, scan):
        return scan.get("name", [])

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

    def _get_selected_ion(self, scan):
        pinfo_dict = scan["precursorList"]['precursor'][0]["selectedIonList"]['selectedIon'][0]
        return pinfo_dict

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
        pinfo_dict = self._get_selected_ion(scan)
        try:
            precursor_scan_id = scan["precursorList"]['precursor'][0]['spectrumRef']
        except KeyError:
            precursor_scan_id = None
            last_index = self._scan_index(scan) - 1
            current_level = self._ms_level(scan)
            i = 0
            while last_index > 0 and i < 100:
                prev_scan = self.get_scan_by_index(last_index)
                if prev_scan.ms_level >= current_level:
                    last_index -= 1
                else:
                    precursor_scan_id = self._scan_id(prev_scan)
                    break
                i += 1
        pinfo = PrecursorInformation(
            mz=pinfo_dict['selected ion m/z'],
            intensity=pinfo_dict.get('peak intensity', 0.0),
            charge=pinfo_dict.get('charge state', ChargeNotProvided),
            precursor_scan_id=precursor_scan_id,
            source=self,
            product_scan_id=self._scan_id(scan))
        return pinfo

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
        try:
            return scan["spectrum title"]
        except KeyError:
            return scan["id"]

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
        return scan["id"]

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
        return scan['index']

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
        return int(scan['ms level'])

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
            return scan['scanList']['scan'][0]['scan start time']
        except KeyError:
            return 0.0

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
        return "profile spectrum" in scan or "profile spectrum" in self._stray_cvs(scan)

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
        if "positive scan" in scan or "positive scan" in self._stray_cvs(scan):
            return 1
        elif "negative scan" in scan or "negative scan" in self._stray_cvs(scan):
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
            struct = dict(scan['precursorList']['precursor'][0]['activation'])
            activation = None
            for key in tuple(struct):
                if key in ActivationInformation.dissociation_methods:
                    activation = key
                    struct.pop(key)
                    break
            else:
                activation = "unknown dissociation method"
            energy = struct.pop("collision energy", -1)
            return ActivationInformation(activation, energy, struct)
        except KeyError:
            return None

    def _isolation_window(self, scan):
        try:
            struct = dict(scan['precursorList']['precursor'][0]['isolationWindow'])
        except KeyError:
            return None
        lower = struct.get("isolation window lower offset")
        target = struct.get('isolation window target m/z')
        upper = struct.get("isolation window upper offset")
        if lower is upper is target is None:
            upper = struct.get("isolation window upper limit")
            lower = struct.get("isolation window lower limit")
            if upper is lower is None:
                return None
            else:
                target = (upper - lower) / 2 + lower
                upper = upper - target
                lower = target - lower
        elif lower is upper is None:
            lower = upper = 0.0
        return IsolationWindow(lower, target, upper)

    def _acquisition_information(self, scan):
        scan_info = {}
        scan_list_struct = scan['scanList']
        combination = "unknown"
        if "no combination" in scan_list_struct:
            combination = "no combination"
        elif "sum of spectra" in scan_list_struct:
            combination = "sum of spectra"
        elif "median of spectra" in scan_list_struct:
            combination = "median of spectra"
        elif "mean of spectra" in scan_list_struct:
            combination = "mean of spectra"
        scan_info['combination'] = combination
        scan_info_scan_list = []
        for scan in scan_list_struct.get("scan", []):
            struct = {}
            try:
                struct['start_time'] = scan['scan start time']
            except KeyError:
                struct['start_time'] = 0
            try:
                struct['drift_time'] = scan['ion mobility drift time']
            except KeyError:
                struct['drift_time'] = 0
            windows = []
            for window in scan.get("scanWindowList", {}).get("scanWindow", []):
                windows.append(ScanWindow(
                    window['scan window lower limit'],
                    window['scan window upper limit']))
            struct['window_list'] = windows
            scan_info_scan_list.append(ScanEventInformation(**struct))
        scan_info['scan_list'] = scan_info_scan_list
        return ScanAcquisitionInformation(**scan_info)


def _find_section(source, section):
    value = next(source.iterfind(section))
    source.reset()
    return value


class _MzMLMetadataLoader(object):
    def file_description(self):
        return _find_section(self._source, "fileDescription")

    def instrument_configuration(self):
        return _find_section(self._source, "instrumentConfigurationList")

    def data_processing(self):
        return _find_section(self._source, "dataProcessingList")

    def samples(self):
        return _find_section(self._source, "sampleList")


class MzMLLoader(MzMLDataInterface, XMLReaderBase, _MzMLMetadataLoader):
    """Reads scans from PSI-HUPO mzML XML files. Provides both iterative and
    random access.

    Attributes
    ----------
    source_file: str
        Path to file to read from.
    source: pyteomics.mzml.MzML
        Underlying scan data source
    """

    @staticmethod
    def prebuild_byte_offset_file(path):
        return _MzMLParser.prebuild_byte_offset_file(path)

    def __init__(self, source_file, use_index=True):
        self.source_file = source_file
        self._source = _MzMLParser(source_file, read_schema=True, iterative=True, use_index=use_index)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._use_index = use_index

    def _validate(self, scan):
        return "m/z array" in scan._data

    def _yield_from_index(self, scan_source, start):
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
            yield scan_source.get_by_id(key)
