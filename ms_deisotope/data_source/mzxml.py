'''**mzXML** is a standard XML-format for raw mass spectrometry data storage created
by the Institute for Systems Biology, intended to be replaced with **mzML**.
This module provides :class:`MzXMLLoader`, a :class:`~.RandomAccessScanSource`
implementation.

The parser is based on :mod:`pyteomics.mzxml`.
'''
from six import string_types as basestring

import numpy as np
from pyteomics import mzxml
from .common import (
    PrecursorInformation, ScanDataSource, ChargeNotProvided,
    ActivationInformation, IsolationWindow, ScanAcquisitionInformation,
    ScanEventInformation, ScanWindow,
    ComponentGroup, component, InstrumentInformation,
    FileInformation, ScanFileMetadataBase)
from .metadata import data_transformation, file_information
from .xml_reader import (
    XMLReaderBase, iterparse_until)


class _MzXMLParser(mzxml.MzXML):
    pass


scan_number_only_id_format = file_information.id_format("MS:1000776")


class _MzXMLMetadataLoader(ScanFileMetadataBase):
    def file_description(self):
        """Read the file provenance from the ``<parentFile>`` tags
        if any are present.

        This returns no information about the file's contents as this
        was not part of the mzXML schema

        Returns
        -------
        FileInformation
            The description of the file's  sources
        """
        file_info = map(self.source._get_info_smart, iterparse_until(self.source, "parentFile", "scan"))
        self.source.reset()
        file_info = list(file_info)
        fi = FileInformation({}, [])
        for parent in file_info:
            path = parent.get("fileName")
            if path is None:
                continue
            path = path.replace("file:///", '')
            fi.add_file(path, check=False)
        return fi

    def instrument_configuration(self):
        """Read the instrument configurations settings from the
        ``<msInstrument>`` elements.

        Returns
        -------
        list of InstrumentConfiguration
            A list of different instrument states that scans may be acquired under
        """
        instrument_configuration = map(self.source._get_info_smart, iterparse_until(
            self.source, "msInstrument", "scan"))
        self.source.reset()
        instrument_configuration = [
            self._convert_instrument(ic) for ic in instrument_configuration
        ]
        return instrument_configuration

    def _convert_instrument(self, configuration):
        try:
            detector = configuration.get('msDetector', {}).get('value')
        except AttributeError:
            detector = None
        if detector is not None:
            detector = component(detector)
        try:
            ionisation = configuration.get('msIonisation', {}).get('value')
        except AttributeError:
            ionisation = None
        if ionisation is not None:
            ionisation = component(ionisation)
        try:
            analyzer = configuration.get('msMassAnalyzer', {}).get('value')
        except AttributeError:
            analyzer = None
        if analyzer is not None:
            analyzer = component(analyzer)
        parts = [
            ComponentGroup("source", [ionisation], 1),
            ComponentGroup("analyzer", [analyzer], 2),
            ComponentGroup("detector", [detector], 3)
        ]
        return InstrumentInformation(configuration.get('msInstrumentID', 1), parts)

    def data_processing(self):
        data_processing = map(self.source._get_info_smart, iterparse_until(
            self.source, "dataProcessing", "scan"))
        self.source.reset()
        operation_groups = []
        for i, group in enumerate(data_processing):
            software_id = group.get("software", {}).get("name", "")
            method_group = data_transformation.ProcessingMethod(order=1, software_id=software_id)
            operations = group.get("processingOperation", [])
            if isinstance(operations, list):
                for op in operations:
                    method_group.add(op)
            else:
                method_group.add(operations)
            operation_groups.append(
                data_transformation.DataProcessingInformation([method_group], id=i))
        return operation_groups

    @property
    def id_format(self):
        return scan_number_only_id_format


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
            return (scan['m/z array'], scan["intensity array"])
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
            precursor_scan_id = pinfo_dict.get('precursorScanNum')
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
            return -2

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
        try:
            return not bool(int(scan['centroided']))
        except KeyError:
            return True

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
        try:
            if scan['polarity'] == '+':
                return 1
            else:
                return -1
        except KeyError:
            return None

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

    def _isolation_window(self, scan):
        try:
            pinfo_dict = scan['precursorMz'][0]
            target = float(pinfo_dict['precursorMz'])
            width = float(pinfo_dict['windowWideness'])
            lower = width / 2
            upper = width / 2
            return IsolationWindow(lower, target, upper)
        except KeyError:
            try:
                pinfo_dict = scan['precursorMz'][0]
                target = float(pinfo_dict['precursorMz'])
                return IsolationWindow.make_empty(target)
            except KeyError:
                return None

    def _instrument_configuration(self, scan):
        try:
            return self._instrument_config[scan['msInstrumentID']]
        except KeyError:
            return None

    def _acquisition_information(self, scan):
        scan_event = ScanEventInformation(
            scan['retentionTime'],
            window_list=[
                ScanWindow(scan.get("lowMz"), scan.get("highMz"))
            ])
        return ScanAcquisitionInformation("no combination", [scan_event])


class MzXMLLoader(MzXMLDataInterface, XMLReaderBase, _MzXMLMetadataLoader):
    """Reads scans from mzXML files. Provides both iterative and
    random access.

    Attributes
    ----------
    source_file: str
        Path to file to read from.
    source: pyteomics.mzxml.MzXML
        Underlying scan data source
    """

    _parser_cls = _MzXMLParser


    def __init__(self, source_file, use_index=True, **kwargs):
        self.source_file = source_file
        self._source = _MzXMLParser(source_file, read_schema=True, iterative=True,
                                    huge_tree=True, use_index=use_index)
        self.initialize_scan_cache()
        self._use_index = use_index
        self._scan_index_lookup = None
        if self._use_index:
            self._build_scan_index_lookup()
        self._instrument_config = {
            k.id: k for k in self.instrument_configuration()
        }
        self.reset()
        self.make_iterator()

    @property
    def index(self):
        return self._source.index['scan']

    def _get_scan_by_id_raw(self, scan_id):
        return self._source.get_by_id(scan_id, "num")

    def _build_scan_index_lookup(self):
        if not self._use_index:
            raise ValueError("Must index the entire file before sequential indices may computed.")
        index = dict()
        i = 0
        for scan, _offset in self.index.items():
            index[scan] = i
            i += 1
        self._scan_index_lookup = index

    def _validate(self, scan):
        return "m/z array" in scan._data

    def _yield_from_index(self, scan_source, start=None):
        offset_provider = scan_source._offset_index['scan']
        keys = list(offset_provider.keys())
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
