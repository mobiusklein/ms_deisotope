'''mzML is a standard rich XML-format for raw mass spectrometry data storage.
This module provides :class:`MzMLLoader`, a :class:`~.RandomAccessScanSource`
implementation.

The parser is based on :mod:`pyteomics.mzml`.
'''

from six import string_types as basestring

import numpy as np
from pyteomics import mzml
from .common import (
    PrecursorInformation, ScanDataSource,
    ChargeNotProvided, ActivationInformation,
    ScanAcquisitionInformation, ScanEventInformation,
    ScanWindow, IsolationWindow,
    FileInformation, SourceFile, MultipleActivationInformation,
    ScanFileMetadataBase)
from .metadata.activation import (
    supplemental_energy, UnknownDissociation)
from .metadata.instrument_components import (
    InstrumentInformation, ComponentGroup, component,
    instrument_models)
from .metadata.software import Software
from .metadata import file_information
from .metadata import data_transformation
from .metadata.sample import Sample
from .xml_reader import (
    XMLReaderBase, iterparse_until,
    get_tag_attributes, _find_section, in_minutes)


class _MzMLParser(mzml.MzML):
    # we do not care about chromatograms
    _indexed_tags = {'spectrum', }

    def _handle_param(self, element, **kwargs):
        try:
            element.attrib["value"]
        except KeyError:
            element.attrib["value"] = ""
        return super(_MzMLParser, self)._handle_param(element, **kwargs)

    def _determine_array_dtype(self, info):
        dtype = None
        types = {'32-bit float': np.float32, '64-bit float': np.float64,
                 '32-bit integer': np.int32, '64-bit integer': np.int64,
                 'null-terminated ASCII string': np.uint8}
        for t, code in types.items():
            if t in info:
                dtype = code
                del info[t]
                break
        # sometimes it's under 'name'
        else:
            if 'name' in info:
                for t, code in types.items():
                    if t in info['name']:
                        dtype = code
                        info['name'].remove(t)
                        break
        return dtype


def _find_arrays(data_dict, decode=False):
    arrays = dict()
    for key, value in data_dict.items():
        if " array" in key:
            arrays[key] = value.decode() if decode else value
    return arrays


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
            decode = not self._decode_binary
        except AttributeError:
            decode = False
        arrays = _find_arrays(scan, decode=decode)
        try:
            return arrays.pop('m/z array'), arrays.pop("intensity array"), arrays
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
            # only attempt to scan if there are supposed to be MS1 scans in the file
            if self._has_ms1_scans() and self._use_index:
                last_index = self._scan_index(scan) - 1
                current_level = self._ms_level(scan)
                i = 0
                while last_index > 0 and i < 100:
                    prev_scan = self.get_scan_by_index(last_index)
                    if prev_scan.ms_level >= current_level:
                        last_index -= 1
                    else:
                        precursor_scan_id = self._scan_id(prev_scan._data)
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
            time = scan['scanList']['scan'][0]['scan start time']
            time = in_minutes(time)
            return time
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
            activation_methods = []
            for key in tuple(struct):
                if key in ActivationInformation.dissociation_methods:
                    activation_methods.append(key)
                    struct.pop(key)

            if activation_methods:
                activation = activation_methods[0]
                activation_methods = activation_methods[1:]
            else:
                activation = UnknownDissociation
            energy = struct.pop("collision energy", -1)
            if energy == -1:
                energy = struct.pop("activation energy", -1)

            supplemental_energy_ = struct.pop(supplemental_energy, -1)

            if supplemental_energy_ != -1:
                energies = [energy, supplemental_energy_]
            else:
                energies = [energy]

            if len(activation_methods) == 0:
                return ActivationInformation(activation, energy, struct)
            else:
                return MultipleActivationInformation(
                    [activation] + activation_methods, energies, struct)
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
                if target is not None:
                    return IsolationWindow.make_empty(target)
                else:
                    return None
            else:
                target = (upper - lower) / 2 + lower
                upper = upper - target
                lower = target - lower
        elif target is None:
            target = self._precursor_information(scan).mz
        if lower is None:
            lower = 0.0
        if upper is None:
            upper = 0.0
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
            scan = scan.copy()
            struct = {}
            struct['start_time'] = scan.pop('scan start time', 0)
            struct['drift_time'] = scan.pop('ion mobility drift time', 0)
            struct['injection_time'] = scan.pop("ion injection time", 0)
            windows = []
            for window in scan.pop("scanWindowList", {}).get("scanWindow", []):
                windows.append(ScanWindow(
                    window['scan window lower limit'],
                    window['scan window upper limit']))
            struct['window_list'] = windows
            scan.pop("instrumentConfigurationRef", None)
            struct['traits'] = scan
            scan_info_scan_list.append(ScanEventInformation(**struct))
        scan_info['scan_list'] = scan_info_scan_list
        return ScanAcquisitionInformation(**scan_info)

    def _instrument_configuration(self, scan):
        try:
            scan_list_struct = scan['scanList']
            reference = None
            for scan in scan_list_struct.get("scan", []):
                reference = scan.get("instrumentConfigurationRef")
                if reference is None:
                    continue
            if reference is None:
                reference = self._run_information['defaultInstrumentConfigurationRef']
            config = self._instrument_config[reference]
            return config
        except KeyError:
            return None

    def _annotations(self, scan):
        annot = dict()
        annot['base peak intensity'] = scan.get('base peak intensity')
        annot['base peak m/z'] = scan.get('base peak m/z')
        annot['lowest observed m/z'] = scan.get('lowest observed m/z')
        annot['highest observed m/z'] = scan.get('highest observed m/z')
        try:
            annot['filter string'] = scan['scanList']['scan'][0].get('filter string')
        except KeyError:
            annot['filter string'] = None
        annot['total ion current'] = scan.get('total ion current')
        return annot


class _MzMLMetadataLoader(ScanFileMetadataBase):

    def _collect_reference_groups(self):
        params = next(iterparse_until(self._source, "referenceableParamGroupList", "run"))
        self.reset()
        if params is None:
            return {}
        params = self._source._get_info_smart(params)
        param_groups = params.get("referenceableParamGroup", [])
        by_id = {
            g.pop('id'): g for g in param_groups
        }
        return by_id

    def file_description(self):
        """Read the file metadata and provenance from the ``<fileDescription>`` tag
        if it is present.

        Returns
        -------
        FileInformation
            The description of the file's contents and its sources
        """
        desc = _find_section(self._source, "fileDescription")
        contents = desc.get("fileContent", {})
        if not isinstance(contents, dict):
            contents = {k: '' for k in contents}
        fi = FileInformation(contents, [])
        for sf_data in desc.get('sourceFileList', {}).get("sourceFile", []):
            sf_data = sf_data.copy()
            sf_name = sf_data.pop('name', '')
            sf_location = sf_data.pop('location', '')
            sf_id = sf_data.pop('id', '')
            # incorrectly merged parameters may contaminate the "name" attribute
            if isinstance(sf_name, list):
                temp = sf_name
                sf_name = temp[0]
                for d in temp[1:]:
                    sf_data[d] = ''
            sf = SourceFile(
                sf_name, sf_location, sf_id,
                parameters=sf_data)
            fi.add_file(sf)
        return fi

    def _convert_instrument(self, configuration):
        group_collection = []
        if "componentList" in configuration:
            for category, groups in configuration['componentList'].items():
                if category == 'count':
                    continue
                for group in groups:
                    group = group.copy()
                    order = group.pop('order')
                    parts = [component(key) for key in group]
                    group_collection.append(ComponentGroup(category, parts, order))
        conf_id = configuration['id']
        serial_number = configuration.get("instrument serial number")
        potential_models = [k for k in configuration if instrument_models.get(k) is not None]
        if potential_models:
            instrument_model = potential_models[0]
        else:
            instrument_model = None
        config = InstrumentInformation(conf_id, group_collection, serial_number=serial_number, model=instrument_model)
        return config

    def instrument_configuration(self):
        """Read the instrument configurations settings from the
        ``<instrumentConfigurationList>``

        Returns
        -------
        list of InstrumentConfiguration
            A list of different instrument states that scans may be acquired under
        """
        instrument_info_list = _find_section(self._source, "instrumentConfigurationList").get(
            "instrumentConfiguration", [])
        out = []
        for instrument_info in instrument_info_list:
            if "referenceableParamGroupRef" in instrument_info:
                reference_params = self._collect_reference_groups()
                for group in instrument_info.pop('referenceableParamGroupRef', []):
                    param_group = reference_params[group['ref']]
                    instrument_info.update(param_group)
            out.append(self._convert_instrument(instrument_info))
        return out

    def software_list(self):
        softwares = _find_section(self._source, "softwareList")
        software_list = []
        for software in softwares.get('software', []):
            software_list.append(
                Software(**software))
        return software_list

    def data_processing(self):
        data_processing_list = _find_section(self._source, "dataProcessingList")
        processing_list = []
        for processing_group in data_processing_list.get("dataProcessing", []):
            dpinfo = data_transformation.DataProcessingInformation()
            dpinfo.id = processing_group['id']
            for method_group in processing_group['processingMethod']:
                order = int(method_group.pop("order", 0))
                software_id = method_group.pop("softwareRef", '')
                method = data_transformation.ProcessingMethod(order=order, software_id=software_id)
                for key, value in method_group.items():
                    if key == 'name':
                        key = value
                        value = ''
                    method.add(key, value)
                dpinfo.methods.append(method)
            processing_list.append(dpinfo)
        return processing_list

    def samples(self):
        """Describe the sample(s) used to generate the mass spectrometry
        data contained in this file.

        Returns
        -------
        :class:`list` of :class:`~.Sample`
        """
        sample_list = _find_section(self._source, "sampleList")
        result = []
        for sample_ in sample_list.get("sample", []):
            name = sample_.pop("sampleName", None)
            sample_.setdefault("name", name)
            result.append(Sample(**sample_))
        return result

    def _get_run_attributes(self):
        return get_tag_attributes(self.source, "run")


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

    _parser_cls = _MzMLParser


    def __init__(self, source_file, use_index=True, decode_binary=True, **kwargs):
        self.source_file = source_file
        self._source = self._parser_cls(source_file, read_schema=True, iterative=True,
                                        huge_tree=True, decode_binary=decode_binary,
                                        use_index=use_index)
        self.initialize_scan_cache()
        self._use_index = use_index
        self._decode_binary = decode_binary
        self._run_information = self._get_run_attributes()
        self._instrument_config = {
            k.id: k for k in self.instrument_configuration()
        }
        self._file_description = self.file_description()
        self.make_iterator()

    def _has_msn_scans(self):
        return file_information.MS_MSn_Spectrum in self._file_description

    def _has_ms1_scans(self):
        return file_information.MS_MS1_Spectrum in self._file_description

    def has_msn_scans(self):
        return self._has_msn_scans()

    def has_ms1_scans(self):
        return self._has_ms1_scans()

    @property
    def index(self):
        return self._source.index['spectrum']

    def _validate(self, scan):
        return "m/z array" in scan._data

    def _yield_from_index(self, scan_source, start):
        offset_provider = scan_source._offset_index['spectrum']
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
            yield scan_source.get_by_id(key)

    def __reduce__(self):
        return self.__class__, (self.source_file, self._use_index, self._decode_binary)
