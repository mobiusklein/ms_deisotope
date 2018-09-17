import os
from collections import OrderedDict
try:
    from collections import Sequence, Mapping
except ImportError:
    from collections.abc import Sequence, Mapping
from uuid import uuid4
import warnings

import numpy as np

from six import string_types as basestring

from ms_peak_picker import PeakIndex, PeakSet, FittedPeak


try:
    from psims.mzml import writer
except ImportError:
    writer = None

try:
    WindowsError
    on_windows = True
except NameError:
    on_windows = False

from .common import ScanSerializerBase, ScanDeserializerBase
from .text_utils import (envelopes_to_array, decode_envelopes)
from ms_deisotope import peak_set
from ms_deisotope.utils import Base
from ms_deisotope.averagine import neutral_mass
from ms_deisotope.data_source.common import PrecursorInformation, ScanBunch, ChargeNotProvided
from ms_deisotope.data_source.metadata import data_transformation
from ms_deisotope.data_source.mzml import MzMLLoader
from ms_deisotope.feature_map import ExtendedScanIndex


class SpectrumDescription(Sequence):

    def __init__(self, attribs=None):
        self.descriptors = list(attribs or [])

    def __getitem__(self, i):
        if isinstance(i, int):
            return self.descriptors[i]
        else:
            for d in self:
                if i == d.get('name'):
                    return d.get('value')
            raise KeyError(i)

    def __len__(self):
        return len(self.descriptors)

    def append(self, desc):
        i = len(self)
        self.descriptors.append(desc)
        return i

    def __repr__(self):
        return "{self.__class__.__name__}({self.descriptors})".format(self=self)

    @classmethod
    def from_peak_set(cls, peak_list):
        descriptors = cls()
        try:
            base_peak = max(peak_list, key=lambda x: x.intensity)
        except ValueError:
            base_peak = None
        descriptors.append({
            "name": "base peak m/z",
            "value": base_peak.mz if base_peak else 0,
        })
        descriptors.append({
            "name": "base peak intensity",
            "value": base_peak.intensity if base_peak else 0,
            "unit_name": writer.DEFAULT_INTENSITY_UNIT
        })
        descriptors.append({
            "name": "total ion current",
            "value": sum(p.intensity for p in peak_list),
        })
        peaks_mz_order = sorted(peak_list, key=lambda x: x.mz)
        try:
            descriptors.append({
                "name": "lowest observed m/z",
                "value": peaks_mz_order[0].mz
            })
            descriptors.append({
                "name": "highest observed m/z",
                "value": peaks_mz_order[-1].mz
            })
        except IndexError:
            pass
        return descriptors

    @classmethod
    def from_arrays(cls, arrays):
        descriptors = cls()
        try:
            base_peak_i = np.argmax(arrays.intensity)
        except ValueError:
            base_peak_i = None

        descriptors.append({
            "name": "base peak m/z",
            "value": arrays.mz[base_peak_i] if base_peak_i else 0
        })
        descriptors.append({
            "name": "base peak intensity",
            "value": arrays.intensity[base_peak_i] if base_peak_i else 0,
            "unit_name": writer.DEFAULT_INTENSITY_UNIT
        })
        descriptors.append({
            "name": "total ion current",
            "value": arrays.intensity.sum(),
        })
        try:
            descriptors.append({
                "name": "lowest observed m/z",
                "value": arrays.mz[0]
            })
            descriptors.append({
                "name": "highest observed m/z",
                "value": arrays.mz[-1]
            })
        except IndexError:
            pass
        return descriptors


class SampleRun(Base):

    def __init__(self, name, uuid, completed=True, sample_type=None, **kwargs):
        self.name = name
        self.uuid = uuid
        self.sample_type = sample_type
        self.completed = completed
        self.parameters = kwargs


class MzMLSerializer(ScanSerializerBase):

    """Write :mod:`ms_deisotope` data structures to a file in mzML format.

    Attributes
    ----------
    base_peak_chromatogram_tracker : :class:`OrderedDict`
        Accumulated mapping of scan time to base peak intensity. This is
        used to write the *base peak chromatogram*.
    chromatogram_queue : :class:`list`
        Accumulate chromatogram data structures which will be written out
        after all spectra have been written to file.
    compression : :class:`str`
        The compression type to use for binary data arrays. Should be one of
        :obj:`"zlib"`, :obj:`"none"`, or :obj:`None`
    data_encoding : :class:`dict` or :class:`int` or :obj:`numpy.dtype` or :class:`str`
        The encoding specification to specify the binary encoding of numeric data arrays
        that is passed to :meth:`~.MzMLWriter.write_spectrum` and related methods.
    data_processing_list : :class:`list`
        List of packaged :class:`~.DataProcessingInformation` to write out
    deconvoluted : bool
        Indicates whether the translation should include extra deconvolution information
    file_contents_list : :class:`list`
        List of terms to include in the :obj:`<fileContents>` tag
    handle : file-like
        The file-like object being written to
    indexer : :class:`~.ExtendedScanIndex`
        The external index builder
    instrument_configuration_list : :class:`list`
        List of packaged :class:`~.InstrumentInformation` to write out
    n_spectra : int
        The number of spectra to provide a size for in the :obj:`<spectrumList>`
    processing_parameters : :class:`list`
        List of additional terms to include in a newly created :class:`~.DataProcessingInformation`
    sample_list : :class:`list`
        List of :class:`~.SampleRun` objects to write out
    sample_name : :class:`str`
        Default sample name
    sample_run : :class:`~.SampleRun`
        Description
    software_list : :class:`list`
        List of packaged :class:`~.Software` objects to write out
    source_file_list : :class:`list`
        List of packaged :class:`~.SourceFile` objects to write out
    total_ion_chromatogram_tracker : :class:`OrderedDict`
        Accumulated mapping of scan time to total intensity. This is
        used to write the *total ion chromatogram*.
    writer : :class:`~psims.mzml.writer.MzMLWriter`
        The lower level writer implementation
    """

    def __init__(self, handle, n_spectra=2e4, compression=None,
                 deconvoluted=True, sample_name=None, build_extra_index=True,
                 data_encoding=None):
        if data_encoding is None:
            data_encoding = {
                writer.MZ_ARRAY: np.float64,
                writer.INTENSITY_ARRAY: np.float32,
                writer.CHARGE_ARRAY: np.int32,
            }
        if writer is None:
            raise ImportError(
                "Cannot write mzML without psims. Please install psims to use this feature.")
        if compression is None:
            compression = writer.COMPRESSION_ZLIB
        self.handle = handle
        self.writer = writer.MzMLWriter(handle)
        self.n_spectra = n_spectra
        self.compression = compression
        self.data_encoding = data_encoding
        self._has_started_writing_spectra = False

        self.writer.__enter__()
        self._run_tag = None
        self._spectrum_list_tag = None
        self._chromatogram_list_tag = None

        self.writer.controlled_vocabularies()
        self.deconvoluted = deconvoluted

        self._initialize_description_lists()
        self._init_sample(sample_name)

        self.total_ion_chromatogram_tracker = OrderedDict()
        self.base_peak_chromatogram_tracker = OrderedDict()
        self.chromatogram_queue = []

        self.indexer = None
        if build_extra_index:
            self.indexer = ExtendedScanIndex()

    def _init_sample(self, sample_name, **kwargs):
        self.sample_name = sample_name
        self.sample_run = SampleRun(name=sample_name, uuid=str(uuid4()))
        self.add_sample({
            "name": self.sample_run.name,
            "id": "sample_1",
            "params": [
                {"name": "SampleRun-UUID", "value": self.sample_run.uuid},
            ]})

    def _initialize_description_lists(self):
        self.file_contents_list = []
        self.software_list = []
        self.source_file_list = []
        self.data_processing_list = []
        self.instrument_configuration_list = []
        self.sample_list = []

        self.processing_parameters = []

    def add_instrument_configuration(self, configuration):
        """Add an :class:`~.InstrumentInformation` object to
        the output document.

        Parameters
        ----------
        configuration: :class:`~.InstrumentInformation`
            The instrument configuration to add
        """
        component_list = []
        for group in configuration.groups:
            tag = None
            if group.type == 'source':
                tag = self.writer.Source
            elif group.type == 'analyzer':
                tag = self.writer.Analyzer
            elif group.type == 'detector':
                tag = self.writer.Detector
            else:
                continue
            component_list.append(
                tag(order=group.order, params=[g.name for g in group]))
        config_element = self.writer.InstrumentConfiguration(
            configuration.id, component_list)
        self.instrument_configuration_list.append(config_element)

    def add_software(self, software_description):
        """Add a :class:`~.Software` object to the output document.

        Parameters
        ----------
        software_description : :class:`~.Software`
            The software description to add
        """
        self.software_list.append(software_description)

    def add_file_information(self, file_information):
        for key, value in file_information.contents.items():
            if value is None:
                value = ''
            self.add_file_contents({str(key): value})
        for source_file in file_information.source_files:
            self.add_source_file(source_file)

    def add_file_contents(self, file_contents):
        """Add a key to the resulting :obj:`<fileDescription>`
        of the output document.

        Parameters
        ----------
        file_contents: :class:`str` or :class:`Mapping`
            The parameter to add
        """
        self.file_contents_list.append(file_contents)

    def remove_file_contents(self, name):
        for i, content in enumerate(self.file_contents_list):
            if isinstance(content, Mapping):
                if 'name' in content:
                    content = content['name']
                elif len(content) == 1:
                    content = list(content.keys())[0]
                else:
                    continue
            if content == name:
                break
        else:
            raise KeyError(name)
        self.file_contents_list.pop(i)

    def add_source_file(self, source_file):
        """Add the :class:`~.SourceFile` to the output document

        Parameters
        ----------
        source_file : :class:`~.SourceFile`
            The source fil to add
        """
        unwrapped = {
            "name": source_file.name,
            "location": source_file.location,
            "id": source_file.id,
            "params": []
        }
        unwrapped['params'].extend([(getattr(key, 'accession', str(key)), value)
                                    for key, value in source_file.parameters.items()])
        if source_file.id_format:
            unwrapped['params'].append(str(source_file.id_format))
        if source_file.file_format:
            unwrapped['params'].append(str(source_file.file_format))
        self.source_file_list.append(unwrapped)

    def add_data_processing(self, data_processing_description):
        """Add a new :class:`~.DataProcessingInformation` or :class:`~ProcessingMethod`
        to the output document as a new :obj:`<dataProcessing>` entry describing one or
        more :obj:`<processingMethod>`s for a single referenced :class:`~.Software`
        instance.

        Parameters
        ----------
        data_processing_description : :class:`~.DataProcessingInformation` or :class:`~.ProcessingMethod`
            Data manipulation sequence to add to the document
        """
        if isinstance(data_processing_description, data_transformation.DataProcessingInformation):
            methods = []
            for method in data_processing_description:
                content = []
                for op, val in method:
                    content.append({'name': op.name, 'value': val})
                method_descr = {
                    'software_reference': method.software_id,
                    'order': method.order,
                    'params': content
                }
                methods.append(method_descr)
            payload = {
                'id': data_processing_description.id,
                'processing_methods': methods
            }
            self.data_processing_list.append(payload)
        elif isinstance(data_processing_description, data_transformation.ProcessingMethod):
            content = []
            for op, val in data_processing_description:
                content.append({"name": op.name, 'value': val})
            payload = {
                'id': "data_processing_%d" % len(self.data_processing_list),
                'processing_methods': [{
                    'software_reference': data_processing_description.software_id,
                    'order': data_processing_description.order,
                    'params': content
                }]
            }
            self.data_processing_list.append(payload)
        else:
            self.data_processing_list.append(data_processing_description)

    def add_processing_parameter(self, name, value=None):
        """Add a new processing method to the writer's own
        :obj:`<dataProcessing>` element.

        Parameters
        ----------
        name : str
            The processing technique's name
        value : obj
            The processing technique's value, if any
        """
        self.processing_parameters.append({"name": name, "value": value})

    def add_sample(self, sample):
        self.sample_list.append(sample)

    def copy_metadata_from(self, reader):
        try:
            description = reader.file_description()
            self.add_file_information(description)
        except AttributeError:
            pass

        try:
            instrument_configs = reader.instrument_configuration()
        except AttributeError:
            instrument_configs = []
        for config in instrument_configs:
            self.add_instrument_configuration(config)

        try:
            software_list = reader.software_list()
        except AttributeError:
            software_list = []
        for software in software_list:
            self.add_software(software)

        try:
            data_processing_list = reader.data_processing()
        except AttributeError:
            data_processing_list = []
        for data_processing_ in data_processing_list:
            self.add_data_processing(data_processing_)

    def _create_file_description(self):
        self.writer.file_description(
            self.file_contents_list, self.source_file_list)

    def _create_software_list(self):
        software_list = []
        for sw in self.software_list:
            d = {
                'id': sw.id,
                'version': sw.version
            }
            if sw.is_name(sw.name):
                d[sw.name] = ''
            else:
                d['MS:1000799'] = sw.name
            d['params'] = list(sw.options.items())
            software_list.append(d)

        software_list.append({
            "id": "ms_deisotope_1",
            'MS:1000799': "ms_deisotope"
        })
        self.writer.software_list(software_list)

    def _create_sample_list(self):
        self.writer.sample_list(self.sample_list)

    def build_processing_method(self, order=1, picked_peaks=True, smoothing=True,
                                baseline_reduction=True, additional_parameters=tuple(),
                                software_id=None, data_processing_id=None):
        if software_id is None:
            software_id = "ms_deisotope_1"
        if data_processing_id is None:
            data_processing_id = 'ms_deisotope_processing_%d' % len(
                self.data_processing_list)

        method = data_transformation.ProcessingMethod(software_id=software_id)
        if self.deconvoluted:
            method.add("deisotoping")
            method.add("charge deconvolution")
            method.add("precursor recalculation")

        if picked_peaks:
            method.add("peak picking")
        if smoothing:
            method.add("smoothing")
        if baseline_reduction:
            method.add("baseline reduction")

        method.add("Conversion to mzML")
        method.update(additional_parameters)
        method.update(self.processing_parameters)
        method.order = order
        data_processing_info = data_transformation.DataProcessingInformation(
            [method], data_processing_id)
        # self.add_data_processing(data_processing_info)
        return data_processing_info

    def _create_data_processing_list(self):
        self.writer.data_processing_list(self.data_processing_list)

    def _create_instrument_configuration(self):
        self.writer.instrument_configuration_list(
            self.instrument_configuration_list)

    def _add_spectrum_list(self):
        self._create_file_description()
        self._create_sample_list()
        self._create_software_list()
        self._create_instrument_configuration()
        self._create_data_processing_list()

        self._run_tag = self.writer.run(
            id=self.sample_name or 1,
            sample='sample_1')
        self._run_tag.__enter__()
        self._spectrum_list_tag = self.writer.spectrum_list(
            count=self.n_spectra)
        self._spectrum_list_tag.__enter__()

    def has_started_writing_spectra(self):
        return self._has_started_writing_spectra

    def _pack_activation(self, activation_information):
        """Pack :class:`~.ActivationInformation` into a :class:`dict` structure
        which that :class:`~psims.mzml.writer.MzMLWriter` expects.

        Parameters
        ----------
        activation_information: :class:`~.ActivationInformation`

        Returns
        -------
        :class:`dict`
        """
        params = []
        params.append({
            "name": str(activation_information.method),
        })
        if activation_information.is_multiple_dissociation():
            for method in activation_information.methods[1:]:
                params.append({"name": str(method)})
        # NOTE: Only correct for CID/HCD spectra with absolute collision energies, but that is all I have
        # to test with.
        params.append({
            "name": "collision energy",
            "value": activation_information.energy,
            "unitName": "electron volt"
        })
        if activation_information.is_multiple_dissociation():
            for energy in activation_information.energies[1:]:
                params.append({
                    "name": "collision energy",
                    "value": energy,
                    "unitName": "electron volt"
                })

        for key, val in activation_information.data.items():
            arg = {
                "name": key,
                "value": val
            }
            try:
                arg['unitName'] = val.unit_info
            except AttributeError:
                pass
            params.append(arg)
        return params

    def _pack_precursor_information(self, precursor_information, activation_information=None,
                                    isolation_window=None):
        """Repackage the :class:`~.PrecursorInformation`, :class:`~.ActivationInformation`,
        and :class:~.IsolationWindow` into the nested :class:`dict` structure that
        :class:`~psims.mzml.writer.MzMLWriter` expects.

        Parameters
        ----------
        precursor_information : :class:`~.PrecursorInformation`
        activation_information : :class:`~.ActivationInformation`, optional
        isolation_window : :class:`~.IsolationWindow`, optional

        Returns
        -------
        :class:`dict`
        """
        # If the scan bunch has been fully deconvoluted and it's PrecursorInformation
        # filled in, its extracted fields will be populated and should be used, otherwise
        # use the default read values.
        if precursor_information.extracted_neutral_mass != 0:
            package = {
                "mz": precursor_information.extracted_mz,
                "intensity": precursor_information.extracted_intensity,
                "charge": precursor_information.extracted_charge,
                "scan_id": precursor_information.precursor_scan_id,
                "params": [
                    {"ms_deisotope:defaulted": precursor_information.defaulted},
                    {"ms_deisotope:orphan": precursor_information.orphan}
                ]
            }
        else:
            package = {
                "mz": precursor_information.mz,
                "intensity": precursor_information.intensity,
                "charge": precursor_information.charge,
                "scan_id": precursor_information.precursor_scan_id
            }
        if package['charge'] == ChargeNotProvided:
            package["charge"] = None
        if activation_information is not None:
            package['activation'] = self._pack_activation(
                activation_information)
        if isolation_window is not None:
            package['isolation_window_args'] = {
                "lower": isolation_window.lower,
                "target": isolation_window.target,
                "upper": isolation_window.upper
            }
        return package

    def _prepare_extra_arrays(self, scan):
        extra_arrays = []
        if self.deconvoluted:
            score_array = [
                peak.score for peak in scan.deconvoluted_peak_set
            ]
            extra_arrays.append(("deconvolution score array", score_array))
            envelope_array = envelopes_to_array(
                [peak.envelope for peak in scan.deconvoluted_peak_set])
            extra_arrays.append(("isotopic envelopes array", envelope_array))
        return extra_arrays

    def _get_annotations(self, scan):
        skip = {'filter string', 'base peak intensity', 'base peak m/z', 'lowest observed m/z',
                'highest observed m/z', 'total ion current', }
        annotations = []
        for key, value in scan.annotations.items():
            if key in skip:
                continue
            annotations.append({
                key: value
            })
        return annotations

    def save_scan(self, scan, **kwargs):
        """Write a :class:`~.Scan` to the output document
        as a collection of related :obj:`<spectrum>` tags.

        .. note::

            If no spectra have been written to the output document
            yet, this method will call :meth:`_add_spectrum_list` and
            writes all of the metadata lists out. After this point,
            no new document-level metadata can be added.

        Parameters
        ----------
        scan: :class:`~.Scan`
            The scan to write.
        deconvoluted: :class:`bool`
            Whether the scan to write out should include deconvolution information
        """
        if not self._has_started_writing_spectra:
            self._add_spectrum_list()
            self._has_started_writing_spectra = True

        deconvoluted = kwargs.get("deconvoluted", self.deconvoluted)
        if deconvoluted:
            centroided = True
            precursor_peaks = scan.deconvoluted_peak_set
        elif scan.peak_set:
            centroided = True
            precursor_peaks = scan.peak_set
        else:
            centroided = False
            precursor_peaks = scan.arrays
        polarity = scan.polarity
        if deconvoluted:
            charge_array = [p.charge for p in precursor_peaks]
        else:
            charge_array = None

        if centroided:
            descriptors = SpectrumDescription.from_peak_set(precursor_peaks)
            mz_array = [p.mz for p in precursor_peaks]
            intensity_array = [p.intensity for p in precursor_peaks]
        else:
            descriptors = SpectrumDescription.from_arrays(precursor_peaks)
            mz_array = precursor_peaks.mz
            intensity_array = precursor_peaks.intensity

        instrument_config = scan.instrument_configuration
        if instrument_config is None:
            instrument_config_id = None
        else:
            instrument_config_id = instrument_config.id

        scan_parameters, scan_window_list = self.extract_scan_event_parameters(
            scan)

        if scan.precursor_information:
            precursor_information = self._pack_precursor_information(
                scan.precursor_information,
                scan.activation,
                scan.isolation_window)
        else:
            precursor_information = None

        spectrum_params = [
            {"name": "ms level", "value": scan.ms_level},
            {"name": "MS1 spectrum"} if scan.ms_level == 1 else {"name": "MSn spectrum"},
        ] + list(descriptors)

        spectrum_params.extend(self._get_annotations(scan))

        self.writer.write_spectrum(
            mz_array, intensity_array,
            charge_array,
            id=scan.id, params=spectrum_params,
            centroided=centroided,
            polarity=polarity,
            scan_start_time=scan.scan_time,
            compression=self.compression,
            other_arrays=self._prepare_extra_arrays(scan),
            instrument_configuration_id=instrument_config_id,
            precursor_information=precursor_information,
            scan_params=scan_parameters,
            scan_window_list=scan_window_list,
            encoding=self.data_encoding)

        self.total_ion_chromatogram_tracker[
            scan.scan_time] = (descriptors["total ion current"])
        self.base_peak_chromatogram_tracker[
            scan.scan_time] = (descriptors["base peak intensity"])

    def save_scan_bunch(self, bunch, **kwargs):
        """Write a :class:`~.ScanBunch` to the output document
        as a collection of related :obj:`<spectrum>` tags.

        .. note::

            If no spectra have been written to the output document
            yet, this method will call :meth:`_add_spectrum_list` and
            writes all of the metadata lists out. After this point,
            no new document-level metadata can be added.

        Parameters
        ----------
        bunch : :class:`~.ScanBunch`
            The scan set to write.
        """
        if bunch.precursor is not None:
            self.save_scan(bunch.precursor)

        for prod in bunch.products:
            self.save_scan(prod)

        if self.indexer is not None:
            self.indexer.add_scan_bunch(bunch)

    def extract_scan_event_parameters(self, scan):
        """Package :class:`~.ScanAcquisitionInformation` into a pair of
        :class:`list`s that :class:`~psims.mzml.writer.MzMLWriter` expects.

        Parameters
        ----------
        scan : :class:`~.Scan`

        Returns
        -------
        scan_parameters: :class:`list`
            Parameters qualifying the scan event (:class:`dict`)
        scan_window_list: :class:`list`
            Packed pairs of scan windows (:class:`list`)
        """
        scan_parameters = []
        scan_window_list = []
        acquisition_info = scan.acquisition_information
        filter_line = scan.annotations.get("filter_line")
        if filter_line is not None:
            scan_parameters.append({"name": "filter line", "value": filter_line})
        if acquisition_info is not None and len(acquisition_info) > 0:
            scan_event = acquisition_info[0]
            if scan_event.has_ion_mobility():
                scan_parameters.append({
                    "name": "ion mobility drift time",
                    "value": scan_event.drift_time,
                    "unit_name": "millisecond",
                    'unit_cv_ref': "UO",
                    "unit_accession": 'UO:0000028'
                })
            if scan_event.injection_time is not None:
                scan_parameters.append({
                    "accession": 'MS:1000927', "value": scan_event.injection_time,
                    "unit_name": getattr(scan_event.injection_time, 'unit_info', None),
                })
            traits = scan_event.traits.items()
            for name, value in traits:
                param = {"name": name, "value": value, 'unit_name': getattr(value, 'unit_info', None)}
                scan_parameters.append(param)
            scan_window_list = list(scan_event)
        return scan_parameters, scan_window_list

    def save_chromatogram(self, chromatogram_dict, chromatogram_type, params=None, **kwargs):
        time_array, intensity_array = zip(*chromatogram_dict.items())
        self.writer.write_chromatogram(
            time_array, intensity_array, id=kwargs.get('id'),
            chromatogram_type=chromatogram_type, compression=self.compression,
            params=params)

    def _make_default_chromatograms(self):
        d = dict(
            chromatogram=self.total_ion_chromatogram_tracker,
            chromatogram_type='total ion current chromatogram',
            id='TIC')
        if len(self.total_ion_chromatogram_tracker) > 0:
            self.chromatogram_queue.append(d)

        d = dict(
            chromatogram=self.base_peak_chromatogram_tracker,
            chromatogram_type="basepeak chromatogram",
            id='BPC')
        if len(self.base_peak_chromatogram_tracker) > 0:
            self.chromatogram_queue.append(d)

    def write_chromatograms(self):
        self._chromatogram_list_tag = self.writer.chromatogram_list(
            count=len(self.chromatogram_queue))
        with self._chromatogram_list_tag:
            for chromatogram in self.chromatogram_queue:
                self.save_chromatogram(
                    chromatogram.pop("chromatogram"),
                    **chromatogram)

    def complete(self):
        """Finish writing to the output document.

        This closes the open list tags, empties the chromatogram accumulator,
        and closes the :obj:`<mzML>` tag, and attempts to flush the output file.
        """
        self._spectrum_list_tag.__exit__(None, None, None)
        self._make_default_chromatograms()
        self.write_chromatograms()
        self._run_tag.__exit__(None, None, None)
        self.writer.__exit__(None, None, None)
        if self.indexer is not None:
            try:
                name = self.handle.name
            except AttributeError:
                name = "_detatched_mzml_index"
            try:
                with open(ExtendedScanIndex.index_file_name(name), 'w') as ixfile:
                    self.indexer.serialize(ixfile)
            except IOError as e:
                warnings.warn(
                    "Could not write extended index file due to error %r" % (e,))

        try:
            self.writer.outfile.flush()
        except (IOError, AttributeError, ValueError):
            pass

    def format(self):
        """Formats the output mzML document to be well indented
        and wraps it in an :obj:`<indexedmzML>` structure.

        This method calls :meth:`psims.mzml.writer.MzMLWriter.format`
        which reads the entire document into memory, which may be
        prohibitive for large documents.
        """
        try:
            self.writer.format()
        except OSError as e:
            if on_windows and e.errno == 32:
                pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.complete()
        if hasattr(self.handle, "closed"):
            if not self.handle.closed:
                try:
                    self.handle.close()
                except AttributeError:
                    pass
        self.format()


MzMLScanSerializer = MzMLSerializer


def deserialize_deconvoluted_peak_set(scan_dict):
    envelopes = decode_envelopes(scan_dict["isotopic envelopes array"])
    peaks = []
    mz_array = scan_dict['m/z array']
    intensity_array = scan_dict['intensity array']
    charge_array = scan_dict['charge array']
    score_array = scan_dict['deconvolution score array']
    n = len(scan_dict['m/z array'])
    for i in range(n):
        mz = mz_array[i]
        charge = charge_array[i]
        peak = peak_set.DeconvolutedPeak(
            neutral_mass(mz, charge), intensity_array[i], charge=charge, signal_to_noise=score_array[i],
            index=0, full_width_at_half_max=0, a_to_a2_ratio=0, most_abundant_mass=0,
            average_mass=0, score=score_array[i], envelope=envelopes[i], mz=mz
        )
        peaks.append(peak)
    peaks = peak_set.DeconvolutedPeakSet(peaks)
    peaks._reindex()
    return peaks


def deserialize_peak_set(scan_dict):
    mz_array = scan_dict['m/z array']
    intensity_array = scan_dict['intensity array']
    n = len(scan_dict['m/z array'])
    peaks = []
    for i in range(n):
        peak = FittedPeak(
            mz_array[i], intensity_array[i], 1, i, i,
            0, intensity_array[i], 0, 0)
        peaks.append(peak)
    peak_set = PeakSet(peaks)
    peak_set.reindex()
    return PeakIndex(np.array([]), np.array([]), peak_set)


class ProcessedMzMLDeserializer(MzMLLoader, ScanDeserializerBase):
    """Extends :class:`.MzMLLoader` to support deserializing preprocessed data
    and to take advantage of additional indexing information.

    Attributes
    ----------
    extended_index: :class:`~.ExtendedIndex`
        Holds the additional indexing information
        that may have been generated with the data
        file being accessed.
    sample_run: :class:`SampleRun`

    """

    def __init__(self, source_file, use_index=True):
        MzMLLoader.__init__(self, source_file, use_index=use_index)
        self.extended_index = None
        self._scan_id_to_rt = dict()
        self._sample_run = None
        if self._use_index:
            try:
                if self.has_index_file():
                    self.read_index_file()
                else:
                    self.build_extended_index()
            except IOError:
                pass
            except ValueError:
                pass
            self._build_scan_id_to_rt_cache()

    def read_index_file(self):
        with open(self._index_file_name) as handle:
            self.extended_index = ExtendedScanIndex.deserialize(handle)

    def deserialize_deconvoluted_peak_set(self, scan_dict):
        return deserialize_deconvoluted_peak_set(scan_dict)

    def deserialize_peak_set(self, scan_dict):
        return deserialize_peak_set(scan_dict)

    def has_index_file(self):
        return os.path.exists(self._index_file_name)

    def _make_sample_run(self):
        samples = self.samples()
        sample = samples['sample'][0]
        return SampleRun(name=sample['name'], uuid=sample['SampleRun-UUID'])

    @property
    def sample_run(self):
        if self._sample_run is None:
            self._sample_run = self._make_sample_run()
        return self._sample_run

    def _validate(self, scan):
        return bool(scan.deconvoluted_peak_set) or bool(scan.peak_set)

    def get_index_information_by_scan_id(self, scan_id):
        try:
            try:
                return self.extended_index.msn_ids[scan_id]
            except KeyError:
                return self.extended_index.ms1_ids[scan_id]
        except Exception:
            return {}

    def iter_scan_headers(self, iterator=None):
        self.reset()
        if iterator is None:
            iterator = iter(self._source)
        precursor_scan = None
        product_scans = []

        _make_scan = super(ProcessedMzMLDeserializer, self)._make_scan
        _validate = super(ProcessedMzMLDeserializer, self)._validate

        current_level = 1
        for scan in iterator:
            packed = _make_scan(scan)
            if not _validate(packed):
                continue
            if scan['ms level'] == 2:
                if current_level < 2:
                    current_level = 2
                product_scans.append(packed)
            elif scan['ms level'] == 1:
                if current_level > 1:
                    precursor_scan.product_scans = list(product_scans)
                    yield ScanBunch(precursor_scan, product_scans)
                else:
                    if precursor_scan is not None:
                        precursor_scan.product_scans = list(product_scans)
                        yield ScanBunch(precursor_scan, product_scans)
                precursor_scan = packed
                product_scans = []
            else:
                raise Exception(
                    "This object is not able to handle MS levels higher than 2")
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)
        self.reset()

    def get_scan_header_by_id(self, scan_id):
        """Retrieve the scan object for the specified scan id. If the
        scan object is still bound and in memory somewhere, a reference
        to that same object will be returned. Otherwise, a new object will
        be created.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        Scan
        """
        try:
            packed = super(ProcessedMzMLDeserializer, self)._make_scan(
                self._source.get_by_id(scan_id))
            return packed
        except AttributeError as ae:
            raise AttributeError("Could not read attribute (%s) while looking up scan %s" % (
                ae, scan_id))

    @property
    def _index_file_name(self):
        if isinstance(self.source_file, basestring):
            return ExtendedScanIndex.index_file_name(self.source_file)
        else:
            return ExtendedScanIndex.index_file_name(self.source_file.name)

    def build_extended_index(self, header_only=True):
        self.reset()
        indexer = ExtendedScanIndex()
        iterator = self
        if header_only:
            iterator = self.iter_scan_headers()
        for bunch in iterator:
            indexer.add_scan_bunch(bunch)
        self.reset()
        self.extended_index = indexer
        try:
            with open(self._index_file_name, 'w') as handle:
                indexer.serialize(handle)
        except (IOError, OSError, AttributeError) as err:
            print(err)

    def _make_scan(self, data):
        scan = super(ProcessedMzMLDeserializer, self)._make_scan(data)
        if scan.precursor_information:
            scan.precursor_information.default()
            selected_ion_dict = self._get_selected_ion(data)
            scan.precursor_information.orphan = selected_ion_dict.get(
                "ms_deisotope:orphan") == "true"
            scan.precursor_information.defaulted = selected_ion_dict.get(
                "ms_deisotope:defaulted") == "true"
            scan.annotations['precursor purity'] = data.get(
                'precursor purity', 0)
        if "isotopic envelopes array" in data:
            scan.peak_set = PeakIndex(np.array([]), np.array([]), PeakSet([]))
            scan.deconvoluted_peak_set = self.deserialize_deconvoluted_peak_set(
                data)
            if scan.id in self.extended_index.ms1_ids:
                chosen_indices = self.extended_index.ms1_ids[
                    scan.id]['msms_peaks']
                for ix in chosen_indices:
                    scan.deconvoluted_peak_set[ix].chosen_for_msms = True
        else:
            scan.peak_set = self.deserialize_peak_set(data)
            scan.deconvoluted_peak_set = None
        packed = scan.pack()
        return packed

    def convert_scan_id_to_retention_time(self, scan_id):
        try:
            time = self._scan_id_to_rt[scan_id]
            return time
        except KeyError:
            header = self.get_scan_header_by_id(scan_id)
            return header.scan_time

    def _build_scan_id_to_rt_cache(self):
        if self.extended_index:
            for key in self.extended_index.ms1_ids:
                self._scan_id_to_rt[key] = self.extended_index.ms1_ids[
                    key]['scan_time']
            for key in self.extended_index.msn_ids:
                self._scan_id_to_rt[key] = self.extended_index.msn_ids[
                    key]['scan_time']

    # LC-MS/MS Database API

    def precursor_information(self):
        out = []
        for key, info in self.extended_index.msn_ids.items():
            mz = info['mz']
            neutral_mass = info['neutral_mass']
            charge = info['charge']
            if charge == 'ChargeNotProvided':
                charge = ChargeNotProvided
            intensity = info['intensity']
            precursor_scan_id = info['precursor_scan_id']
            product_scan_id = info['product_scan_id']
            orphan = info.get('orphan', False)
            defaulted = info.get('defaulted', False)
            pinfo = PrecursorInformation(
                mz, intensity, charge, precursor_scan_id,
                self, neutral_mass, charge, intensity,
                product_scan_id=product_scan_id, orphan=orphan,
                defaulted=defaulted)
            out.append(pinfo)
        return out

    def ms1_peaks_above(self, mass_threshold=500, intensity_threshold=1000.):
        accumulate = []
        for ms1_id in self.extended_index.ms1_ids:
            scan = self.get_scan_by_id(ms1_id)
            for peak in scan.deconvoluted_peak_set:
                if peak.intensity < intensity_threshold or peak.neutral_mass < mass_threshold:
                    continue
                accumulate.append((ms1_id, peak, id(peak)))
        return accumulate

    def ms1_scan_times(self):
        times = sorted(
            [bundle['scan_time'] for bundle in
             self.extended_index.ms1_ids.values()])
        return np.array(times)

    def extract_total_ion_current_chromatogram(self):
        current = []
        for scan_id in self.extended_index.ms1_ids:
            header = self.get_scan_header_by_id(scan_id)
            current.append(header.arrays[1].sum())
        return np.array(current)

    def msms_for(self, neutral_mass, mass_error_tolerance=1e-5, start_time=None, end_time=None):
        m = neutral_mass
        w = neutral_mass * mass_error_tolerance
        lo = m - w
        hi = m + w
        out = []
        for pinfo in self.precursor_information():
            if lo <= pinfo.neutral_mass <= hi:
                valid = True
                product = None
                if start_time is not None or end_time is not None:
                    product = self.get_scan_header_by_id(pinfo.product_scan_id)
                if start_time is not None and product.scan_time < start_time:
                    valid = False
                elif end_time is not None and product.scan_time > end_time:
                    valid = False
                if valid:
                    out.append(pinfo)
        return out


try:
    has_c = True
    _deserialize_deconvoluted_peak_set = deserialize_deconvoluted_peak_set
    from ms_deisotope._c.utils import deserialize_deconvoluted_peak_set
except ImportError:
    has_c = False
