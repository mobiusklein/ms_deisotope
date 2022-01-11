'''
Writing mzML
------------

Using the :mod:`psims` library, :mod:`ms_deisotope.output.mzml` can write an mzML
file with all associated metadata, including deconvoluted peak arrays, chromatograms,
and data transformations. The :class:`~.MzMLSerializer` class handles all facets of
this process.

This module also contains a specialized version of :class:`~.MzMLLoader`,
:class:`~.ProcessedMzMLLoader`, which can directly reconstruct each
deconvoluted peak list and provides fast access to an extended index of
metadata that :class:`~.MzMLSerializer` writes to an external file.


.. code:: python

    import ms_deisotope
    from ms_deisotope.test.common import datafile
    from ms_deisotope.output.mzml import MzMLSerializer

    reader = ms_deisotope.MSFileLoader(datafile("small.mzML"))
    with open("small.deconvoluted.mzML", 'wb') as fh:
        writer = MzMLSerializer(fh, n_spectra=len(reader))

        writer.copy_metadata_from(reader)
        for bunch in reader:
            bunch.precursor.pick_peaks()
            bunch.precursor.deconvolute()
            for product in bunch.products:
                product.pick_peaks()
                product.deconvolute()
            writer.save(bunch)

        writer.close()

'''
import hashlib
import array
import warnings


from contextlib import contextmanager
from collections import OrderedDict
try:
    from collections.abc import Sequence, Mapping
except ImportError:
    from collections import Sequence, Mapping
from uuid import uuid4, UUID

import numpy as np

from ms_peak_picker import PeakIndex, PeakSet, FittedPeak

try:
    from psims.mzml import writer
except ImportError:
    writer = None


from ms_deisotope import version as lib_version
from ms_deisotope.peak_set import (DeconvolutedPeak, DeconvolutedPeakSet, Envelope, IonMobilityDeconvolutedPeak)
from ms_deisotope.averagine import neutral_mass
from ms_deisotope.qc.isolation import CoIsolation

from ms_deisotope.data_source.common import (
    ChargeNotProvided,
    _SingleScanIteratorImpl,
    _InterleavedGroupedScanIteratorImpl)
from ms_deisotope.data_source.metadata import data_transformation
from ms_deisotope.data_source.metadata.software import (Software, software_name)
from ms_deisotope.data_source.mzml import MzMLLoader
from ms_deisotope.data_source.scan.mobility_frame import Generic3DIonMobilityFrameSource

from ms_deisotope.feature_map import ExtendedScanIndex
from ms_deisotope.feature_map.feature_fit import DeconvolutedLCMSFeature, DeconvolutedLCMSFeatureTreeNode
from ms_deisotope.feature_map.feature_map import DeconvolutedLCMSFeatureMap, LCMSFeature, LCMSFeatureMap


from .common import ScanSerializerBase, ScanDeserializerBase, SampleRun, LCMSMSQueryInterfaceMixin
from .text_utils import (envelopes_to_array, decode_envelopes)


class SpectrumDescription(Sequence):
    '''A helper class to calculate properties of a spectrum derived from
    their peak data or raw signal
    '''
    def __init__(self, attribs=None):
        Sequence.__init__(self)
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
        """Add the descriptor `desc` to the collection

        Adds a descriptor which conforms to any of :mod:`psims`'s
        cvParam specification patterns.

        Parameters
        ----------
        desc : :class:`dict`, :class:`tuple`, or :class:`str`
            The descriptor to add.

        Returns
        -------
        int:
            The current size of the description list
        """
        i = len(self)
        self.descriptors.append(desc)
        return i

    def __repr__(self):
        return "{self.__class__.__name__}({self.descriptors})".format(self=self)

    @classmethod
    def from_peak_set(cls, peak_list):
        """Calculate the spectrum's descriptors from a :class:`Sequence` of :class:`~.PeakLike`
        objects.

        Parameters
        ----------
        peak_list : :class:`Sequence`
            The peaks to calculate properties from.

        Returns
        -------
        :class:`SpectrumDescription`
        """
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
        """Calculate the spectrum's descriptors from a :class:`RawDataArrays`
        instance.

        Parameters
        ----------
        arrays : :class:`RawDataArrays`
            The signal to calculate properties from.

        Returns
        -------
        :class:`SpectrumDescription`
        """
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

    try:
        default_compression = writer.COMPRESSION_ZLIB
    except AttributeError:
        default_compression = None

    file_extensions = {
        "mzml",
        "mzml.gz",
    }

    default_data_encoding = {
        writer.MZ_ARRAY: np.float64,
        writer.INTENSITY_ARRAY: np.float32,
        writer.CHARGE_ARRAY: np.int32,
    }

    _format_conversion_term = "Conversion to mzML"

    def __init__(self, handle, n_spectra=int(2e5), compression=None,
                 deconvoluted=True, sample_name=None, build_extra_index=True,
                 data_encoding=None, include_software_entry=True):
        if data_encoding is None:
            data_encoding = self.default_data_encoding
        if writer is None:
            raise ImportError(
                "Cannot write mzML without psims. Please install psims to use this feature.")
        if compression is None:
            compression = self.default_compression
        super(MzMLSerializer, self).__init__()
        self.handle = handle
        self.n_spectra = n_spectra
        self.compression = compression
        self.data_encoding = data_encoding
        self._has_started_writing_spectra = False

        self.writer = self._make_writer(handle)
        self.writer.begin()
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
        self._include_software_entry = include_software_entry
        self._this_software = None

    def _make_writer(self, handle):
        return writer.MzMLWriter(handle)

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
        '''Add the information of a :class:`~.FileInformation` to the
        output document.

        Parameters
        ----------
        file_information: :class:`~.FileInformation`
            The information to add.
        '''
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
        """Remove a key to the resulting :obj:`<fileDescription>`
        of the output document.

        Parameters
        ----------
        file_contents: :class:`str` or :class:`Mapping`
            The parameter to remove

        Raises
        ------
        KeyError:
            When the content is not found.
        """
        i = None
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
        if i is None:
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

        for key, value in source_file.parameters.items():
            accession = getattr(key, 'accession', str(key))
            if accession is None:
                accession = str(key)
            unwrapped['params'].append((accession, value))

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
        """Copies the file-level metadata from an instance of :class:`~.ScanFileMetadataBase`
        into the metadata of the file to be written

        Parameters
        ----------
        reader : :class:`~.ScanFileMetadataBase`
            The source to copy metadata from
        """
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
        if self._include_software_entry:
            if self._this_software is None:
                self._this_software = this_software = self._make_software_entry()
            else:
                this_software = self._this_software
            self.software_list.append(this_software)
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

        self.writer.software_list(software_list)

    def _make_software_entry(self):
        ms_deisotope_entries = []
        for sw in self.software_list:
            if 'ms_deisotope' in str(sw.id):
                ms_deisotope_entries.append(str(sw.id))
        for i in range(1, 100):
            query = 'ms_deisotope_%d' % i
            if query in ms_deisotope_entries:
                continue
            else:
                new_entry_id = query
                break
        else:
            new_entry_id = 'ms_deisotope_%s' % str(uuid4())
        self._this_software = inst = Software("ms_deisotope", new_entry_id, lib_version)
        return inst

    def _create_sample_list(self):
        self.writer.sample_list(self.sample_list)

    def build_processing_method(self, order=1, picked_peaks=True, smoothing=True,
                                baseline_reduction=True, additional_parameters=tuple(),
                                software_id=None, data_processing_id=None):
        if software_id is None:
            if self._this_software is None:
                self._make_software_entry()
            software_id = self._this_software.id
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

        method.add(self._format_conversion_term)
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
            "unit_name": "electronvolt",
            'unit_accession': 'UO:0000266'
        })
        if activation_information.is_multiple_dissociation():
            energies = activation_information.energies[1:]
            supplemental_energy = None
            if activation_information.has_supplemental_dissociation():
                supplemental_energy = energies[-1]
                energies = energies[:-1]
            for energy in energies:
                params.append({
                    "name": "collision energy",
                    "value": energy,
                    "unit_name": "electronvolt",
                    'unit_accession': 'UO:0000266'
                })
            if supplemental_energy is not None:
                params.append({
                    "name": 'supplemental collision energy',
                    "value": supplemental_energy,
                    "unit_name": "electronvolt",
                    'unit_accession': 'UO:0000266'
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
        package = {}
        # If the scan bunch has been fully deconvoluted and it's PrecursorInformation
        # filled in, its extracted fields will be populated and should be used, otherwise
        # use the default read values.
        if precursor_information is not None:
            extracted_neutral_mass = precursor_information.extracted_neutral_mass
            if (extracted_neutral_mass != 0):
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
                    "scan_id": precursor_information.precursor_scan_id,
                    "params": []
                }
            # This implicitly captures ion mobility which is stored as an annotation key-value pair.
            for key, value in precursor_information.annotations.items():
                package['params'].append({
                    key: value
                })
            if precursor_information.coisolation:
                for p in precursor_information.coisolation:
                    package['params'].append({
                        "name": "ms_deisotope:coisolation",
                        "value": "%f %f %d" % (p.neutral_mass, p.intensity, p.charge)
                    })
        else:
            package['mz'] = None
            package["charge"] = None
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

    def _prepare_extra_arrays(self, scan, **kwargs):
        deconvoluted = kwargs.get("deconvoluted", self.deconvoluted)
        extra_arrays = []
        if deconvoluted:
            score_array = [
                peak.score for peak in scan.deconvoluted_peak_set
            ]
            extra_arrays.append(("deconvolution score array", score_array))
            envelope_array = envelopes_to_array(
                [peak.envelope for peak in scan.deconvoluted_peak_set])
            extra_arrays.append(("isotopic envelopes array", envelope_array))
        else:
            if scan.arrays.data_arrays:
                extra_arrays.extend(sorted(scan.arrays.data_arrays.items()))
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

    def _get_peak_data(self, scan, kwargs):
        deconvoluted = kwargs.get("deconvoluted", self.deconvoluted)
        if deconvoluted:
            centroided = True
            peak_data = scan.deconvoluted_peak_set
        elif scan.peak_set is not None:
            centroided = True
            peak_data = scan.peak_set
        else:
            centroided = not scan.is_profile
            peak_data = scan.arrays
        if deconvoluted:
            charge_array = [p.charge for p in peak_data]
        else:
            charge_array = None

        if centroided:
            descriptors = SpectrumDescription.from_peak_set(peak_data)
            mz_array = [p.mz for p in peak_data]
            intensity_array = [p.intensity for p in peak_data]
        else:
            descriptors = SpectrumDescription.from_arrays(peak_data)
            mz_array = peak_data.mz
            intensity_array = peak_data.intensity

        other_arrays = self._prepare_extra_arrays(scan, deconvoluted=deconvoluted)

        return (centroided, descriptors, mz_array, intensity_array,
                charge_array, other_arrays)

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

        (centroided, descriptors, mz_array, intensity_array,
         charge_array, other_arrays) = self._get_peak_data(scan, kwargs)
        polarity = scan.polarity
        try:
            instrument_config = scan.instrument_configuration
        except AttributeError:
            instrument_config = None
        if instrument_config is None:
            instrument_config_id = None
        else:
            instrument_config_id = instrument_config.id

        scan_parameters, scan_window_list = self.extract_scan_event_parameters(scan)

        if (scan.precursor_information or scan.isolation_window or scan.activation):
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

        try:
            spectrum_params.extend(self._get_annotations(scan))
        except AttributeError:
            pass

        self.writer.write_spectrum(
            mz_array, intensity_array,
            charge_array,
            id=scan.id, params=spectrum_params,
            centroided=centroided,
            polarity=polarity,
            scan_start_time=scan.scan_time,
            compression=self.compression,
            other_arrays=other_arrays,
            instrument_configuration_id=instrument_config_id,
            precursor_information=precursor_information,
            scan_params=scan_parameters,
            scan_window_list=scan_window_list,
            encoding=self.data_encoding)

        self.total_ion_chromatogram_tracker[
            scan.scan_time] = (descriptors["total ion current"])
        self.base_peak_chromatogram_tracker[
            scan.scan_time] = (descriptors["base peak intensity"])

        if kwargs.get("_include_in_index", True) and self.indexer is not None:
            self.indexer.add_scan(scan)

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
            self.save_scan(bunch.precursor, _include_in_index=False, **kwargs)

        for prod in bunch.products:
            self.save_scan(prod, _include_in_index=False, **kwargs)

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
        filter_string = scan.annotations.get("filter_string")
        if filter_string is not None:
            scan_parameters.append({"name": "filter string", "value": filter_string})
        if acquisition_info is not None and len(acquisition_info) > 0:
            scan_event = acquisition_info[0]
            if scan_event.injection_time is not None:
                scan_parameters.append({
                    "accession": 'MS:1000927', "value": scan_event.injection_time,
                    "unit_name": getattr(scan_event.injection_time, 'unit_info', None),
                })
            traits = scan_event.traits.items()
            for name, value in traits:
                param = {"name": str(name), "value": value, 'unit_name': getattr(value, 'unit_info', None)}
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
        if self._spectrum_list_tag is not None:
            self._spectrum_list_tag.__exit__(None, None, None)
        if self._run_tag is not None:
            self._make_default_chromatograms()
            self.write_chromatograms()
        if self._run_tag is not None:
            self._run_tag.__exit__(None, None, None)
        self.writer.__exit__(None, None, None)
        if self.indexer is not None:
            try:
                name = self.handle.name
                try:
                    with open(ExtendedScanIndex.index_file_name(name), 'w') as ixfile:
                        self.indexer.serialize(ixfile)
                except IOError as e:
                    warnings.warn(
                        "Could not write extended index file due to error %r" % (e,))
            except AttributeError:
                warnings.warn("Could not determine name to write extended index file to")

        try:
            self.writer.outfile.flush()
        except (IOError, AttributeError, ValueError):
            pass

    def format(self):
        """This method is no longer needed.
        """
        pass

    def close(self):
        self.complete()
        if hasattr(self.handle, "closed"):
            if not self.handle.closed:
                try:
                    self.handle.close()
                except AttributeError:
                    pass
        else:
            try:
                self.handle.close()
            except (AttributeError, ValueError, TypeError, OSError):
                pass


MzMLScanSerializer = MzMLSerializer


def deserialize_deconvoluted_peak_set(scan_dict, include_envelopes=True):
    if include_envelopes:
        envelopes = decode_envelopes(scan_dict["isotopic envelopes array"])
    else:
        envelopes = None
    peaks = []
    mz_array = scan_dict['m/z array']
    intensity_array = scan_dict['intensity array']
    charge_array = scan_dict['charge array']
    score_array = scan_dict['deconvolution score array']
    n = len(scan_dict['m/z array'])
    for i in range(n):
        mz = mz_array[i]
        charge = charge_array[i]
        peak = DeconvolutedPeak(
            neutral_mass(mz, charge), intensity_array[i], charge=charge, signal_to_noise=score_array[i],
            index=0, full_width_at_half_max=0, a_to_a2_ratio=0, most_abundant_mass=0,
            average_mass=0, score=score_array[i], envelope=envelopes[i] if include_envelopes else None, mz=mz
        )
        peaks.append(peak)
    peaks = DeconvolutedPeakSet(peaks)
    peaks.reindex()
    return peaks


def deserialize_deconvoluted_ion_mobility_peak_set(scan_dict):
    envelopes = decode_envelopes(scan_dict["isotopic envelopes array"])
    peaks = []
    mz_array = scan_dict['m/z array']
    intensity_array = scan_dict['intensity array']
    charge_array = scan_dict['charge array']
    score_array = scan_dict['deconvolution score array']
    drift_time_array = scan_dict['mean drift time array']
    n = len(scan_dict['m/z array'])
    for i in range(n):
        mz = mz_array[i]
        charge = charge_array[i]
        peak = IonMobilityDeconvolutedPeak(
            neutral_mass(mz, charge), intensity_array[i], charge=charge, signal_to_noise=score_array[i],
            index=0, full_width_at_half_max=0, a_to_a2_ratio=0, most_abundant_mass=0,
            average_mass=0, score=score_array[i], envelope=envelopes[i], mz=mz, drift_time=drift_time_array[i],
        )
        peaks.append(peak)
    peaks = DeconvolutedPeakSet(peaks)
    peaks.reindex()
    return peaks


def deserialize_external_deconvoluted_peaks(scan_dict, fill_envelopes=True, averagine=None):
    if averagine is None:
        from ms_deisotope import peptide
        averagine = peptide
    peaks = []
    mz_array = scan_dict['m/z array']
    intensity_array = scan_dict['intensity array']
    charge_array = scan_dict['charge array']
    for i in range(len(mz_array)):
        mz = mz_array[i]
        charge = charge_array[i]
        mass = neutral_mass(mz, charge)
        intensity = intensity_array[i]
        if fill_envelopes:
            envelope = averagine.isotopic_cluster(mz, charge, truncate_after=0.8)
            envelope = Envelope([(p.mz, p.intensity) for p in envelope.scale_raw(intensity)])
        else:
            envelope = Envelope([])
        peak = DeconvolutedPeak(
            mass, intensity, charge=charge, signal_to_noise=intensity,
            index=0, full_width_at_half_max=0, a_to_a2_ratio=0, most_abundant_mass=0,
            average_mass=0, score=intensity, envelope=envelope, mz=mz
        )
        peaks.append(peak)
    peaks = DeconvolutedPeakSet(peaks)
    peaks.reindex()
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


class PeakSetDeserializingMixin(object):
    parse_peaks = True
    parse_envelopes = True

    def deserialize_deconvoluted_peak_set(self, scan_dict):
        try:
            if "mean drift time array" in scan_dict:
                return deserialize_deconvoluted_ion_mobility_peak_set(scan_dict)
            return deserialize_deconvoluted_peak_set(scan_dict)
        except KeyError as err:
            if "charge array" in scan_dict and "isotopic envelopes array" not in scan_dict:
                return self.deserialize_external_deconvoluted_peak_set(scan_dict)
            else:
                raise err

    def deserialize_external_deconvoluted_peak_set(self, scan_dict):
        return deserialize_external_deconvoluted_peaks(scan_dict)

    def deserialize_peak_set(self, scan_dict):
        return deserialize_peak_set(scan_dict)

    def _validate(self, scan):
        if self.parse_peaks:
            return bool(scan.deconvoluted_peak_set) or bool(scan.peak_set)
        else:
            return True

    @contextmanager
    def toggle_peak_loading(self):
        self.parse_peaks = False
        self.decode_binary = False
        yield self
        self.reset()
        self.parse_peaks = True
        self.decode_binary = True

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
        precursor = super(PeakSetDeserializingMixin,
                          self)._precursor_information(scan)
        if precursor is None:
            return None
        precursor.orphan = precursor.annotations.pop(
            "ms_deisotope:orphan", None) == "true"
        precursor.defaulted = precursor.annotations.pop(
            "ms_deisotope:defaulted", None) == "true"
        coisolation_params = precursor.annotations.pop(
            "ms_deisotope:coisolation", [])
        if not isinstance(coisolation_params, list):
            coisolation_params = [coisolation_params]
        coisolation = []
        for entry in coisolation_params:
            mass, intensity, charge = entry.split(" ")
            coisolation.append(
                CoIsolation(float(mass), float(intensity), int(charge)))
        precursor.coisolation = coisolation
        return precursor

    def _make_scan(self, data):
        scan = super(PeakSetDeserializingMixin, self)._make_scan(data)
        try:
            precursor_information = scan.precursor_information
            if precursor_information:
                scan.precursor_information.default()
                selected_ion_dict = self._get_selected_ion(data)
                scan.precursor_information.orphan = selected_ion_dict.get(
                    "ms_deisotope:orphan") == "true"
                scan.precursor_information.defaulted = selected_ion_dict.get(
                    "ms_deisotope:defaulted") == "true"
                scan.annotations['precursor purity'] = data.get(
                    'precursor purity', 0)
        except KeyError:
            pass
        if self.parse_peaks:
            if "m/z array" not in data:
                warnings.warn("No m/z array found for scan %r" % (scan.id, ))
                scan.peak_set = PeakIndex(np.array([]), np.array([]), PeakSet([]))
                scan.deconvoluted_peak_set = DeconvolutedPeakSet([])
            elif "charge array" in data:
                scan.peak_set = PeakIndex(np.array([]), np.array([]), PeakSet([]))
                scan.deconvoluted_peak_set = self.deserialize_deconvoluted_peak_set(
                    data)
                if self.has_extended_index() and scan.id in self.extended_index.ms1_ids:
                    chosen_indices = self.extended_index.ms1_ids[
                        scan.id]['msms_peaks']
                    for ix in chosen_indices:
                        scan.deconvoluted_peak_set[ix].chosen_for_msms = True
            else:
                scan.peak_set = self.deserialize_peak_set(data)
                scan.deconvoluted_peak_set = None
        else:
            scan.peak_set = None
            scan.deconvoluted_peak_set = None
        packed = scan.pack(bind=True)
        return packed

    # These methods have to be replicated over from LCMSQueryInterface because of the extra
    # round of overloading `_make_scan`

    def iter_scan_headers(self, iterator=None, grouped=True):
        try:
            if not self._has_ms1_scans():
                grouped = False
        except Exception:
            pass
        self.reset()
        if iterator is None:
            iterator = iter(self._source)

        _make_scan = super(PeakSetDeserializingMixin, self)._make_scan
        _validate = super(PeakSetDeserializingMixin, self)._validate

        if grouped:
            impl = _InterleavedGroupedScanIteratorImpl(
                iterator, _make_scan, _validate)
        else:
            impl = _SingleScanIteratorImpl(iterator, _make_scan, _validate)

        for x in impl:
            yield x

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
            packed = super(PeakSetDeserializingMixin, self)._make_scan(
                self._source.get_by_id(scan_id))
            return packed
        except AttributeError as ae:
            raise AttributeError("Could not read attribute (%s) while looking up scan %s" % (
                ae, scan_id))


class ProcessedMzMLLoader(PeakSetDeserializingMixin, MzMLLoader, ScanDeserializerBase, LCMSMSQueryInterfaceMixin):
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

    file_extensions = {
        "mzml",
        "mzml.gz",
    }

    def __init__(self, source_file, use_index=True, use_extended_index=True):
        super(ProcessedMzMLLoader, self).__init__(source_file, use_index=use_index, decode_binary=True)
        self.extended_index = None
        self._scan_id_to_rt = dict()
        self._sample_run = None
        self._use_extended_index = use_extended_index
        if self._use_index:
            if self._use_extended_index:
                self.require_extended_index()

    def _dispose(self):
        self._scan_id_to_rt.clear()
        self.extended_index.clear()
        super(ProcessedMzMLLoader, self)._dispose()

    def __reduce__(self):
        return self.__class__, (self.source_file, self._use_index, self._use_extended_index)

    def _make_sample_run(self):
        samples = self.samples()
        sample = samples[0]
        uuid_from_sample = None
        try:
            uuid_from_sample = sample['SampleRun-UUID']
        except KeyError:
            try:
                source_name = self.source_file_name
            except AttributeError:
                source_name = str(self)
            uuid_from_sample = str(UUID(hashlib.new('md5', source_name).hexdigest()))
        return SampleRun(name=sample.name, uuid=uuid_from_sample, **dict(sample.items()))

    @property
    def sample_run(self):
        if self._sample_run is None:
            self._sample_run = self._make_sample_run()
        return self._sample_run


ProcessedMzMLDeserializer = ProcessedMzMLLoader

try:
    has_c = True
    _deserialize_deconvoluted_peak_set = deserialize_deconvoluted_peak_set
    _deserialize_peak_set = deserialize_peak_set
    from ms_deisotope._c.utils import deserialize_deconvoluted_peak_set, deserialize_peak_set
except ImportError:
    has_c = False


def extracted_features_to_3d_arrays(features):
    mz_array = array.array('d')
    intensity_array = array.array('d')
    ion_mobility_array = array.array('d')
    feature_id_array = array.array('L')
    for i, feature in enumerate(features):
        for node in feature:
            time = node.time
            for peak in node.members:
                ion_mobility_array.append(time)
                mz_array.append(peak.mz)
                intensity_array.append(peak.intensity)
                feature_id_array.append(i)
    mz_array = np.array(mz_array, copy=False)
    intensity_array = np.array(intensity_array, copy=False)
    ion_mobility_array = np.array(ion_mobility_array, copy=False)
    feature_id_array = np.array(feature_id_array, copy=False)
    mask = np.lexsort(np.stack((ion_mobility_array, mz_array)))
    return (mz_array[mask], intensity_array[mask], ion_mobility_array[mask], feature_id_array[mask])


def deserialize_features(scan_dict, ion_mobility_array_name='raw ion mobility array'):
    mz_array = scan_dict['m/z array']
    intensity_array = scan_dict['intensity array']
    drift_time_array = scan_dict[ion_mobility_array_name]
    feature_id_array = scan_dict['feature id array']
    n = len(mz_array)
    features = dict()
    for i in range(n):
        mz = mz_array[i]
        im = drift_time_array[i]
        inten = intensity_array[i]
        peak = FittedPeak(
            mz, inten, signal_to_noise=inten,
            index=0, full_width_at_half_max=0)
        try:
            feature = features[feature_id_array[i]]
        except KeyError:
            feature = features[feature_id_array[i]] = LCMSFeature([])
        feature.insert(peak, im)
    feature_map = LCMSFeatureMap(list(features.values()))
    return feature_map


def deconvoluted_features_to_3d_arrays(features):
    mz_array = array.array('d')
    intensity_array = array.array('d')
    charge_array = array.array('i')
    score_array = array.array('d')
    ion_mobility_array = array.array('d')
    envelopes = []
    feature_id_array = array.array('L')
    point_count = 0
    for i, feature in enumerate(features):
        for node in feature:
            time = node.time
            for peak in node.members:
                ion_mobility_array.append(time)
                mz_array.append(peak.mz)
                intensity_array.append(peak.intensity)
                score_array.append(peak.score)
                charge_array.append(peak.charge)
                envelopes.append(peak.envelope)
                point_count += (len(peak.envelope) + 1) * 2
                feature_id_array.append(i)
    mz_array = np.array(mz_array, copy=False)
    intensity_array = np.array(intensity_array, copy=False)
    charge_array = np.array(charge_array, copy=False)
    score_array = np.array(score_array, copy=False)
    ion_mobility_array = np.array(ion_mobility_array, copy=False)
    feature_id_array = np.array(feature_id_array, copy=False)
    mask = np.lexsort(np.stack((ion_mobility_array, mz_array)))
    envelope_array = np.zeros(point_count, dtype=np.float32)
    k = 0
    for j in mask:
        for point in envelopes[j]:
            envelope_array[k] = point.mz
            envelope_array[k+1] = point.intensity
            k += 2
        k += 2
    return (mz_array[mask], intensity_array[mask], charge_array[mask],
            score_array[mask], ion_mobility_array[mask], envelope_array,
            feature_id_array[mask])


def deserialize_deconvoluted_features(scan_dict, ion_mobility_array_name='raw ion mobility array'):
    envelopes = decode_envelopes(scan_dict["isotopic envelopes array"])
    mz_array = scan_dict['m/z array']
    intensity_array = scan_dict['intensity array']
    drift_time_array = scan_dict[ion_mobility_array_name]
    feature_id_array = scan_dict['feature id array']
    charge_array = scan_dict['charge array']
    score_array = scan_dict['deconvolution score array']
    n = len(mz_array)
    features = dict()
    for i in range(n):
        mz = mz_array[i]
        charge = charge_array[i]
        im = drift_time_array[i]
        peak = DeconvolutedPeak(
            neutral_mass(mz, charge), intensity_array[i], charge=charge, signal_to_noise=score_array[i],
            index=0, full_width_at_half_max=0, a_to_a2_ratio=0, most_abundant_mass=0,
            average_mass=0, score=score_array[i], envelope=envelopes[i], mz=mz
        )
        try:
            feature = features[feature_id_array[i]]
        except KeyError:
            feature = features[feature_id_array[i]
                               ] = DeconvolutedLCMSFeature([], charge)
        feature.insert(peak, im)
    feature_map = DeconvolutedLCMSFeatureMap(list(features.values()))
    return feature_map


class IonMobilityAware3DMzMLSerializer(MzMLSerializer):
    default_data_encoding = MzMLSerializer.default_data_encoding.copy()
    default_data_encoding.update({
        "feature id array": np.int32,
    })
    def _get_peak_data(self, scan, kwargs):
        deconvoluted = kwargs.get("deconvoluted", self.deconvoluted)
        if deconvoluted:
            centroided = True
            peak_data = scan.deconvoluted_features
        elif scan.features:
            centroided = True
            peak_data = scan.features
        else:
            centroided = False
            peak_data = scan.arrays

        if centroided:
            descriptors = SpectrumDescription.from_peak_set(peak_data)
        else:
            descriptors = SpectrumDescription.from_arrays(peak_data)

        if deconvoluted:
            (mz_array, intensity_array, charge_array, score_array,
             ion_mobility_array, envelope_array, feature_id_array) = deconvoluted_features_to_3d_arrays(peak_data)
            other_arrays = [
                ('raw ion mobility array', ion_mobility_array),
                ('isotopic envelopes array', envelope_array),
                ('deconvolution score array', score_array),
                ('feature id array', feature_id_array),
            ]
        elif centroided:
            (mz_array, intensity_array, ion_mobility_array,
             feature_id_array) = extracted_features_to_3d_arrays(peak_data)
            charge_array = None
            other_arrays = [
                ('raw ion mobility array', ion_mobility_array),
                ('feature id array', feature_id_array),
            ]
        else:
            mz_array = peak_data.mz
            intensity_array = peak_data.intensity
            charge_array = None
            ion_mobility_array = peak_data.ion_mobility
            other_arrays = [
                ('raw ion mobility array', ion_mobility_array),
            ]
        return (centroided, descriptors, mz_array, intensity_array,
                charge_array, other_arrays)


class ProcessedGeneric3DIonMobilityFrameSource(Generic3DIonMobilityFrameSource):
    def _make_frame(self, data):
        frame = super(ProcessedGeneric3DIonMobilityFrameSource, self)._make_frame(data)
        if 'feature id array' in frame._data:
            if 'charge array' in frame._data:
                frame.deconvoluted_features = deserialize_deconvoluted_features(frame._data)
            else:
                frame.features = deserialize_features(frame._data)
        return frame
