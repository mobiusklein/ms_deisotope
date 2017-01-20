import os
import json
from collections import OrderedDict

import numpy as np

from ms_peak_picker import PeakIndex, PeakSet
from psims.mzml import writer

from .common import ScanSerializerBase, ScanDeserializerBase
from .text_utils import (envelopes_to_array, decode_envelopes)
from ms_deisotope import peak_set
from ms_deisotope.averagine import neutral_mass
from ms_deisotope.data_source.common import PrecursorInformation, ScanBunch
from ms_deisotope.data_source.mzml import MzMLLoader


def _base_peak_from_descriptors(descriptors):
    return descriptors[1]['value']


def _total_intensity_from_descriptors(descriptors):
    return descriptors[2]['value']


def describe_spectrum(peak_list):
    descriptors = []
    base_peak = max(peak_list, key=lambda x: x.intensity)
    descriptors.append({
        "name": "base peak m/z",
        "value": base_peak.mz,
    })
    descriptors.append({
        "name": "base peak intensity",
        "value": base_peak.intensity
    })
    descriptors.append({
        "name": "total ion current",
        "value": sum(p.intensity for p in peak_list)
    })
    peaks_mz_order = sorted(peak_list, key=lambda x: x.mz)
    descriptors.append({
        "name": "lowest observed m/z",
        "value": peaks_mz_order[0].mz
    })
    descriptors.append({
        "name": "highest observed m/z",
        "value": peaks_mz_order[-1].mz
    })
    return descriptors


class ExtendedScanIndex(object):
    def __init__(self, ms1_ids=None, msn_ids=None):
        if ms1_ids is None:
            ms1_ids = {}
        if msn_ids is None:
            msn_ids = {}
        self.ms1_ids = OrderedDict(ms1_ids)
        self.msn_ids = OrderedDict(msn_ids)

    def _package_precursor_information(self, product):
        precursor_information = product.precursor_information
        if precursor_information.extracted_neutral_mass != 0:
            package = {
                "neutral_mass": precursor_information.extracted_neutral_mass,
                "mz": precursor_information.extracted_mz,
                "intensity": precursor_information.extracted_intensity,
                "charge": precursor_information.extracted_charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "product_scan_id": product.id
            }
        else:
            package = {
                "neutral_mass": neutral_mass(
                    precursor_information.mz, precursor_information.charge),
                "mz": precursor_information.mz,
                "intensity": precursor_information.intensity,
                "charge": precursor_information.charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "product_scan_id": product.id
            }
        return package

    def add_scan_bunch(self, bunch):
        self.ms1_ids[bunch.precursor.id] = bunch.precursor.scan_time
        for product in bunch.products:
            self.msn_ids[product.id] = self._package_precursor_information(product)

    def serialize(self, handle):
        mapping = {
            "ms1_ids": list(self.ms1_ids.items()),
            "msn_ids": list(self.msn_ids.items())
        }
        json.dump(mapping, handle)

    @staticmethod
    def index_file_name(name):
        return name + '-idx.json'

    @classmethod
    def deserialize(cls, handle):
        mapping = json.load(handle)
        return cls(**mapping)


class MzMLScanSerializer(ScanSerializerBase):

    def __init__(self, handle, n_spectra=2e4, compression="zlib", deconvoluted=True,
                 sample_name=None, build_extra_index=True):
        self.handle = handle
        self.writer = writer.MzMLWriter(handle)
        self.n_spectra = n_spectra
        self.compression = compression
        self._has_started_writing_spectra = False
        self.writer.__enter__()
        self.context_stack = [self.writer]
        self.writer.controlled_vocabularies()
        self.deconvoluted = deconvoluted
        self.sample_name = sample_name

        self.file_contents_list = []
        self.software_list = []
        self.source_file_list = []
        self.data_processing_list = []
        self.instrument_configuration_list = []

        self.total_ion_chromatogram_tracker = OrderedDict()
        self.base_peak_chromatogram_tracker = OrderedDict()
        self.indexer = None
        if build_extra_index:
            self.indexer = ExtendedScanIndex()

    def add_software(self, software_description):
        self.software_list.append(software_description)

    def add_file_contents(self, file_contents):
        self.file_contents_list.append(file_contents)

    def add_source_file(self, source_file_description):
        self.source_file_list.append(source_file_description)

    def add_data_processing(self, data_processing_description):
        self.data_processing_list.append(data_processing_description)

    def add_instrument_configuration(self, instrument_description):
        self.instrument_configuration_list.append(instrument_description)

    def _create_file_description(self):
        self.writer.file_description(self.file_contents_list, self.source_file_list)

    def _create_software_list(self):
        self.writer.software_list([{
            "id": "ms_deisotope_1",
            "name": "ms_deisotope"
        }])

    def _add_spectrum_list(self):
        self._create_file_description()
        self._create_software_list()
        self.context_stack.append(self.writer.run(id=self.sample_name))
        self.context_stack[-1].__enter__()
        self.context_stack.append(self.writer.spectrum_list(count=self.n_spectra))
        self.context_stack[-1].__enter__()

    def _pack_activation(self, activation_information):
        params = []
        params.append({
            "name": str(activation_information.method),
        })
        # NOTE: Only correct for CID/HCD spectra with absolute collision energies, but that is all I have
        # to test with.
        params.append({
            "name": "collision energy",
            "value": activation_information.energy,
            "unitName": "electron volts"
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

    def _pack_precursor_information(self, precursor_information, activation_information=None):
        # If the scan bunch has been fully deconvoluted and it's PrecursorInformation
        # filled in, its extracted fields will be populated and should be used, otherwise
        # use the default read values.
        if precursor_information.extracted_neutral_mass != 0:
            package = {
                "mz": precursor_information.extracted_mz,
                "intensity": precursor_information.extracted_intensity,
                "charge": precursor_information.extracted_charge,
                "scan_id": precursor_information.precursor_scan_id
            }
        else:
            package = {
                "mz": precursor_information.mz,
                "intensity": precursor_information.intensity,
                "charge": precursor_information.charge,
                "scan_id": precursor_information.precursor_scan_id
            }
        if activation_information is not None:
            package['activation'] = self._pack_activation(activation_information)
        return package

    def _prepare_extra_arrays(self, scan):
        extra_arrays = []
        if self.deconvoluted:
            score_array = [
                peak.score for peak in scan.deconvoluted_peak_set
            ]
            extra_arrays.append((score_array, "deconvolution score array"))
            envelope_array = envelopes_to_array([peak.envelope for peak in scan.deconvoluted_peak_set])
            extra_arrays.append((envelope_array, "isotopic envelopes array"))
        return extra_arrays

    def save_scan_bunch(self, bunch, **kwargs):
        if not self._has_started_writing_spectra:
            self._add_spectrum_list()
            self._has_started_writing_spectra = True

        if self.deconvoluted:
            precursor_peaks = bunch.precursor.deconvoluted_peak_set
        else:
            precursor_peaks = bunch.precursor.peak_set

        polarity = bunch.precursor.polarity
        if self.deconvoluted:
            charge_array = [p.charge for p in precursor_peaks]
        else:
            charge_array = None

        descriptors = describe_spectrum(precursor_peaks)

        self.writer.write_spectrum(
            [p.mz for p in precursor_peaks], [p.intensity for p in precursor_peaks], charge_array,
            id=bunch.precursor.id, params=[
                {"name": "ms level", "value": bunch.precursor.ms_level},
                {"name": "MS1 spectrum"}] + descriptors,
            polarity=polarity,
            scan_start_time=bunch.precursor.scan_time,
            compression=self.compression,
            other_arrays=self._prepare_extra_arrays(bunch.precursor))

        self.total_ion_chromatogram_tracker[
            bunch.precursor.scan_time] = _total_intensity_from_descriptors(descriptors)
        self.base_peak_chromatogram_tracker[
            bunch.precursor.scan_time] = _base_peak_from_descriptors(descriptors)

        for prod in bunch.products:
            if self.deconvoluted:
                product_peaks = prod.deconvoluted_peak_set
            else:
                product_peaks = prod.peak_set

            descriptors = describe_spectrum(product_peaks)

            self.total_ion_chromatogram_tracker[
                prod.scan_time] = _total_intensity_from_descriptors(descriptors)
            self.base_peak_chromatogram_tracker[
                prod.scan_time] = _base_peak_from_descriptors(descriptors)

            if self.deconvoluted:
                charge_array = [p.charge for p in product_peaks]
            else:
                charge_array = None

            self.writer.write_spectrum(
                [p.mz for p in product_peaks], [p.intensity for p in product_peaks], charge_array,
                id=prod.id, params=[
                    {"name": "ms level", "value": prod.ms_level},
                    {"name": "MSn spectrum"}] + descriptors,
                polarity=prod.polarity,
                scan_start_time=prod.scan_time, precursor_information=self._pack_precursor_information(
                    prod.precursor_information, prod.activation),
                compression=self.compression,
                other_arrays=self._prepare_extra_arrays(prod))

        if self.indexer is not None:
            self.indexer.add_scan_bunch(bunch)

    def complete(self):
        for element in self.context_stack[::-1]:
            element.__exit__(None, None, None)
        self.context_stack = []
        if self.indexer is not None:
            try:
                name = self.handle.name
            except AttributeError:
                name = "_detatched_mzml_index"
            try:
                with open(ExtendedScanIndex.index_file_name(name), 'w') as ixfile:
                    self.indexer.serialize(ixfile)
            except IOError:
                pass


def marshal_deconvoluted_peak_set(scan_dict):
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


class ProcessedMzMLDeserializer(MzMLLoader, ScanDeserializerBase):
    def __init__(self, source_file, use_index=True):
        MzMLLoader.__init__(self, source_file, use_index=use_index)
        self.extended_index = ExtendedScanIndex()
        if self._use_index:
            try:
                if os.path.exists(self._index_file_name):
                    self.extended_index = ExtendedScanIndex.deserialize(
                        open(self._index_file_name))
            except IOError:
                pass
            except ValueError:
                pass

    def iter_scan_headers(self, iterator=None):
        if iterator is None:
            iterator = iter(self._source)
        precursor_scan = None
        product_scans = []

        _make_scan = super(ProcessedMzMLDeserializer, self)._make_scan

        current_level = 1
        for scan in iterator:
            if not self._validate(scan):
                continue
            packed = _make_scan(scan)
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
                raise Exception("This object is not able to handle MS levels higher than 2")
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)

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
            return self._scan_cache[scan_id]
        except KeyError:
            packed = super(ProcessedMzMLDeserializer, self)._make_scan(self._source.get_by_id(scan_id))
            return packed

    @property
    def _index_file_name(self):
        return ExtendedScanIndex.index_file_name(self.source_file)

    def build_extended_index(self):
        self.reset()
        indexer = ExtendedScanIndex()
        for bunch in self.iter_scan_headers():
            indexer.add_scan_bunch(bunch)
        self.reset()
        self.extended_index = indexer
        try:
            with open(self._index_file_name, 'w') as handle:
                indexer.serialize(handle)
        except (IOError, OSError, AttributeError) as err:
            print(err)
            pass

    def _make_scan(self, data):
        scan = super(ProcessedMzMLDeserializer, self)._make_scan(data)
        scan.peak_set = PeakIndex(np.array([]), np.array([]), PeakSet([]))
        scan.deconvoluted_peak_set = marshal_deconvoluted_peak_set(data)
        return scan.pack()

    def convert_scan_id_to_retention_time(self, scan_id):
        header = self.get_scan_header_by_id(scan_id)
        return header.scan_time

    # LC-MS/MS Database API

    def precursor_information(self):
        out = []
        for key, info in self.extended_index.msn_ids.items():
            mz = info['mz']
            neutral_mass = info['neutral_mass']
            charge = info['charge']
            intensity = info['intensity']
            precursor_scan_id = info['precursor_scan_id']
            product_scan_id = info['product_scan_id']
            pinfo = PrecursorInformation(
                mz, intensity, charge, precursor_scan_id,
                self, neutral_mass, charge, intensity,
                product_scan_id=product_scan_id)
            out.append(pinfo)
        return out

    def ms1_peaks_above(self, threshold=1000.):
        accumulate = []
        for ms1_id in self.extended_index.ms1_ids:
            scan = self.get_scan_by_id(ms1_id)
            for peak in scan.deconvoluted_peak_set:
                if peak.intensity < threshold:
                    continue
                accumulate.append((ms1_id, peak, id(peak)))
        return accumulate

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
