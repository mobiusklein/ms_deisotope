from .common import ScanSerializerBase

from psims.mzml import writer


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


class MzMLScanSerializer(ScanSerializerBase):

    def __init__(self, handle, n_spectra=2e4, compression="zlib", deconvoluted=True):
        self.handle = handle
        self.writer = writer.MzMLWriter(handle)
        self.n_spectra = n_spectra
        self.compression = compression
        self._has_started_writing_spectra = False
        self.writer.__enter__()
        self.context_stack = [self.writer]
        self.writer.controlled_vocabularies()
        self.deconvoluted = deconvoluted

    def _add_spectrum_list(self):
        self.context_stack.append(self.writer.run())
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
        for key, val in activation_information.items():
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
            compression=self.compression)

        for prod in bunch.products:
            if self.deconvoluted:
                product_peaks = prod.deconvoluted_peak_set
            else:
                product_peaks = prod.peak_set

            descriptors = describe_spectrum(product_peaks)

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
                compression=self.compression)

    def complete(self):
        for element in self.context_stack[::-1]:
            element.__exit__(None, None, None)
