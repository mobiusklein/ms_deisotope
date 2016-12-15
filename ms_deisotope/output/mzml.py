from .common import ScanSerializerBase

from mzml_writer import writer


class MzMLScanSerializer(ScanSerializerBase):

    def __init__(self, handle, n_spectra=2e4, compression="zlib"):
        self.handle = handle
        self.writer = writer.MzMLWriter(handle)
        self.n_spectra = n_spectra
        self.compression = compression
        self._has_started_writing_spectra = False
        self.writer.__enter__()
        self.element_stack = [self.writer]
        self.writer.controlled_vocabularies()

    def _add_spectrum_list(self):
        self.element_stack.append(self.writer.element("run"))
        self.element_stack[-1].__enter__()
        self.element_stack.append(self.writer.element(
            "spectrumList", count=self.n_spectra))
        self.element_stack[-1].__enter__()

    def save_scan_bunch(self, bunch, **kwargs):
        if not self._has_started_writing_spectra:
            self._add_spectrum_list()
            self._has_started_writing_spectra = True
        peaks = bunch.precursor.deconvoluted_peak_set
        self.writer.write_spectrum(
            [p.mz for p in peaks], [p.intensity for p in peaks], [
                p.charge for p in peaks],
            id=bunch.precursor.id, params=[
                {"name": "ms level", "value": bunch.precursor.ms_level}],
            polarity=peaks[0].charge,
            scan_start_time=bunch.precursor.scan_time,
            compression=self.compression)
        for prod in bunch.products:
            peaks = prod.deconvoluted_peak_set
            self.writer.write_spectrum(
                [p.mz for p in peaks], [p.intensity for p in peaks], [
                    p.charge for p in peaks],
                id=prod.id, params=[{"name": "ms level", "value": prod.ms_level}],
                polarity=peaks[0].charge,
                scan_start_time=prod.scan_time, precursor_information={
                    "mz": prod.precursor_information.extracted_mz,
                    "intensity": prod.precursor_information.extracted_intensity,
                    "charge": prod.precursor_information.extracted_charge,
                    "scan_id": prod.precursor_information.precursor_scan_id
                },
                compression=self.compression)

    def complete(self):
        for element in self.element_stack[::-1]:
            element.__exit__()
