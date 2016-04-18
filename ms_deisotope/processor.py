import warnings
import logging

from ms_peak_picker import pick_peaks

from .deconvolution import deconvolute_peaks
from .data_source.mzml import PyteomicsMzMLLoader
from .data_source.common import ScanBunch


logger = logging.getLogger("deconvolution_scan_processor")


def get_nearest_index(query_mz, peak_list):
    best_index = None
    best_error = float('inf')

    for i, peak in enumerate(peak_list):
        error = abs(peak.mz - query_mz)
        if error < best_error:
            best_error = error
            best_index = i
    return best_index


class ScanProcessor(object):
    def __init__(self, mzml_file, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 pick_only_tandem_envelopes=False, precursor_ion_mass_accuracy=2e-5):
        self.mzml_file = mzml_file
        self.ms1_peak_picking_args = ms1_peak_picking_args or {}
        self.msn_peak_picking_args = msn_peak_picking_args or ms1_peak_picking_args or {}
        self.ms1_deconvolution_args = ms1_deconvolution_args or {}
        self.msn_deconvolution_args = msn_deconvolution_args or {}
        self.pick_only_tandem_envelopes = pick_only_tandem_envelopes
        self.precursor_ion_mass_accuracy = precursor_ion_mass_accuracy

        self._signal_source = PyteomicsMzMLLoader(mzml_file)

    def pick_precursor_scan_peaks(self, precursor_scan):
        logger.info("Picking Precursor Scan Peaks: %r", precursor_scan)
        prec_mz, prec_intensity = self._signal_source.scan_arrays(precursor_scan)
        if not self.pick_only_tandem_envelopes:
            prec_peaks = pick_peaks(prec_mz, prec_intensity, **self.ms1_peak_picking_args)
        else:
            chosen_envelopes = [s.precursor_information for s in precursor_scan.product_scans]
            chosen_envelopes = sorted([(p.mz - 5, p.mz + 10) for p in chosen_envelopes])
            prec_peaks = pick_peaks(prec_mz, prec_intensity, target_envelopes=chosen_envelopes,
                                    **self.ms1_peak_picking_args)

        precursor_scan.peak_set = prec_peaks
        return prec_peaks

    def pick_product_scan_peaks(self, product_scan):
        logger.info("Picking Product Scan Peaks: %r", product_scan)
        product_mz, product_intensity = self._signal_source.scan_arrays(product_scan)
        peaks = pick_peaks(product_mz, product_intensity, **self.msn_peak_picking_args)
        product_scan.peak_set = peaks
        return peaks

    def get_precursor_peak_for_product_scans(self, precursor_scan):
        priorities = []
        peaks = precursor_scan.peak_set
        for scan in precursor_scan.product_scans:
            precursor_ion = scan.precursor_information
            if peaks is None:
                peaks = self.pick_precursor_scan_peaks(precursor_scan)
            peak, err = peaks.get_nearest_peak(precursor_ion.mz)
            precursor_ion.peak = peak
            if abs(err) > 0.5:
                warnings.warn(
                    "Unable to locate a peak for precursor ion %r for tandem scan %s of precursor scan %s" % (
                        precursor_ion, self._signal_source.scan_title(scan),
                        self._signal_source.scan_title(precursor_scan)))
            else:
                priorities.append(peak)
        return priorities

    def process_scan_group(self, precursor_scan, product_scans):
        prec_peaks = None
        priorities = []
        product_peaks = []

        for scan in product_scans:
            precursor_ion = scan.precursor_information
            if prec_peaks is None:
                prec_peaks = self.pick_precursor_scan_peaks(precursor_scan)
            peak, err = prec_peaks.get_nearest_peak(precursor_ion.mz)
            precursor_ion.peak = peak
            if abs(err) > 0.5:
                warnings.warn(
                    "Unable to locate a peak for precursor ion %r for tandem scan %s of precursor scan %s" % (
                        precursor_ion, self._signal_source.scan_title(scan),
                        self._signal_source.scan_title(precursor_scan)))
            else:
                priorities.append(peak)
            self.pick_product_scan_peaks(scan)
            product_peaks.append(scan)
        if prec_peaks is None:
            prec_peaks = self.pick_precursor_scan_peaks(precursor_scan)

        return precursor_scan, priorities, product_peaks

    def deconvolute_precursor_scan(self, precursor_scan, priorities=None):
        if priorities is None:
            priorities = []

        logger.info("Deconvoluting Precursor Scan %r: %r", precursor_scan, priorities)

        dec_peaks, priority_results = deconvolute_peaks(
            precursor_scan.peak_set, priority_list=priorities,
            **self.ms1_deconvolution_args)

        for product_scan in precursor_scan.product_scans:
            precursor_information = product_scan.precursor_information

            i = get_nearest_index(precursor_information.mz, priorities)
            peak = priority_results[i]
            if peak is None:
                warnings.warn("Could not find deconvolution for %r" % precursor_information)
                precursor_information.extracted_neutral_mass = precursor_information.neutral_mass
                precursor_information.extracted_charge = int(precursor_information.charge)
                precursor_information.extracted_peak_height = precursor_information.peak.intensity

                continue
            elif peak.charge == 1:
                precursor_information.extracted_neutral_mass = precursor_information.neutral_mass
                precursor_information.extracted_charge = int(precursor_information.charge)
                precursor_information.extracted_peak_height = precursor_information.peak.intensity
                continue

            precursor_information.extracted_neutral_mass = peak.neutral_mass
            precursor_information.extracted_charge = peak.charge
            precursor_information.extracted_peak_height = peak.intensity
            precursor_information.extracted_peak = peak

        return dec_peaks

    def deconvolute_product_scan(self, product_scan):
        logger.info("Deconvoluting Product Scan %r", product_scan)
        precursor_ion = product_scan.precursor_information
        top_charge_state = precursor_ion.extracted_charge
        deconargs = dict(self.msn_deconvolution_args)
        charge_range = list(deconargs.get("charge_range", [1, top_charge_state]))
        if top_charge_state is not None and top_charge_state != 0:
            charge_range[1] = top_charge_state

        deconargs["charge_range"] = charge_range

        dec_peaks, _ = deconvolute_peaks(product_scan.peak_set, **deconargs)
        return dec_peaks

    def _next(self):
        precursor, products = self._signal_source.next()

        if self.pick_only_tandem_envelopes:
            while len(products) == 0:
                precursor, products = self._signal_source.next()

        return precursor, products

    def _process(self, precursor, products):
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        deconvoluted_precursor_peaks = self.deconvolute_precursor_scan(precursor_scan, priorities)
        precursor_scan.deconvoluted_peak_set = deconvoluted_precursor_peaks

        for product_scan in product_scans:
            peaks = self.deconvolute_product_scan(product_scan)
            product_scan.deconvoluted_peak_set = peaks
            pass

        return precursor_scan, product_scans

    def next(self):
        precursor, products = self._next()
        precursor_scan, product_scans = self._process(precursor, products)
        return ScanBunch(precursor_scan, product_scans)

    def seek_precursor_scan(self, scan_id):
        precursor_scan, product_scans = self._next()
        while precursor_scan.id != scan_id:
            precursor_scan, product_scans = self._next()
        precursor_scan, product_scans = self._process(precursor_scan, product_scans)
        return ScanBunch(precursor_scan, product_scans)

    def seek_while(self, condition):
        precursor_scan, product_scans = self._next()
        while not condition(precursor_scan, product_scans):
            precursor_scan, product_scans = self._next()
        precursor_scan, product_scans = self._process(precursor_scan, product_scans)
        return ScanBunch(precursor_scan, product_scans)

    def pack_next(self):
        precursor, products = self._next()
        precursor_scan, product_scans = self._process(precursor, products)
        return ScanBunch(precursor_scan.pack(), [p.pack() for p in product_scans])


class PrecursorSkimmingProcessor(ScanProcessor):
    def _process(self, precursor, products):
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        deconvoluted_precursor_peaks = self.deconvolute_precursor_scan(precursor_scan, priorities)
        precursor_scan.deconvoluted_peak_set = deconvoluted_precursor_peaks
        return precursor_scan, product_scans
