import warnings
import logging

from ms_peak_picker import pick_peaks

from .deconvolution import deconvolute_peaks
from .data_source.mzml import MzMLLoader
from .data_source.common import ScanBunch


PyteomicsMzMLLoader = MzMLLoader

logger = logging.getLogger("deconvolution_scan_processor")


def identity(x):
    return x


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
    def __init__(self, data_source, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 pick_only_tandem_envelopes=False, precursor_ion_mass_accuracy=2e-5,
                 loader_type=MzMLLoader):
        if loader_type is None:
            loader_type = MzMLLoader

        self.data_source = data_source
        self.ms1_peak_picking_args = ms1_peak_picking_args or {}
        self.msn_peak_picking_args = msn_peak_picking_args or ms1_peak_picking_args or {}
        self.ms1_deconvolution_args = ms1_deconvolution_args or {}
        self.msn_deconvolution_args = msn_deconvolution_args or {}
        self.pick_only_tandem_envelopes = pick_only_tandem_envelopes
        self.precursor_ion_mass_accuracy = precursor_ion_mass_accuracy

        self.loader_type = loader_type

        self._signal_source = self.loader_type(data_source)

    def pick_precursor_scan_peaks(self, precursor_scan):
        logger.info("Picking Precursor Scan Peaks: %r", precursor_scan)
        prec_mz, prec_intensity = precursor_scan.arrays
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
        product_mz, product_intensity = product_scan.arrays
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
                        precursor_ion, scan.title,
                        precursor_scan.title))
            else:
                priorities.append(peak)
        return priorities

    def process_scan_group(self, precursor_scan, product_scans):
        prec_peaks = None
        priorities = []

        for scan in product_scans:
            precursor_ion = scan.precursor_information
            if prec_peaks is None:
                prec_peaks = self.pick_precursor_scan_peaks(precursor_scan)
            peak, err = prec_peaks.get_nearest_peak(precursor_ion.mz)
            precursor_ion.peak = peak
            if abs(err) > 0.5:
                warnings.warn(
                    "Unable to locate a peak for precursor ion %r for tandem scan %s of precursor scan %s" % (
                        precursor_ion, scan.id,
                        precursor_scan.id))
            else:
                priorities.append(peak)
        if prec_peaks is None:
            prec_peaks = self.pick_precursor_scan_peaks(precursor_scan)

        return precursor_scan, priorities, product_scans

    def deconvolute_precursor_scan(self, precursor_scan, priorities=None):
        if priorities is None:
            priorities = []

        logger.info("Deconvoluting Precursor Scan %r", precursor_scan)
        logger.info("Priorities: %r", priorities)

        dec_peaks, priority_results = deconvolute_peaks(
            precursor_scan.peak_set, priority_list=priorities,
            **self.ms1_deconvolution_args)

        for product_scan in precursor_scan.product_scans:
            precursor_information = product_scan.precursor_information

            i = get_nearest_index(precursor_information.mz, priorities)
            if i is None:
                logger.info("Could not find deconvolution for %r (1)" % precursor_information)
                logger.info("%f, %r, %r", precursor_information.mz, priority_results, i)
                precursor_information.default()
                continue

            peak = priority_results[i]
            if peak is None:
                logger.info("Could not find deconvolution for %r (2)" % precursor_information)
                logger.info("%r, %r, %r, %r", precursor_information, priority_results, i, peak)
                precursor_information.default()

                continue
            elif peak.charge == 1:
                precursor_information.default()
                continue

            precursor_information.extract(peak)
        precursor_scan.deconvoluted_peak_set = dec_peaks
        return dec_peaks, priority_results

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
        product_scan.deconvoluted_peak_set = dec_peaks
        return dec_peaks

    def _get_next_scans(self):
        precursor, products = next(self._signal_source)

        if self.pick_only_tandem_envelopes:
            while len(products) == 0:
                precursor, products = next(self._signal_source)

        return precursor, products

    def _process(self, precursor, products):
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        self.deconvolute_precursor_scan(precursor_scan, priorities)

        for product_scan in product_scans:
            self.pick_product_scan_peaks(product_scan)
            self.deconvolute_product_scan(product_scan)

        return precursor_scan, product_scans

    def next(self):
        precursor, products = self._get_next_scans()
        precursor_scan, product_scans = self._process(precursor, products)
        return ScanBunch(precursor_scan, product_scans)

    def __next__(self):
        return self.next()

    def seek_precursor_scan(self, scan_id):
        precursor_scan, product_scans = self._get_next_scans()
        while precursor_scan.id != scan_id:
            precursor_scan, product_scans = self._get_next_scans()
        precursor_scan, product_scans = self._process(precursor_scan, product_scans)
        return ScanBunch(precursor_scan, product_scans)

    def seek_while(self, condition):
        precursor_scan, product_scans = self._get_next_scans()
        while not condition(precursor_scan, product_scans):
            precursor_scan, product_scans = self._get_next_scans()
        precursor_scan, product_scans = self._process(precursor_scan, product_scans)
        return ScanBunch(precursor_scan, product_scans)

    def pack_next(self):
        precursor, products = self._get_next_scans()
        precursor_scan, product_scans = self._process(precursor, products)
        return ScanBunch(precursor_scan.pack(), [p.pack() for p in product_scans])


class PrecursorSkimmingProcessor(ScanProcessor):
    def _process(self, precursor, products):
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        self.deconvolute_precursor_scan(precursor_scan, priorities)
        return precursor_scan, product_scans
