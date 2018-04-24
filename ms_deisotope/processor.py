import logging

from ms_peak_picker import pick_peaks

from .deconvolution import deconvolute_peaks
from .data_source.infer_type import MSFileLoader
from .data_source.common import Scan, ScanBunch, ChargeNotProvided
from .utils import Base, LRUDict
from .peak_dependency_network import NoIsotopicClustersError
from .scoring import InterferenceDetection

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


class PriorityTarget(Base):
    """Represent a targeted envelope deconvolution's parameters and constraints.

    This class is used to tell :func:`ms_deisotope.deconvolution.deconvolute_peaks`
    that the solution produced by this peak should be preferentially extracted.

    Attributes
    ----------
    info : PrecursorInformation
        The associated precursor information block which contains
        the charge state hint.
    peak : FittedPeak
        The peak from which to start the deconvolution
    trust_charge_hint : bool
        Whether or not to force the deconvoluter to only consider
        the charge specified in the hint.
    mz : float
        The m/z of :attr:`peak`
    charge : int
        The charge state hint from :attr:`info`
    isolation_window : IsolationWindow
        The isolation window for this precursor ion. May be `None`
    """

    def __init__(self, peak, info, trust_charge_hint=True, precursor_scan_id=None,
                 product_scan_id=None, isolation_window=None):
        self.peak = peak
        self.info = info
        self.trust_charge_hint = trust_charge_hint
        self.precursor_scan_id = precursor_scan_id
        self.product_scan_id = product_scan_id
        self.isolation_window = isolation_window

    def __iter__(self):
        yield self.peak
        yield self.info

    @property
    def mz(self):
        return self.peak.mz

    @property
    def charge(self):
        return int(self.info.charge)

    def charge_range_hint(self, charge_range):
        """Create an updated charge range for a Deconvoluter to search.

        At the moment, this only amounts to either returning the charge
        range unchanged or returning a charge range that only contains
        the hinted charge state, depending upon whether :attr:`trust_charge_hint`
        is `False` or not.

        Parameters
        ----------
        charge_range : tuple
            The charge range to update

        Returns
        -------
        tuple
            The updated charge range
        """
        if self.trust_charge_hint and self.info.charge is not ChargeNotProvided:
            return (self.charge, self.charge)
        else:
            return charge_range

    def __repr__(self):
        return "PriorityTarget(mz=%0.4f, intensity=%0.4f, charge_hint=%d)" % (
            self.mz, self.peak.intensity, self.charge)


class ScanProcessor(Base):
    """Orchestrates the deconvolution of a `ScanIterator` scan by scan. This process will
    apply different rules for MS1 scans and MSn scans. This type itself is an Iterator,
    consuming (raw) mass spectral data and producing deisotoped and charge deconvolved spectra.

    The algorithms used for each task are independent and can be specified in the appropriate
    attribute dictionary, however there is information sharing between each MS1 scan and its
    MSn scans as the precursor monoisotopic mass is recalibrated according to the MS1 processing
    arguments, and the selected charge state is used to limit the charge range used in the matching
    MSn scan. These are described by :class:`PriorityTarget` objects.

    At the moment, MSn assumes only MS2. Until MS3 data become available for testing, this limit
    will remain.

    Attributes
    ----------
    data_source : str or file-like
        Any valid object to be passed to the `loader_type` callable to produce
        a :class:`ScanIterator` instance. A path to an mzML file will work for the
        default loader. Used to populate :attr:`reader`
    loader_type : callable
        A callable, which when passed `data_source` returns an instance of :class:`ScanIterator`.
        By default, this is :class:`MSFileLoader`. Used to populate :attr:`reader`
    reader: ScanIterator
        Any object implementing the :class:`ScanIterator` interface, produced by calling
        :attr:`loader_type` on :attr:`data_source`.
    ms1_deconvolution_args : dict
        The arguments passed to :func:`ms_deisotope.deconvolution.deconvolute_peaks` for MS1
        scans.
    ms1_peak_picking_args : dict
        The arguments passed to :func:`ms_peak_picker.pick_peaks` for MS1 scans.
    msn_deconvolution_args : dict
        The arguments passed to :func:`ms_deisotope.deconvolution.deconvolute_peaks` for MSn
        scans.
    msn_peak_picking_args : dict
        The arguments passed to :func:`ms_peak_picker.pick_peaks` for MSn scans.
    pick_only_tandem_envelopes : bool
        Whether or not to process whole MS1 scans or just the regions around those peaks
        chosen for MSn
    default_precursor_ion_selection_window : float
        Size of the selection window to use when `pick_only_tandem_envelopes` is `True`
        and the information is not available in the scan.
    trust_charge_hint : bool
        Whether or not to trust the charge provided by the data source when determining
        the charge state of precursor isotopic patterns. Defaults to `True`
    respect_isolation_window: bool
        Whether to use the bounds of the isolation window to reject a monoisotopic peak
        solution
    terminate_on_error: bool
        Whether or not  to stop processing on an error. Defaults to `True`
    """

    def __init__(self, data_source, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None,
                 ms1_deconvolution_args=None,
                 msn_deconvolution_args=None,
                 pick_only_tandem_envelopes=False,
                 default_precursor_ion_selection_window=1.5,
                 trust_charge_hint=True,
                 loader_type=None,
                 envelope_selector=None,
                 terminate_on_error=True,
                 ms1_averaging=0,
                 respect_isolation_window=False):
        if loader_type is None:
            loader_type = MSFileLoader

        self.data_source = data_source
        self.ms1_peak_picking_args = ms1_peak_picking_args or {}
        self.msn_peak_picking_args = msn_peak_picking_args or ms1_peak_picking_args or {}
        self.ms1_deconvolution_args = ms1_deconvolution_args or {}
        self.ms1_deconvolution_args.setdefault("charge_range", (1, 8))
        self.msn_deconvolution_args = msn_deconvolution_args or {}
        self.msn_deconvolution_args.setdefault("charge_range", (1, 8))
        self.pick_only_tandem_envelopes = pick_only_tandem_envelopes

        self.default_precursor_ion_selection_window = default_precursor_ion_selection_window
        self.respect_isolation_window = respect_isolation_window

        self.trust_charge_hint = trust_charge_hint
        self.ms1_averaging = int(ms1_averaging) if ms1_averaging else 0

        self.loader_type = loader_type

        self._signal_source = self.loader_type(data_source)
        self.envelope_selector = envelope_selector
        self.terminate_on_error = terminate_on_error

        self._ms1_index_cache = LRUDict(maxsize=self.ms1_averaging * 2 + 2)

    def _reject_candidate_precursor_peak(self, peak, product_scan):
        isolation = product_scan.isolation_window
        if isolation is None or isolation.is_empty():
            pinfo = product_scan.precursor_information
            err = peak.mz - pinfo.mz
            return abs(err) > self.default_precursor_ion_selection_window
        else:
            return peak.mz not in isolation and self.respect_isolation_window

    @property
    def reader(self):
        return self._signal_source

    def _get_envelopes(self, precursor_scan):
        """Get the m/z intervals to pick peaks from for the
        given MS1 scan

        Parameters
        ----------
        precursor_scan: Scan

        Returns
        -------
        list or None
        """
        if not self.pick_only_tandem_envelopes and self.envelope_selector is None:
            return None
        elif self.envelope_selector is None:
            chosen_envelopes = [s.precursor_information for s in precursor_scan.product_scans]
            chosen_envelopes = sorted([(p.mz - 5, p.mz + 10) for p in chosen_envelopes])
        else:
            chosen_envelopes = self.envelope_selector(precursor_scan)
        return chosen_envelopes

    def _pick_precursor_scan_peaks(self, precursor_scan, chosen_envelopes=None):
        """Pick peaks from the given precursor scan

        Parameters
        ----------
        precursor_scan: Scan
            Scan to pick peaks from
        chosen_envelopes: list, optional
            list of m/z intervals to pick peaks for

        Returns
        -------
        PeakSet
        """
        if precursor_scan.is_profile:
            peak_mode = 'profile'
        else:
            peak_mode = 'centroid'
        prec_mz, prec_intensity = precursor_scan.arrays
        if not self.pick_only_tandem_envelopes and self.envelope_selector is None:
            prec_peaks = pick_peaks(prec_mz, prec_intensity, peak_mode=peak_mode, **self.ms1_peak_picking_args)
        else:
            if chosen_envelopes is None:
                chosen_envelopes = self._get_envelopes(precursor_scan)
            prec_peaks = pick_peaks(prec_mz, prec_intensity, peak_mode=peak_mode,
                                    target_envelopes=chosen_envelopes,
                                    **self.ms1_peak_picking_args)
        return prec_peaks

    def _average_ms1(self, precursor_scan):
        """Average signal from :attr:`self.ms1_averaging` scans from
        before and after ``precursor_scan`` and pick peaks from the
        averaged arrays.

        Parameters
        ----------
        precursor_scan: Scan
            The scan to use as a point of reference

        Returns
        -------
        PeakSet
        """
        new_scan = precursor_scan.average(self.ms1_averaging)
        prec_peaks = pick_peaks(*new_scan.arrays,
                                target_envelopes=self._get_envelopes(precursor_scan),
                                **self.ms1_peak_picking_args)
        return prec_peaks

    def pick_precursor_scan_peaks(self, precursor_scan):
        """Picks peaks for the given ``precursor_scan`` using the
        appropriate strategy.

        If :attr:`ms1_averaging` > 0, then the signal averaging strategy
        is used, otherwise peaks are picked directly.

        Parameters
        ----------
        precursor_scan: Scan

        Returns
        -------
        PeakSet
        """
        logger.info("Picking Precursor Scan Peaks: %r", precursor_scan)
        if self.ms1_averaging > 0:
            prec_peaks = self._average_ms1(precursor_scan)
        else:
            prec_peaks = self._pick_precursor_scan_peaks(precursor_scan)
        precursor_scan.peak_set = prec_peaks
        return prec_peaks

    def pick_product_scan_peaks(self, product_scan):
        if product_scan.is_profile:
            peak_mode = 'profile'
        else:
            peak_mode = 'centroid'
        product_mz, product_intensity = product_scan.arrays
        peaks = pick_peaks(product_mz, product_intensity, peak_mode=peak_mode, **self.msn_peak_picking_args)

        if peaks is None:
            raise EmptyScanError(
                "Could not pick peaks for empty product scan", self)

        product_scan.peak_set = peaks
        return peaks

    def get_precursor_peak_for_product_scans(self, precursor_scan):  # pragma: no cover
        priorities = []
        peaks = precursor_scan.peak_set
        for scan in precursor_scan.product_scans:
            precursor_ion = scan.precursor_information
            if peaks is None:
                peaks = self.pick_precursor_scan_peaks(precursor_scan)
            peak, err = peaks.get_nearest_peak(precursor_ion.mz)
            precursor_ion.peak = peak
            target = PriorityTarget(
                peak,
                precursor_ion,
                self.trust_charge_hint,
                isolation_window=scan.isolation_window)
            if self._reject_candidate_precursor_peak(peak, scan):
                logger.info(
                    "Unable to locate a peak for precursor ion %r for tandem scan %s of precursor scan %s" % (
                        precursor_ion, scan.title,
                        precursor_scan.title))
            else:
                priorities.append(target)
        return priorities

    def process_scan_group(self, precursor_scan, product_scans):
        """Performs the initial extraction of information relating
        `precursor_scan` to `product_scans` and picks peaks for ``precursor_scan``.
        Called by :meth:`process`. May be used separately if doing the process step
        by step.

        Parameters
        ----------
        precursor_scan : Scan
            An MS1 Scan
        product_scans : list of Scan
            A list of MSn Scans related to `precursor_scan`

        Returns
        -------
        precursor_scan: Scan
            As Parameter
        prioritiies: list of PriorityTarget
            list of the peak target windows in `precursor_scan` which
            are related to `product_scans`
        product_scans: list of Scan
            As Parameter
        """
        prec_peaks = self.pick_precursor_scan_peaks(precursor_scan)
        priorities = []

        if prec_peaks is None:
            raise EmptyScanError(
                "Could not pick peaks for empty precursor scan", self)

        for scan in product_scans:
            precursor_ion = scan.precursor_information
            peak, err = prec_peaks.get_nearest_peak(precursor_ion.mz)
            precursor_ion.peak = peak
            target = PriorityTarget(
                peak, precursor_ion, self.trust_charge_hint,
                scan.precursor_information.precursor_scan_id,
                scan.precursor_information.product_scan_id,
                isolation_window=scan.isolation_window)
            if self._reject_candidate_precursor_peak(peak, scan):
                logger.info(
                    "Unable to locate a peak for precursor ion %r for tandem scan %s of precursor scan %s" % (
                        precursor_ion, scan.id,
                        precursor_scan.id))
            else:
                priorities.append(target)

        return precursor_scan, priorities, product_scans

    def _default_all_precursor_information(self, scans):
        for scan in scans:
            if scan.ms_level > 1:
                scan.precursor_information.default(orphan=True)

    def deconvolute_precursor_scan(self, precursor_scan, priorities=None):
        if priorities is None:
            priorities = []

        logger.info("Deconvoluting Precursor Scan %r", precursor_scan)
        logger.info("Priorities: %r", priorities)

        ms1_deconvolution_args = self.ms1_deconvolution_args.copy()

        if precursor_scan.polarity in (1, -1):
            polarity = precursor_scan.polarity
            ms1_deconvolution_args['charge_range'] = tuple(
                polarity * abs(c) for c in ms1_deconvolution_args['charge_range'])
        try:
            decon_result = deconvolute_peaks(
                precursor_scan.peak_set, priority_list=priorities,
                **ms1_deconvolution_args)
        except NoIsotopicClustersError as e:
            e.scan_id = precursor_scan.id
            if self.terminate_on_error:
                raise e
            else:
                logger.warn("No isotopic clusters found in %r" % precursor_scan.id)

        dec_peaks, priority_results = decon_result

        if decon_result.errors:
            logger.error("Errors occurred during deconvolution of %s, %r" % (
                precursor_scan.id, decon_result.errors))

        for pr in priority_results:
            if pr is None:
                continue
            else:
                pr.chosen_for_msms = True

        # `priorities` and `priority_results` are parallel lists. The
        # ith position in `priorities` corresponds to the ith deconvoluted
        # priority result in `priority_results`. The entry in `priority_results`
        # may be `None` if the deconvolution failed, but elements of `priorities`
        # should always be FittedPeak or Peak-like instances

        interference_detector = InterferenceDetection(precursor_scan.peak_set)
        for product_scan in precursor_scan.product_scans:
            precursor_information = product_scan.precursor_information

            isolation_window = product_scan.isolation_window
            if isolation_window.is_empty():
                lower = precursor_information.mz - self.default_precursor_ion_selection_window
                upper = precursor_information.mz + self.default_precursor_ion_selection_window
            else:
                lower = isolation_window.lower_bound
                upper = isolation_window.upper_bound

            # unknown precursor purity
            product_scan.annotations['precursor purity'] = 0.0
            i = get_nearest_index(precursor_information.mz, priorities)

            # If no peak is found in the priority list, it means the priority list is empty.
            # This should never happen in the current implementation. If it did, then we forgot
            # to pass the priority list to this function.
            if i is None:
                logger.warning(
                    "Could not find deconvolution for %r (No nearby peak in the priority list)",
                    precursor_information)
                precursor_information.default(orphan=True)
                continue

            peak = priority_results[i]
            # If the deconvolution result is None, then we have no answer
            if peak is None:
                logger.warning(
                    "Could not find deconvolution for %r (No solution was found for this region)",
                    precursor_information)
                precursor_information.default(orphan=True)

                continue
            elif peak.charge == 1 or (peak.charge != precursor_information.charge and self.trust_charge_hint):
                if precursor_information.charge != ChargeNotProvided:
                    logger.warning(
                        "Could not find deconvolution for %r (Unacceptable solution was proposed: %r)",
                        precursor_information, peak)
                    precursor_information.default()
                    continue

            precursor_purity = 1 - interference_detector.detect_interference(peak.envelope, lower, upper)
            product_scan.annotations['precursor purity'] = precursor_purity
            precursor_information.extract(peak)
        precursor_scan.deconvoluted_peak_set = dec_peaks
        return dec_peaks, priority_results

    def deconvolute_product_scan(self, product_scan):
        logger.info("Deconvoluting Product Scan %r", product_scan)
        precursor_ion = product_scan.precursor_information
        top_charge_state = precursor_ion.extracted_charge
        deconargs = dict(self.msn_deconvolution_args)
        charge_range = list(deconargs.get("charge_range", [1, top_charge_state]))
        if top_charge_state is not None and top_charge_state is not ChargeNotProvided and\
           top_charge_state != 0 and abs(top_charge_state) < abs(charge_range[1]):
            charge_range[1] = top_charge_state

        deconargs["charge_range"] = charge_range

        if product_scan.polarity in (-1, 1):
            polarity = product_scan.polarity
            deconargs["charge_range"] = [
                polarity * abs(c) for c in deconargs["charge_range"]]

        try:
            dec_peaks, _ = deconvolute_peaks(product_scan.peak_set, **deconargs)
        except NoIsotopicClustersError as e:
            logger.info("No Isotopic Clusters found in %r" % product_scan.id)
            e.scan_id = product_scan.id
            if self.terminate_on_error:
                raise e

        product_scan.deconvoluted_peak_set = dec_peaks
        return dec_peaks

    def _get_next_scans(self):
        bunch = next(self.reader)
        try:
            precursor, products = bunch
        except ValueError:
            if isinstance(bunch, Scan):
                if bunch.ms_level == 1:
                    precursor = bunch
                    products = []
                else:
                    precursor = None
                    products = [bunch]

        if self.pick_only_tandem_envelopes:
            while len(products) == 0:
                precursor, products = next(self.reader)

        return precursor, products

    def process(self, precursor, products):
        """Fully preprocesses the `precursor` and `products` scans, performing
        any necessary information sharing.

        Parameters
        ----------
        precursor : Scan
            An MS1 Scan
        products : list of Scan
            A list of MSn Scans related to `precursor`

        Returns
        -------
        precursor_scan: Scan
            The fully processed version of `precursor`
        product_scans: list of Scan
            The fully processed version of `products`
        """
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        if precursor_scan is not None:
            self.deconvolute_precursor_scan(precursor_scan, priorities)
        else:
            self._default_all_precursor_information(product_scans)

        for product_scan in product_scans:
            self.pick_product_scan_peaks(product_scan)
            self.deconvolute_product_scan(product_scan)

        return precursor_scan, product_scans

    def next(self):
        """Fetches the next bunch of scans from :attr:`reader` and
        invokes :meth:`process` on them, picking peaks and deconvoluting them.

        Returns
        -------
        ScanBunch
        """
        precursor, products = self._get_next_scans()
        precursor_scan, product_scans = self.process(precursor, products)
        return ScanBunch(precursor_scan, product_scans)

    def __next__(self):
        """Fetches the next bunch of scans from :attr:`reader` and
        invokes :meth:`process` on them, picking peaks and deconvoluting them.

        Returns
        -------
        ScanBunch
        """
        return self.next()

    def __iter__(self):
        return self

    def pack_next(self):
        """As :meth:`next`, except instead of producing :class:`ScanBunch` of
        :class:`Scan` instances, instead it uses :class:`ProcessedScan` to strip away
        much of the heavy information like the raw data arrays.

        Returns
        -------
        ScanBunch
        """
        precursor, products = self._get_next_scans()
        precursor_scan, product_scans = self.process(precursor, products)
        return ScanBunch(precursor_scan.pack() if precursor_scan else None, [p.pack() for p in product_scans])

    def start_from_scan(self, *args, **kwargs):
        """A wrapper around :meth:`start_from_scan` provided by
        :attr:`reader`.

        Returns
        -------
        self

        See Also
        --------
        ms_deisotope.data_source.mzml.start_from_scan
        """
        self.reader.start_from_scan(*args, **kwargs)
        return self


class EmptyScanError(ValueError):
    def __init__(self, msg, scan_id=None):
        ValueError.__init__(self, msg)
        self.scan_id = scan_id
