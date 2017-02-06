import logging

from ms_peak_picker import pick_peaks

from .deconvolution import deconvolute_peaks
from .data_source.mzml import MzMLLoader
from .data_source.infer_type import MSFileLoader
from .data_source.common import ScanBunch
from .utils import Base
from .peak_dependency_network import Interval, IntervalTreeNode

PyteomicsMzMLLoader = MzMLLoader

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
    """Represent a targeted deconvolution's parameters and constraints.

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
    """

    def __init__(self, peak, info, trust_charge_hint=True):
        self.peak = peak
        self.info = info
        self.trust_charge_hint = trust_charge_hint

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
        if self.trust_charge_hint:
            return (self.charge, self.charge)
        else:
            return charge_range

    def __repr__(self):
        return "PriorityTarget(mz=%0.4f, intensity=%0.4f, charge_hint=%d)" % (
            self.mz, self.peak.intensity, self.charge)


class ScanProcessor(Base):
    """Orchestrates the deconvolution of a `ScanIterator` scan by scan. This process will
    apply different rules for MS^1 scans and MS^n scans. This type itself is an Iterator,
    consuming (raw) mass spectral data and producing deisotoped and charge deconvolved spectra.

    The algorithms used for each task are independent and can be specified in the appropriate
    attribute dictionary, however there is information sharing between each MS^1 scan and its
    MS^n scans as the precursor monoisotopic mass is recalibrated according to the MS^1 processing
    arguments, and the selected charge state is used to limit the charge range used in the matching
    MS^n scan. These are described by :class:`PriorityTarget` objects.

    At the moment, MS^n assumes only MS^2. Until MS^3 data become available for testing, this limit
    will remain.

    Attributes
    ----------
    data_source : str or file-like
        Any valid object to be passed to the `loader_type` callable to produce
        a :class:`ScanIterator` instance. A path to an mzML file will work for the
        default loader. Used to populate :attr:`reader`
    loader_type : callable
        A callable, which when passed `data_source` returns an instance of :class:`ScanIterator`.
        By default, this is :class:`MzMLLoader`. Used to populate :attr:`reader`
    reader: ScanIterator
        Any object implementing the :class:`ScanIterator` interface, produced by calling
        :attr:`loader_type` on :attr:`data_source`.
    ms1_deconvolution_args : dict
        The arguments passed to :func:`ms_deisotope.deconvolution.deconvolute_peaks` for MS^1
        scans.
    ms1_peak_picking_args : dict
        The arguments passed to :func:`ms_peak_picker.pick_peaks` for MS^1 scans.
    msn_deconvolution_args : dict
        The arguments passed to :func:`ms_deisotope.deconvolution.deconvolute_peaks` for MS^n
        scans.
    msn_peak_picking_args : dict
        The arguments passed to :func:`ms_peak_picker.pick_peaks` for MS^n scans.
    pick_only_tandem_envelopes : bool
        Whether or not to process whole MS^1 scans or just the regions around those peaks
        chosen for MS^n
    precursor_selection_window : float
        Size of the selection window to use when `pick_only_tandem_envelopes` is `True`
        and the information is not available in the scan.
    trust_charge_hint : bool
        Whether or not to trust the charge provided by the data source when determining
        the charge state of precursor isotopic patterns. Defaults to `True`
    """

    def __init__(self, data_source, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 pick_only_tandem_envelopes=False, precursor_selection_window=1.5,
                 trust_charge_hint=True,
                 loader_type=None,
                 envelope_selector=None):
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
        self.precursor_selection_window = precursor_selection_window
        self.trust_charge_hint = trust_charge_hint

        self.loader_type = loader_type

        self._signal_source = self.loader_type(data_source)
        self.envelope_selector = envelope_selector

    @property
    def reader(self):
        return self._signal_source

    def pick_precursor_scan_peaks(self, precursor_scan):
        logger.info("Picking Precursor Scan Peaks: %r", precursor_scan)
        if precursor_scan.is_profile:
            peak_mode = 'profile'
        else:
            peak_mode = 'centroid'
        prec_mz, prec_intensity = precursor_scan.arrays
        if not self.pick_only_tandem_envelopes and self.envelope_selector is None:
            prec_peaks = pick_peaks(prec_mz, prec_intensity, peak_mode=peak_mode, **self.ms1_peak_picking_args)
        else:
            if self.envelope_selector is None:
                chosen_envelopes = [s.precursor_information for s in precursor_scan.product_scans]
                chosen_envelopes = sorted([(p.mz - 5, p.mz + 10) for p in chosen_envelopes])
            else:
                chosen_envelopes = self.envelope_selector(precursor_scan)
            prec_peaks = pick_peaks(prec_mz, prec_intensity, peak_mode=peak_mode,
                                    target_envelopes=chosen_envelopes,
                                    **self.ms1_peak_picking_args)

        precursor_scan.peak_set = prec_peaks
        return prec_peaks

    def pick_product_scan_peaks(self, product_scan):
        if product_scan.is_profile:
            peak_mode = 'profile'
        else:
            peak_mode = 'centroid'
        product_mz, product_intensity = product_scan.arrays
        peaks = pick_peaks(product_mz, product_intensity, peak_mode=peak_mode, **self.msn_peak_picking_args)
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
            target = PriorityTarget(peak, precursor_ion, self.trust_charge_hint)
            if abs(err) > self.precursor_selection_window:
                logger.info(
                    "Unable to locate a peak for precursor ion %r for tandem scan %s of precursor scan %s" % (
                        precursor_ion, scan.title,
                        precursor_scan.title))
            else:
                priorities.append(target)
        return priorities

    def process_scan_group(self, precursor_scan, product_scans):
        """Performs the initial extraction of information relating
        `precursor_scan` to `product_scans` and picks peaks for `precursor_scan`.
        Called by :meth:`process`. May be used separately if doing the process step
        by step.

        Parameters
        ----------
        precursor_scan : Scan
            An MS^1 Scan
        product_scans : list of Scan
            A list of MS^n Scans related to `precursor_scan`

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

        for scan in product_scans:
            precursor_ion = scan.precursor_information
            peak, err = prec_peaks.get_nearest_peak(precursor_ion.mz)
            precursor_ion.peak = peak
            target = PriorityTarget(peak, precursor_ion, self.trust_charge_hint)
            if abs(err) > self.precursor_selection_window:
                logger.info(
                    "Unable to locate a peak for precursor ion %r for tandem scan %s of precursor scan %s" % (
                        precursor_ion, scan.id,
                        precursor_scan.id))
            else:
                priorities.append(target)

        return precursor_scan, priorities, product_scans

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

        decon_result = deconvolute_peaks(
            precursor_scan.peak_set, priority_list=priorities,
            **ms1_deconvolution_args)

        dec_peaks, priority_results = decon_result

        if decon_result.errors:
            logger.error("Errors occurred during deconvolution of %s, %r" % (
                precursor_scan.scan_id, decon_result.errors))

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

        for product_scan in precursor_scan.product_scans:
            precursor_information = product_scan.precursor_information

            i = get_nearest_index(precursor_information.mz, priorities)

            # If no peak is found in the priority list, it means the priority list is empty.
            # This should never happen in the current implementation. If it did, then we forgot
            # to pass the priority list to this function.
            if i is None:
                logger.info(
                    "Could not find deconvolution for %r (No nearby peak in the priority list)",
                    precursor_information)
                precursor_information.default()
                continue

            peak = priority_results[i]
            # If the deconvolution result is None, then we have no answer
            if peak is None:
                logger.info(
                    "Could not find deconvolution for %r (No solution was found for this region)",
                    precursor_information)
                precursor_information.default()

                continue
            elif peak.charge == 1 or (peak.charge != precursor_information.charge and self.trust_charge_hint):
                logger.info(
                    "Could not find deconvolution for %r (Unacceptable solution was proposed: %r)",
                    precursor_information, peak)
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
        if top_charge_state is not None and top_charge_state != 0 and abs(
                top_charge_state) < abs(charge_range[1]):
            charge_range[1] = top_charge_state

        deconargs["charge_range"] = charge_range

        if product_scan.polarity in (-1, 1):
            polarity = product_scan.polarity
            deconargs["charge_range"] = [
                polarity * abs(c) for c in deconargs["charge_range"]]

        dec_peaks, _ = deconvolute_peaks(product_scan.peak_set, **deconargs)
        product_scan.deconvoluted_peak_set = dec_peaks
        return dec_peaks

    def _get_next_scans(self):
        precursor, products = next(self.reader)

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
            An MS^1 Scan
        products : list of Scan
            A list of MS^n Scans related to `precursor`

        Returns
        -------
        precursor_scan: Scan
            The fully processed version of `precursor`
        product_scans: list of Scan
            The fully processed version of `products`
        """
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        self.deconvolute_precursor_scan(precursor_scan, priorities)

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
        return ScanBunch(precursor_scan.pack(), [p.pack() for p in product_scans])

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


class PrecursorSkimmingProcessor(ScanProcessor):
    def process(self, precursor, products):
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        self.deconvolute_precursor_scan(precursor_scan, priorities)
        return precursor_scan, product_scans


def extract_intervals(scan_iterator, time_radius=5., mz_lower=2., mz_higher=3.):
    intervals = []
    for scan, products in scan_iterator:
        for product in products:
            intervals.append((
                Interval(max(0, product.precursor_information.mz - mz_lower),
                         product.precursor_information.mz + mz_higher),
                Interval(max(0, product.scan_time - time_radius),
                         product.scan_time + time_radius)))
    return intervals


def combine_intervals(a, b):
    min_mz = min(a[0].start, b[0].start)
    max_mz = max(a[0].end, b[0].end)
    min_rt = min(a[1].start, b[1].start)
    max_rt = max(a[1].end, b[1].end)
    return (Interval(min_mz, max_mz), Interval(min_rt, max_rt))


def overlaps_2d(a, b):
    return a[0].overlaps(b[0]) and a[1].overlaps(b[1])


def merge_interval_set(intervals, minimum_overlap_size=0.3):
    merged_intervals = []
    for interval in intervals:
        for i, candidate in enumerate(merged_intervals):
            if overlaps_2d(interval, candidate) and interval[0].overlap_size(
                    candidate[0]) > (interval[0].end - interval[0].start * minimum_overlap_size):
                merged_intervals[i] = combine_intervals(candidate, interval)
                break
        else:
            merged_intervals.append(interval)
    return merged_intervals


def make_rt_tree(intervals):
    temp = []
    for interval in intervals:
        mz, rt = interval
        rt.members.append(mz)
        temp.append(rt)
    tree = IntervalTreeNode.build(temp)
    return tree


class ScanIntervalTree(object):
    @classmethod
    def build(cls, scan_iterator, time_radius=5., mz_lower=2., mz_higher=3.):
        intervals = extract_intervals(scan_iterator, time_radius=time_radius, mz_lower=mz_lower, mz_higher=mz_higher)
        # merged_intervals = merge_interval_set(intervals)
        merged_intervals = intervals
        return cls(make_rt_tree(merged_intervals), intervals)

    def __init__(self, rt_tree, original_intervals=None):
        self.rt_tree = rt_tree
        self.original_intervals = original_intervals

    def get_mz_intervals_for_rt(self, rt_point):
        intervals = [
            i.members[0] for i in self.rt_tree.contains_point(rt_point)
        ]
        intervals = sorted([[i.start, i.end] for i in intervals])
        out = []
        if len(intervals) == 0:
            return []
        last = intervals[0]
        for interval in intervals[1:]:
            if interval[0] < last[1]:
                last[1] = interval[1]
            else:
                out.append(last)
                last = interval
        out.append(last)
        return out

    def __call__(self, scan):
        intervals = self.get_mz_intervals_for_rt(scan.scan_time)
        return intervals
