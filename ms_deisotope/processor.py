'''
Deconvolution Pipeline
----------------------

The deconvolution process can be pipelined from start to finish
using the :class:`~.ScanProcessor` class. This includes precursor recalculation and
coisolation detection. :class:`~.ScanProcessor` is intended to be used either as a
replacement for :class:`~.ScanIterator` when deconvolution is desired, or that it's
:meth:`ScanProcessor.process` method will be used to handle individual :class:`~.ScanBunch`
objects/pairs of precursor :class:`~.Scan` objects and a :class:`list` of product :class:`~.Scan`
objects.

The free-function :func:`~.process` is a convenience wrapper around :class:`~.ScanProcessor`,
with fewer configurable parameters.

.. code:: python

    from ms_deisotope import ScanProcessor, glycopeptide, peptide
    from ms_deisotope.scoring import PenalizedMSDeconVFitter, MSDeconVFitter
    from ms_deisotope.test.common import datafile

    # Locate example dataset
    path = datafile("20150710_3um_AGP_001_29_30.mzML.gz")

    proc = processor.ScanProcessor(path, ms1_deconvolution_args={
        "averagine": glycopeptide,
        "scorer": PenalizedMSDeconVFitter(20., 2.),
        "truncate_after": 0.95
    }, msn_deconvolution_args={
        "averagine": peptide,
        "scorer": MSDeconVFitter(10.),
        "truncate_after": 0.8
    })

    bunch = next(proc)
    print(bunch)
    print(bunch.precursor.deconvoluted_peak_set)
'''
import logging

from six import string_types as basestring

import numpy as np

from ms_peak_picker import pick_peaks, PeakSet, PeakIndex
from ms_peak_picker.scan_filter import FTICRBaselineRemoval

from ms_deisotope import constants
from .averagine import AveragineCache, peptide, PROTON
from .scoring import PenalizedMSDeconVFitter, MSDeconVFitter
from .deconvolution import deconvolute_peaks
from .data_source import MSFileLoader, ScanIterator
from .data_source.common import Scan, ScanBunch, ChargeNotProvided
from .utils import Base
from .peak_dependency_network import NoIsotopicClustersError
from .qc.isolation import PrecursorPurityEstimator
from .task import LogUtilsMixin

logger = logging.getLogger("deconvolution_scan_processor")
logger.addHandler(logging.NullHandler())


def _get_nearest_index(query_mz, peak_list):
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
    info : :class:`~.PrecursorInformation`
        The associated precursor information block which contains
        the charge state hint.
    peak : :class:`~.FittedPeak`
        The peak from which to start the deconvolution
    trust_charge_hint : bool
        Whether or not to force the deconvoluter to only consider
        the charge specified in the hint.
    mz : float
        The m/z of :attr:`peak`
    charge : int
        The charge state hint from :attr:`info`
    isolation_window : :class:`~.IsolationWindow`
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
        '''
        The m/z of the matched peak

        Returns
        -------
        float
        '''
        return self.peak.mz

    @property
    def charge(self):
        '''
        The charge state of the precursor ion reported by the source.

        Returns
        -------
        int
        '''
        try:
            return int(self.info.charge)
        except TypeError:
            return 0

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


def _loader_creator(specification):
    if isinstance(specification, basestring):
        return MSFileLoader(specification)
    elif isinstance(specification, ScanIterator):
        return specification
    else:
        raise ValueError("Cannot determine how to get a ScanIterator from %r" % (specification,))


def _simplify_peak_set(peaks, bin_width=5.0):
    bin_edges = np.arange(0, peaks[-1].mz + bin_width, bin_width)
    bins = []
    for i, bin_edge in enumerate(bin_edges, 1):
        if i == len(bin_edges):
            next_edge = bin_edges[-1] + bin_width
        else:
            next_edge = bin_edges[i]
        subset = peaks.between(bin_edge, next_edge)
        bins.append(subset)

    thresholds = []
    reduced_subsets = {}
    k = 0
    for b in bins:
        if len(b) > 0:
            bin_intensities = np.array([p.intensity for p in b])
            thresholds.append(np.max(bin_intensities) / 3.)
            for p in b:
                if p.intensity > thresholds[-1]:
                    reduced_subsets[p.peak_count] = p
            k += (bin_intensities > thresholds[-1]).sum()
        else:
            thresholds.append(0.0)
    subset_peaks = PeakSet(
        sorted(reduced_subsets.values(), key=lambda x: x.mz)).clone()
    subset_peaks.reindex()
    return PeakIndex(np.array([]), np.array([]), subset_peaks)


class ScanProcessor(Base, LogUtilsMixin):
    """Orchestrates the deconvolution of a :class:`~.ScanIterator` scan by scan. This process will
    apply different rules for MS1 scans and MSn scans. This type itself mimics a :class:`~.ScanIterator`,
    consuming (raw) mass spectral data and producing deisotoped and charge deconvolved spectra.

    The algorithms used for each task are independent and can be specified in the appropriate
    attribute dictionary, however there is information sharing between each MS1 scan and its
    MSn scans as the precursor monoisotopic mass is recalibrated according to the MS1 processing
    arguments, and the selected charge state is used to limit the charge range used in the matching
    MSn scan. These are described by :class:`PriorityTarget` objects.

    If an averagine-based deconvoluter is used, the averagine cache will be pre-populated.

    At the moment, MSn assumes only MS2. Until MS3 data become available for testing, this limit
    will remain.

    Attributes
    ----------
    data_source : :class:`str`, :class:`~.ScanIterator` or file-like
        Any valid object to be passed to the `loader_type` callable to produce
        a :class:`~.ScanIterator` instance. A path to a mass spectrometry data file,
        a file-like object, or an instance of :class:`~.ScanIterator`. Used to populate :attr:`reader`
    loader_type : callable
        A callable, which when passed :attr:`data_source` returns an instance of :class:`~.ScanIterator`.
        By default, this is :func:`~.MSFileLoader`. Used to populate :attr:`reader`
    reader: ScanIterator
        Any object implementing the :class:`~.ScanIterator` interface, produced by calling
        :attr:`loader_type` on :attr:`data_source`.
    ms1_deconvolution_args : :class:`dict`
        The arguments passed to :func:`~ms_deisotope.deconvolution.deconvolute_peaks` for MS1
        scans.
    ms1_peak_picking_args : :class:`dict`
        The arguments passed to :func:`ms_peak_picker.pick_peaks` for MS1 scans.
    msn_deconvolution_args : :class:`dict`
        The arguments passed to :func:`~ms_deisotope.deconvolution.deconvolute_peaks` for MSn
        scans.
    msn_peak_picking_args : :class:`dict`
        The arguments passed to :func:`ms_peak_picker.pick_peaks` for MSn scans.
    pick_only_tandem_envelopes : :class:`bool`
        Whether or not to process whole MS1 scans or just the regions around those peaks
        chosen for MSn
    default_precursor_ion_selection_window : :class:`float`
        Size of the selection window to use when :attr:`pick_only_tandem_envelopes` is `True`
        and the information is not available in the scan.
    trust_charge_hint : :class:`bool`
        Whether or not to trust the charge provided by the data source when determining
        the charge state of precursor isotopic patterns. Defaults to `True`
    respect_isolation_window: :class:`bool`
        Whether to use the bounds of the isolation window to reject a monoisotopic peak
        solution
    terminate_on_error: :class:`bool`
        Whether or not  to stop processing on an error. Defaults to `True`
    ms1_averaging: :class:`int`
        The number of adjacent MS1 scans to average prior to picking peaks.
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
                 respect_isolation_window=False,
                 too_many_peaks_threshold=7000):
        if loader_type is None:
            loader_type = _loader_creator

        self.data_source = data_source
        self.ms1_peak_picking_args = ms1_peak_picking_args or {}
        self.msn_peak_picking_args = msn_peak_picking_args or ms1_peak_picking_args or {}
        self.ms1_deconvolution_args = ms1_deconvolution_args or {}
        self.ms1_deconvolution_args.setdefault("charge_range", (1, 8))
        self.msn_deconvolution_args = msn_deconvolution_args or {}
        self.msn_deconvolution_args.setdefault("charge_range", (1, 8))
        self.pick_only_tandem_envelopes = pick_only_tandem_envelopes

        self.too_many_peaks_threshold = too_many_peaks_threshold

        self.default_precursor_ion_selection_window = default_precursor_ion_selection_window
        self.respect_isolation_window = respect_isolation_window

        self.trust_charge_hint = trust_charge_hint
        self.ms1_averaging = int(ms1_averaging) if ms1_averaging else 0

        self.loader_type = loader_type

        self._signal_source = self.loader_type(data_source)
        self.envelope_selector = envelope_selector
        self.terminate_on_error = terminate_on_error
        self._prepopulate_averagine_cache()

    def _prepopulate_averagine_cache(self):
        if 'averagine' in self.ms1_deconvolution_args:
            averagine = self.ms1_deconvolution_args['averagine']
            ms1_truncate_after = self.ms1_deconvolution_args.get(
                'truncate_after', constants.TRUNCATE_AFTER)
            ms1_ignore_below = self.ms1_deconvolution_args.get(
                'ignore_below', constants.IGNORE_BELOW)
            ms1_charge_range = self.ms1_deconvolution_args.get('charge_range', (1, 8))
            ms1_charge_carrier = self.ms1_deconvolution_args.get(
                'charge_carrier', PROTON)
            if isinstance(averagine, (list, tuple)):
                averagine = [
                    AveragineCache(a).populate(
                        truncate_after=ms1_truncate_after,
                        ignore_below=ms1_ignore_below,
                        min_charge=ms1_charge_range[0],
                        max_charge=ms1_charge_range[1],
                        charge_carrier=ms1_charge_carrier)
                    for a in averagine]
            else:
                averagine = AveragineCache(averagine).populate(
                    truncate_after=ms1_truncate_after,
                    ignore_below=ms1_ignore_below,
                    min_charge=ms1_charge_range[0],
                    max_charge=ms1_charge_range[1],
                    charge_carrier=ms1_charge_carrier)
            self.ms1_deconvolution_args['averagine'] = averagine
        if 'averagine' in self.msn_deconvolution_args:
            averagine = self.msn_deconvolution_args['averagine']
            msn_truncate_after = self.msn_deconvolution_args.get(
                'truncate_after', constants.TRUNCATE_AFTER)
            msn_ignore_below = self.msn_deconvolution_args.get(
                'ignore_below', constants.IGNORE_BELOW)
            msn_charge_range = self.msn_deconvolution_args.get(
                'charge_range', (1, 8))
            msn_charge_carrier = self.msn_deconvolution_args.get(
                'charge_carrier', PROTON)
            if isinstance(averagine, (list, tuple)):
                averagine = [
                    AveragineCache(a).populate(
                        truncate_after=msn_truncate_after,
                        ignore_below=msn_ignore_below,
                        min_charge=msn_charge_range[0],
                        max_charge=msn_charge_range[1],
                        charge_carrier=msn_charge_carrier
                    ) for a in averagine]
            else:
                averagine = AveragineCache(averagine).populate(
                    truncate_after=msn_truncate_after,
                    ignore_below=msn_ignore_below,
                    min_charge=msn_charge_range[0],
                    max_charge=msn_charge_range[1],
                    charge_carrier=msn_charge_carrier)
            self.msn_deconvolution_args['averagine'] = averagine

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
        '''The :class:`~.ScanIterator` which generates the raw scans that will
        be processed.

        Returns
        -------
        :class:`~.ScanIterator`
        '''
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
        # averaged scans are always profile mode
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
        self.log("Picking Precursor Scan Peaks: %r" % (precursor_scan, ))
        if self.ms1_averaging > 0:
            prec_peaks = self._average_ms1(precursor_scan)
        else:
            prec_peaks = self._pick_precursor_scan_peaks(precursor_scan)
        n_peaks = len(prec_peaks)
        if n_peaks > self.too_many_peaks_threshold:
            self.log("%d peaks found for %r, applying local intensity threshold." % (n_peaks, precursor_scan))
            prec_peaks = _simplify_peak_set(prec_peaks)
        precursor_scan.peak_set = prec_peaks
        return prec_peaks

    def pick_product_scan_peaks(self, product_scan):
        """Pick the peaks of product scan

        Parameters
        ----------
        product_scan: :class:`~.Scan`
            The scan to pick peaks from.

        Returns
        -------
        PeakSet
        """
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
        """A utility method to obtain :class:`PriorityTarget` objects for
        each product scan of `precursor_scan`.

        Parameters
        ----------
        precursor_scan: :class:`~.Scan`
            The scan to extract.

        Returns
        -------
        :class:`list` of :class:`PriorityTarget`
        """
        priorities = []
        peaks = precursor_scan.peak_set
        for scan in precursor_scan.product_scans:
            precursor_ion = scan.precursor_information
            if peaks is None:
                peaks = self.pick_precursor_scan_peaks(precursor_scan)
            peak, _ = peaks.get_nearest_peak(precursor_ion.mz)
            precursor_ion.peak = peak
            target = PriorityTarget(
                peak,
                precursor_ion,
                self.trust_charge_hint,
                isolation_window=scan.isolation_window)
            if self._reject_candidate_precursor_peak(peak, scan):
                self.log(
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
        precursor_scan : :class:`~.Scan`
            An MS1 Scan
        product_scans : :class:`list` of :class:`~.Scan`
            A :class:`list` of MSn Scans related to `precursor_scan`

        Returns
        -------
        precursor_scan: :class:`~.Scan`
            As Parameter
        prioritiies: :class:`list` of :class:`~PriorityTarget`
            :class:`list` of the peak target windows in `precursor_scan` which
            are related to `product_scans`
        product_scans: :class:`list` of :class:`~.Scan`
            As Parameter
        """
        prec_peaks = self.pick_precursor_scan_peaks(precursor_scan)
        priorities = []

        if prec_peaks is None:
            raise EmptyScanError(
                "Could not pick peaks for empty precursor scan", self)

        for scan in product_scans:
            precursor_ion = scan.precursor_information
            if precursor_ion is None:
                continue
            peak = prec_peaks.has_peak(precursor_ion.mz)
            if peak is not None:
                err = abs(peak.mz - precursor_ion.mz)
            else:
                peak, err = prec_peaks.get_nearest_peak(precursor_ion.mz)
            self.debug("For Precursor at %0.4f, found Peak at %0.4f with error %0.4f" % (
                precursor_ion.mz, peak.mz, err))
            precursor_ion.peak = peak
            target = PriorityTarget(
                peak, precursor_ion, self.trust_charge_hint,
                scan.precursor_information.precursor_scan_id,
                scan.precursor_information.product_scan_id,
                isolation_window=scan.isolation_window)
            if self._reject_candidate_precursor_peak(peak, scan):
                self.log(
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

    def deconvolute_precursor_scan(self, precursor_scan, priorities=None, product_scans=None):
        """Deconvolute the given precursor scan, giving priority to its product ions,
        correcting the :attr:`precursor_information` attributes of priority targets,
        as well as calculating the degree of precursor purity and coisolating ions.

        Parameters
        ----------
        precursor_scan : :class:`~.Scan`
            The precursor scan to deconvolute
        priorities : :class:`list` of :class:`PriorityTarget`, optional
            The priority targets for the product ions derived from `precursor_scan`
        product_scans: :class:`list` of :class:`~.Scan`
            The product ion scans of `precursor_scan`.

        Returns
        -------
        :class:`~DeconvolutedPeakSet`
            The deconvoluted peaks of ``precursor_scan``
        :class:`list` of :class:`PriorityTarget`
            The precursor ions selected, with updated mass and charge information

        Raises
        ------
        Exception
            Any errors which are thrown during the deconvolution process may be thrown
            if :attr:`terminate_on_error` is `True`.
        """
        if priorities is None:
            priorities = []
        if product_scans is None:
            product_scans = []

        self.log("Deconvoluting Precursor Scan %r" % precursor_scan)
        self.log("Priorities: %r" % priorities)

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
                self.log("No isotopic clusters found in %r" % precursor_scan.id)

        dec_peaks, priority_results = decon_result
        if decon_result.errors:
            self.error("Errors occurred during deconvolution of %s, %r" % (
                precursor_scan.id, decon_result.errors))
        precursor_scan.deconvoluted_peak_set = dec_peaks
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

        coisolation_detection = PrecursorPurityEstimator(default_width=self.default_precursor_ion_selection_window)
        self.debug("Priority Targets for %s: %r" % (
            precursor_scan.id, [
                (p.mz, p.charge) if p is not None else None for p in priorities
            ]))
        if priorities:
            if not product_scans:
                self.debug("Priority targets were passed without product scans")
        for product_scan in product_scans:
            precursor_information = product_scan.precursor_information
            if precursor_information is None:
                continue
            # unknown precursor purity
            product_scan.annotations['precursor purity'] = 0.0
            i = _get_nearest_index(precursor_information.mz, priorities)

            # If no peak is found in the priority list, it means the priority list is empty.
            # This should never happen in the current implementation. If it did, then we forgot
            # to pass the priority list to this function.
            if i is None:
                self.log(
                    "Could not find deconvolution for %r (No nearby peak in the priority list)" %
                    precursor_information)
                precursor_information.default(orphan=True)
                continue

            peak = priority_results[i]
            # If the deconvolution result is None, then we have no answer
            if peak is None:
                self.log(
                    "Could not find deconvolution for %r (No solution was found for this region)" %
                    precursor_information)
                precursor_information.default(orphan=True)
                coisolation = coisolation_detection.coisolation(
                    precursor_scan, None, product_scan.isolation_window, 0.0)
                precursor_information.coisolation = coisolation
                continue
            elif peak.charge == 1 or (peak.charge != precursor_information.charge and self.trust_charge_hint):
                if precursor_information.charge != ChargeNotProvided:
                    self.log(
                        "Could not find deconvolution for %r (Unacceptable solution was proposed: %r)" %
                        (precursor_information, peak))
                    precursor_information.default()
                    continue

            precursor_purity = -1.0
            if peak is not None:
                precursor_purity, coisolation = coisolation_detection(
                    precursor_scan,
                    peak,
                    product_scan.isolation_window)
                precursor_information.coisolation = coisolation
                self.debug(
                    "Precursor m/z %f\nExperimental = %r\nTheoretical = %r" % (
                        peak.mz,
                        ', '.join(["(%0.4f, %0.1f)" % (p.mz, p.intensity) for p in peak.envelope]),
                        ', '.join(["(%0.4f, %0.1f)" % (p.mz, p.intensity) for p in peak.fit.theoretical]))
                )
            else:
                coisolation = coisolation_detection.coisolation(
                    precursor_scan, None, product_scan.isolation_window, 0.0)
                precursor_information.coisolation = coisolation

            product_scan.annotations['precursor purity'] = precursor_purity
            precursor_information.extract(peak)
        return dec_peaks, priority_results

    def deconvolute_product_scan(self, product_scan):
        """Deconvolute the peaks of `product_scan`.

        This method will override the upper limit "charge_range" of
        :attr:`msn_deconvolution_args` to the charge information of
        the precursor ion.

        This method sets the :attr:`~.Scan.deconvoluted_peak_set` of
        `product_scan`.

        Parameters
        ----------
        product_scan : :class:`~.Scan`
            The scan to deconvolute.

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`

        Raises
        ------
        Exception
            Any errors which are thrown during the deconvolution process may be thrown
            if :attr:`terminate_on_error` is `True`.
        """
        self.log("Deconvoluting Product Scan %r" % (product_scan, ))
        precursor_ion = product_scan.precursor_information
        deconargs = dict(self.msn_deconvolution_args)
        if precursor_ion is not None:
            top_charge_state = precursor_ion.extracted_charge
            if not top_charge_state:
                top_charge_state = precursor_ion.charge
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
            self.log("No Isotopic Clusters found in %r" % product_scan.id)
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

        This method may be used to process scans from other sources not from the
        wrapped :class:`~.ScanIterator`.

        Parameters
        ----------
        precursor : :class:`~.Scan`
            An MS1 Scan
        products : :class:`list` of :class:`~.Scan`
            A list of MSn Scans related to `precursor`

        Returns
        -------
        :class:`~.ScanBunch`
            The fully processed version of `precursor` and `products`
        """
        precursor_scan, priorities, product_scans = self.process_scan_group(precursor, products)
        if precursor_scan is not None:
            self.deconvolute_precursor_scan(
                precursor_scan, priorities, product_scans)
        else:
            self._default_all_precursor_information(product_scans)

        for product_scan in product_scans:
            self.pick_product_scan_peaks(product_scan)
            self.deconvolute_product_scan(product_scan)

        return ScanBunch(precursor_scan, product_scans)

    def next(self):
        """Fetches the next bunch of scans from :attr:`reader` and
        invokes :meth:`process` on them, picking peaks and deconvoluting them.

        Returns
        -------
        ScanBunch
        """
        precursor, products = self._get_next_scans()
        bunch = self.process(precursor, products)
        return bunch

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
        """A wrapper around :meth:`~.RandomAccessScanSource.start_from_scan` provided by
        :attr:`reader`, if available.

        Returns
        -------
        self

        See Also
        --------
        :meth:`~.RandomAccessScanSource.start_from_scan`
        """
        self.reader.start_from_scan(*args, **kwargs)
        return self


ScanProcessor.log_with_logger(logger)


def process(data_source, ms1_averagine=peptide, msn_averagine=peptide,
            ms1_score_threshold=20, msn_score_threshold=5, denoise=False,
            ms1_max_missed_peaks=1, pick_only_tandem_envelopes=False,
            trust_charge_hint=True, envelope_selector=None, terminate_on_error=True,
            ms1_averaging=0, respect_isolation_window=False, use_quick_charge=True):
    """Construct a deconvolution pipeline for common applications.

    Parameters
    ----------
    data_source : :class:`str` or :class:`~.ScanIterator` or file-like object
        The scan data source to read raw spectra from
    ms1_averagine : :class:`~.Averagine` or :class:`~.AveragineCache`, optional
        The :class:`~.Averagine` model to use for MS1 scans. Defaults to
        :data:`ms_deisotope.averagine.peptide`.
    msn_averagine : :class:`~.Averagine` or :class:`~.AveragineCache`, optional
        The :class:`~.Averagine` model to use for MSn scans. Defaults to
        :data:`ms_deisotope.averagine.peptide`.
    ms1_score_threshold : float, optional
        The score threshold to use to reject isotopic pattern fits for MS1 scans.
        The default is 20.0.
    msn_score_threshold : float, optional
        The score threshold to use to reject isotopic pattern fits for MS1 scans.
        The default is 5.0.
    denoise : :class:`bool` or :class:`float`, optional
        Whether to denoise MS1 scans. If the value is not false-y, it may either be
        a float to set the scale of the denoising process, or 5.0 if the value is
        :const:`True`.
    ms1_max_missed_peaks : :class:`int`, optional
        The maximum number of missed peaks to permit for MS1 scans. The default is 1.
    pick_only_tandem_envelopes : :class:`bool`
        Whether or not to process whole MS1 scans or just the regions around those peaks
        chosen for MSn
    default_precursor_ion_selection_window : :class:`float`
        Size of the selection window to use when an explicit isolation window width is not
        available in the scan.
    trust_charge_hint : :class:`bool`
        Whether or not to trust the charge provided by the data source when determining
        the charge state of precursor isotopic patterns. Defaults to `True`
    terminate_on_error: :class:`bool`
        Whether or not  to stop processing on an error. Defaults to `True`
    ms1_averaging: :class:`int`
        The number of adjacent MS1 scans to average prior to picking peaks.
    respect_isolation_window: :class:`bool`
        Whether to use the bounds of the isolation window to reject a monoisotopic peak
        solution
    use_quick_charge : :class:`bool`, optional
        Whether or not to used the QuickCharge algorithm for expediting charge calculation.

    Returns
    -------
    :class:`ScanProcessor`
    """
    if denoise:
        ms1_peak_picking_args = {
            "filters": [
                FTICRBaselineRemoval(
                    scale=denoise if denoise is not True else 5., window_length=2)
            ]
        }
    else:
        ms1_peak_picking_args = None

    ms1_deconvolution_args = {
        "averagine": ms1_averagine,
        "scorer": PenalizedMSDeconVFitter(ms1_score_threshold, 2.0),
        "use_quick_charge": use_quick_charge,
        "max_missed_peaks": ms1_max_missed_peaks,
        "truncate_after": 0.95,
    }
    msn_deconvolution_args = {
        "averagine": msn_averagine,
        "scorer": MSDeconVFitter(msn_score_threshold),
        "use_quick_charge": use_quick_charge,
        "truncate_after": 0.8,
    }
    processor = ScanProcessor(
        data_source, ms1_peak_picking_args,
        None, ms1_deconvolution_args,
        msn_deconvolution_args, pick_only_tandem_envelopes,
        trust_charge_hint=trust_charge_hint,
        envelope_selector=envelope_selector,
        terminate_on_error=terminate_on_error,
        ms1_averaging=ms1_averaging,
        respect_isolation_window=respect_isolation_window)
    return processor


class EmptyScanError(ValueError):
    """A sub-type of :class:`ValueError` which is used to indicate
    that a spectrum is empty and could not be manipulated.
    """
    def __init__(self, msg, scan_id=None):
        ValueError.__init__(self, msg)
        self.scan_id = scan_id
