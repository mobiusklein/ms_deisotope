'''
High Level Deconvolution API
----------------------------

A high-level wrapper around the deconvolution machinery, orchestrating the
process of constructing a deconvoluter instance, performing deconvolution, and
extracting targets of interest from the result.
'''

from ms_peak_picker import FittedPeak

from ms_deisotope.peak_set import merge
from ms_deisotope.utils import DeconvolutionProcessResult
from ms_deisotope.constants import (
    ERROR_TOLERANCE,
    TRUNCATE_AFTER,
    MAX_ITERATION,
    SCALE_METHOD)

from ms_deisotope.averagine import PROTON

from .averagine_based import AveraginePeakDependenceGraphDeconvoluter
from .utils import logger, prepare_peaklist


def deconvolute_peaks(peaklist, decon_config=None,
                      charge_range=(1, 8), error_tolerance=ERROR_TOLERANCE,
                      priority_list=None, left_search_limit=1, right_search_limit=0,
                      left_search_limit_for_priorities=None, right_search_limit_for_priorities=None,
                      verbose_priorities=False, verbose=False, charge_carrier=PROTON,
                      truncate_after=TRUNCATE_AFTER, iterations=MAX_ITERATION,
                      deconvoluter_type=AveraginePeakDependenceGraphDeconvoluter,
                      retention_strategy=None,
                      use_quick_charge=False,
                      **kwargs):
    """Deconvolute a centroided mass spectrum

    This function constructs a deconvoluter object using the ``deconvoluter_type`` argument
    and deconvolutes the input ``peaklist`` by calling its :meth:`deconvolute` method.

    If ``priority_list`` is not :const:`None`, it is expected to be an iterable of either
    tuples of (:class:`~.FittedPeak`, ``(min charge, max charge)``) pairs, or instances of
    :class:`~.PriorityTarget`. These will be passed to :meth:`targeted_deconvolution` of
    the deconvoluter.

    Parameters
    ----------
    peaklist : :class:`~.PeakSet` or list of Peak-like objects, see :func:`~.prepare_peaklist`
        The centroided mass spectrum to deconvolute.
    decon_config : dict, optional
        Parameters to use to initialize the deconvoluter instance produced by
        ``deconvoluter_type``
    charge_range : tuple of integers, optional
        The range of charge states to consider. The range is inclusive.
    error_tolerance : float, optional
        PPM error tolerance to use to match experimental to theoretical peaks
    priority_list : list, optional
        The set of peaks to target for deconvolution to be able to enforce external
        constraints on, such as selected precursors for fragmentation.
    left_search_limit : int, optional
        The maximum number of neutron shifts to search to the left  (decrease) from
        each query peak
    right_search_limit : int, optional
        The maximum number of neutron shifts to search to the right (increase) from
        each query peak
    left_search_limit_for_priorities : int, optional
        The maximum number of neutron shifts to search to the left (decrease) from
        each query peak for priority targets
    right_search_limit_for_priorities : int, optional
        The maximum number of neutron shifts to search to the right (increase) from
        each query peak for priority targets
    verbose_priorities : bool, optional
        Whether to turn on verbose mode for priority targets
    verbose : bool, optional
        Passed to the deconvoluter to enable verbose mode globally
    charge_carrier : float, optional
        The mass of the charge carrier. Defaults to |PROTON|
    truncate_after : float, optional
        The percentage of the isotopic pattern to include. Defaults to |TRUNCATE_AFTER|
    deconvoluter_type : type or callable, optional
        A callable returning a deconvoluter. Defaults to :class:`~.AveraginePeakDependenceGraphDeconvoluter`
    retention_strategy: :class:`~.PeakRetentionStrategyBase` or callable, optional
        A callable that may compute additional peaks to include in the output deconvoluted peaks.
    use_quick_charge: :class:`bool`
        Whether or not to use the :title-reference:`QuickCharge` algorithm to quickly filter
        theoretical charge states to consider for each peak.
    **kwargs
        Additional keywords included in ``decon_config``


    Notes
    -----
    If speed is an issue, consider setting `use_quick_charge` :const:`True`. This will pre-screen charge
    states that are are missing peaks supporting the :term:`monoisotopic peak` or the :term:`A+1 peak`.
    Alternatively, you may set the charge range upper bound to something reasonable for your data, such
    as the precursor ion's charge when considering a product ion spectrum.

    When you do not expect a complete isotopic pattern for large ions, as is often the case for low
    abundance FT-MS/MS it may be useful to shrink the `truncate_after` parameter from the default (|TRUNCATE_AFTER|)
    to a slightly smaller value.

    Returns
    -------
    :class:`~.DeconvolutionProcessResult`
    """
    if priority_list is None:
        priority_list = []
    if left_search_limit_for_priorities is None:
        left_search_limit_for_priorities = left_search_limit
    if right_search_limit_for_priorities is None:
        right_search_limit_for_priorities = right_search_limit

    decon_config = decon_config or {}
    decon_config.update(kwargs)
    decon_config.setdefault("use_subtraction", True)
    decon_config.setdefault("scale_method", SCALE_METHOD)
    decon_config.setdefault("use_quick_charge", use_quick_charge)
    decon = deconvoluter_type(peaklist=peaklist, **decon_config)

    peaklist = prepare_peaklist(peaklist)

    if verbose_priorities or verbose:
        decon.verbose = True

    priority_list_results = []
    for p in priority_list:
        try:
            target_info = p
            p = target_info.peak
            hinted_charge_range = target_info.charge_range_hint(charge_range)
        except AttributeError:
            hinted_charge_range = charge_range
        if not isinstance(p, FittedPeak):
            p = decon.peaklist.has_peak(p, error_tolerance)
        priority_result = decon.targeted_deconvolution(
            p, error_tolerance=error_tolerance,
            charge_range=hinted_charge_range,
            left_search_limit=left_search_limit_for_priorities,
            right_search_limit=right_search_limit_for_priorities,
            charge_carrier=charge_carrier,
            truncate_after=truncate_after)
        priority_list_results.append(priority_result)

    if verbose_priorities and not verbose:
        decon.verbose = False

    deconvoluted_peaks = decon.deconvolute(
        error_tolerance=error_tolerance, charge_range=charge_range,
        left_search_limit=left_search_limit, right_search_limit=right_search_limit,
        charge_carrier=charge_carrier, truncate_after=truncate_after,
        iterations=iterations)

    acc = []
    errors = []
    for pr in priority_list_results:
        try:
            result = pr.get()
        except Exception as e:
            result = None
            errors.append(e)
            logger.error("Could not extract a solution for %r",
                         pr.query_peak, exc_info=True)
        acc.append(result)

    priority_list_results = acc

    if retention_strategy is not None:
        retained = retention_strategy(
            decon.peaklist, peaklist, charge_range, deconvoluted_peaks)
        deconvoluted_peaks = merge(deconvoluted_peaks, retained)

    return DeconvolutionProcessResult(
        decon, deconvoluted_peaks, priority_list_results, errors)
