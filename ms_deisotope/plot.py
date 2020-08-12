'''A collection of tools for drawing and annotating mass spectra
'''
# pragma: no cover
import math
import itertools

import numpy as np
try:
    from matplotlib import pyplot as plt, gridspec
    from ms_peak_picker.plot import draw_peaklist, draw_raw
    has_plot = True
except ImportError as err:
    import warnings
    import traceback
    warnings.warn("Could not import matplotlib, plotting tools will not work\n%s" % (err, ))
    pyplot = None
    gridspec = None
    has_plot = False


def _default_color_cycle():
    c = plt.rcParams['axes.prop_cycle']
    colors = c.by_key().get("color")
    if not colors:
        colors = [
            '#1f77b4',
            '#ff7f0e',
            '#2ca02c',
            '#d62728',
            '#9467bd',
            '#8c564b',
            '#e377c2',
            '#7f7f7f',
            '#bcbd22',
            '#17becf'
        ]
    return colors


def annotate_scan(scan, products, nperrow=4, ax=None, label=True):
    '''Given an MS1 :class:`~.ScanBase` ``scan`` and a :class:`~.Sequence` of
    :class:`~.ScanBase` product scans, draw the MS1 spectrum in full profile,
    and then in a subplot grid below it, draw a zoomed-in view of the MS1 spectrum
    surrounding the area around each precursor ion that gave rise to the scans in
    ``products``, with monoisotopic peaks and isolation windows marked.

    .. plot::
        :include-source:

        import ms_deisotope
        from ms_deisotope import plot
        from ms_deisotope.test.common import datafile

        reader = ms_deisotope.MSFileLoader(datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
        bunch = next(reader)

        bunch.precursor.pick_peaks()
        bunch.precursor.deconvolute(
            scorer=ms_deisotope.PenalizedMSDeconVFitter(20., 2.0),
            averagine=ms_deisotope.glycopeptide, use_quick_charge=True)

        ax = plot.annotate_scan(bunch.precursor, bunch.products, nperrow=2)
        ax.figure.set_figwidth(12)
        ax.figure.set_figheight(16)


    Parameters
    ----------
    scan: ScanBase
        The precursor MS1 scan
    products: :class:`~.Sequence` of :class:`~.ScanBase`
        The collection of MSn scans based upon ``scan``
    nperrow: int
        The number of precursor ion subplots to draw per row
        of the grid. Defaults to :const:`4`.
    ax: :class:`matplotlib._axes.Axes`
        An :class:`~.Axes` object to use to find the figure to
        draw the plot on.

    Returns
    -------
    :class:`matplotlib._axes.Axes`:
        The axes of the full MS1 profile plot
    '''
    if ax is None:
        figure = plt.figure()
    else:
        figure = ax.figure
    n = len(products)
    if n > 0:
        gs = gridspec.GridSpec(1 + max(int(math.ceil(n / float(nperrow))), 1), nperrow)
    else:
        gs = gridspec.GridSpec(1 + (n / nperrow), nperrow)
    ax = figure.add_subplot(gs[0, :])
    if scan.is_profile:
        draw_raw(scan.arrays, ax=ax)
    if scan.peak_set is None:
        scan.pick_peaks()
    draw_peaklist(scan.peak_set, ax=ax, lw=0.5, alpha=0.75)
    if scan.deconvoluted_peak_set:
        draw_peaklist(scan.deconvoluted_peak_set, ax=ax, lw=0.5, alpha=0.75)
    ax.set_title(scan.id)
    k = -1
    # for the ith row of the grid
    for i in range(int(n // nperrow) + 1):
        # for the jth column of the row
        for j in range(nperrow):
            # the kth MS^n scan
            k += 1
            try:
                product_scan = products[k]
            except IndexError:
                # the number of MS^n scans is not a multiple of nperrow
                # so the last row has fewer than nperrow
                break
            # add the subplot from the grid spec using i + 1 instead of i
            # because the MS1 scan is drawn across the row at i = 0
            ax = figure.add_subplot(gs[i + 1, j])
            if scan.is_profile:
                draw_raw(scan.arrays, ax=ax, alpha=0.8)
            # obtain the interval around the precursor
            pinfo = product_scan.precursor_information
            if not product_scan.isolation_window.is_empty():
                lower, upper = (product_scan.isolation_window.lower_bound - 2,
                                product_scan.isolation_window.upper_bound + 2)
            else:
                lower = pinfo.mz - 4
                upper = pinfo.mz + 4
            try:
                peak = max(scan.peak_set.between(lower - 1.2, upper + 1.2), key=lambda x: x.intensity)
                local_intensity = peak.intensity
            except ValueError:
                if scan.deconvoluted_peak_set:
                    try:
                        peak = max(scan.deconvoluted_peak_set.between(lower - 1.2, upper + 1.2, use_mz=True),
                                   key=lambda x: x.intensity)
                        local_intensity = peak.intensity
                    except ValueError:
                        local_intensity = 1e3
                else:
                    local_intensity = 1e3

            # get the monoisotopic peak for the precursor, or the isolation center depending
            # upon whether the precursor has been deconvoluted and what the instrument reports
            if pinfo.extracted_charge != 0:
                target_mz = pinfo.extracted_mz
            else:
                target_mz = pinfo.mz

            draw_peaklist(scan.peak_set, ax=ax, alpha=0.5, lw=0.5)
            if scan.deconvoluted_peak_set:
                draw_peaklist(
                    scan.deconvoluted_peak_set.between(
                        lower - 1.2, upper + 1.2, use_mz=True),
                    ax=ax, alpha=0.9, color='blue')

            if label:
                label_peaks(scan, lower, upper, ax=ax,
                            is_deconvoluted=bool(scan.deconvoluted_peak_set))
            ax.set_ylim(0, local_intensity * 1.25)
            ax.set_xlim(lower, upper)
            upper_intensity = local_intensity

            # draw the precursor isolation annotations
            ax.vlines(target_mz, 0, upper_intensity * 1.5, alpha=0.50,
                      color='red', lw=1)
            if product_scan.isolation_window.lower != 0:
                ax.vlines(product_scan.isolation_window.lower_bound, 0,
                          upper_intensity * 1.5, linestyle='--', alpha=0.5)
            if product_scan.isolation_window.upper != 0:
                ax.vlines(product_scan.isolation_window.upper_bound, 0,
                          upper_intensity * 1.5, linestyle='--', alpha=0.5)
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.set_ylabel("")
            ax.set_xlabel("")
            if pinfo.extracted_charge == 0:
                if pinfo.charge == "ChargeNotProvided":
                    charge = '?'
                else:
                    charge = str(pinfo.charge)
            else:
                charge = str(pinfo.extracted_charge)
            ax.set_title("%0.3f @ %s" % (target_mz, charge))
    fig = figure
    # this should probably be a function of figwidth and number of rows
    fig.set_figheight(fig.get_figheight() * 2)
    fig.tight_layout()
    return ax


def annotate_scan_single(scan, product_scan, ax=None, label=True, standalone=True):
    '''Draw a zoomed-in view of the MS1 spectrum ``scan`` surrounding the
    area around each precursor ion that gave rise to ``product_scan``
    with monoisotopic peaks and isolation windows marked.

    .. plot::
        :include-source:

        import ms_deisotope
        from ms_deisotope import plot
        from ms_deisotope.test.common import datafile

        reader = ms_deisotope.MSFileLoader(datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
        bunch = next(reader)

        bunch.precursor.pick_peaks()
        bunch.precursor.deconvolute(
            scorer=ms_deisotope.PenalizedMSDeconVFitter(20., 2.0),
            averagine=ms_deisotope.glycopeptide, use_quick_charge=True)

        ax = plot.annotate_scan_single(bunch.precursor, bunch.products[0])
        ax.figure.set_figwidth(12)



    Parameters
    ----------
    scan: ScanBase
        The MS1 scan to annotate
    product_scan: ScanBase
        The product scan to annotate the precursor ion of
    ax: :class:`matplotlib._axes.Axes`
        An :class:`~.Axes` object to draw the plot on

    Returns
    -------
    :class:`matplotlib._axes.Axes`
    '''
    if ax is None:
        _, ax = plt.subplots(1)

    if scan.is_profile:
        draw_raw(scan.arrays, ax=ax)
    if scan.peak_set is None:
        scan.pick_peaks()
    draw_peaklist(scan.peak_set, ax=ax, lw=0.5, alpha=0.75)
    if scan.deconvoluted_peak_set:
        draw_peaklist(scan.deconvoluted_peak_set, ax=ax, lw=0.5, alpha=0.75)

    pinfo = product_scan.precursor_information
    if not product_scan.isolation_window.is_empty():
        lower, upper = (product_scan.isolation_window.lower_bound - 2,
                        product_scan.isolation_window.upper_bound + 2)
    else:
        lower = pinfo.mz - 4
        upper = pinfo.mz + 4

    peak_set = scan.peak_set
    try:
        peak = max(peak_set.between(lower - 1.2, upper + 1.2), key=lambda x: x.intensity)
        local_intensity = peak.intensity
    except ValueError:
        if scan.deconvoluted_peak_set:
            try:
                peak = max(scan.deconvoluted_peak_set.between(lower - 1.2, upper + 1.2, use_mz=True),
                           key=lambda x: x.intensity)
                local_intensity = peak.intensity
            except ValueError:
                local_intensity = 1e3
        else:
            local_intensity = 1e3

    # get the monoisotopic peak for the precursor, or the isolation center depending
    # upon whether the precursor has been deconvoluted and what the instrument reports
    if pinfo.extracted_charge != 0:
        target_mz = pinfo.extracted_mz
    else:
        target_mz = pinfo.mz

    draw_peaklist(scan.peak_set, ax=ax, alpha=0.5, lw=0.5)
    if scan.deconvoluted_peak_set:
        draw_peaklist(
            scan.deconvoluted_peak_set.between(
                lower - 1.2, upper + 1.2, use_mz=True),
            ax=ax, alpha=0.9, color='blue')

    if label:
        label_peaks(scan, lower, upper, ax=ax,
                    is_deconvoluted=bool(scan.deconvoluted_peak_set))

    ax.set_ylim(0, local_intensity * 1.25)
    ax.set_xlim(lower, upper)
    upper_intensity = local_intensity

    # draw the precursor isolation annotations
    ax.vlines(target_mz, 0, upper_intensity * 1.5, alpha=0.50,
              color='red', lw=1)
    if product_scan.isolation_window.lower != 0:
        ax.vlines(product_scan.isolation_window.lower_bound, 0,
                  upper_intensity * 1.5, linestyle='--', alpha=0.5)
    if product_scan.isolation_window.upper != 0:
        ax.vlines(product_scan.isolation_window.upper_bound, 0,
                  upper_intensity * 1.5, linestyle='--', alpha=0.5)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_ylabel("")
    ax.set_xlabel("")
    if pinfo.extracted_charge == 0:
        if pinfo.charge == "ChargeNotProvided":
            charge = '?'
        else:
            charge = str(pinfo.charge)
    else:
        charge = str(pinfo.extracted_charge)
    ax.set_title("%0.3f @ %s" % (target_mz, charge))

    return ax


def annotate_isotopic_peaks(scan, ax=None, color_cycle=None, **kwargs):
    '''Mark distinct isotopic peaks from the :class:`~.DeconvolutedPeakSet`
    in ``scan``.

    .. plot::
        :include-source:

        import ms_deisotope
        from ms_deisotope import plot
        from ms_deisotope.test.common import datafile

        reader = ms_deisotope.MSFileLoader(datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
        bunch = next(reader)

        bunch.precursor.pick_peaks()
        bunch.precursor.deconvolute(
            scorer=ms_deisotope.PenalizedMSDeconVFitter(20., 2.0),
            averagine=ms_deisotope.glycopeptide, use_quick_charge=True)

        ax = plot.draw_peaklist(bunch.precursor, color='black')
        ax = plot.annotate_isotopic_peaks(bunch.precursor, ax=ax)
        ax.set_xlim(1160, 1165)
        ax.figure.set_figwidth(12)


    Parameters
    ----------
    scan: ScanBase
        The scan to annotate
    color_cycle: :class:`~.Iterable`
        An iterable to draw isotopic cluster colors from
    ax: :class:`matplotlib._axes.Axes`
        An :class:`~.Axes` object to draw the plot on

    Returns
    -------
    :class:`matplotlib._axes.Axes`
    '''
    from .peak_set import DeconvolutedPeakSet

    if ax is None:
        _, ax = plt.subplots(1)
    if color_cycle is None:
        color_cycle = _default_color_cycle()
    color_cycle = itertools.cycle(color_cycle)
    if isinstance(scan, DeconvolutedPeakSet):
        peaks = scan
    else:
        peaks = getattr(scan, "deconvoluted_peak_set", [])
    peaks = sorted(peaks, key=lambda x: x.mz)
    for peak in peaks:
        color = next(color_cycle)
        draw_peaklist(peak.envelope, ax=ax, color=color, alpha=0.75, **kwargs)
        ax.scatter(*zip(*peak.envelope), color=color, alpha=0.75)
    return ax


def label_peaks(scan, min_mz=None, max_mz=None, ax=None, is_deconvoluted=None, threshold=None, **kwargs):
    """Label a region of the peak list, marking centroids with their m/z or mass. If the peaks
    of `scan` have been deconvoluted, the most abundant peak will be annotated with
    "<neutral mass> (<charge>)", otherwise just "<m/z>".

    Parameters
    ----------
    scan : :class:`~.ScanBase`
        The scan to annotate
    min_mz : float, optional
        The minimum m/z to annotate
    max_mz : float, optional
        The maximum m/z to annotate
    ax: :class:`matplotlib._axes.Axes`
        An :class:`~.Axes` object to draw the plot on
    is_deconvoluted : bool, optional
        Whether or not to always use :attr:`Scan.deconvoluted_peak_set`
    threshold : float, optional
        The intensity threshold under which peaks will be ignored

    Returns
    -------
    ax: :class:`matplotlib._axes.Axes`
        The axes the plot was drawn on
    annotations: :class:`list` of :class:`matplotlib.text.Text`
        The list of :class:`matplotlib.text.Text` annotations
    """
    if ax is None:
        # when no axes are provided, draw all the peak data
        _, ax = plt.subplots(1)
        if scan.is_profile:
            draw_raw(scan, ax)
        if scan.peak_set is not None:
            draw_peaklist(scan.peak_set, ax)
        if scan.deconvoluted_peak_set is not None:
            draw_peaklist(scan.deconvoluted_peak_set, ax)
    if is_deconvoluted is None:
        is_deconvoluted = scan.deconvoluted_peak_set is not None
    if min_mz is None:
        min_mz = 0
    if max_mz is None:
        if is_deconvoluted:
            max_mz = max(peak.mz for peak in scan.deconvoluted_peak_set)
        else:
            max_mz = max(peak.mz for peak in scan.peak_set)
    annotations = []
    # select the peak sub range
    if is_deconvoluted:
        subset = scan.deconvoluted_peak_set.between(
            min_mz, max_mz, use_mz=True)
    else:
        subset = scan.peak_set.between(
            min_mz, max_mz)
    if not subset:
        return
    # guess the threshold
    if threshold is None:
        threshold = 0.0
        if is_deconvoluted:
            threshold_list = ([max(i.intensity for i in p.envelope)
                               for p in subset])
        else:
            threshold_list = ([p.intensity for p in subset])
        if threshold_list:
            threshold = np.mean(threshold_list)
        threshold_list = [v > threshold for v in threshold_list]
        if threshold_list:
            threshold = np.mean(threshold_list)
    # draw the actual labels
    kwargs.setdefault("clip_on", True)
    kwargs.setdefault("fontsize", 10)
    if is_deconvoluted:
        for peak in subset:
            if peak.intensity > threshold:
                label = '%0.2f (%d)' % (peak.neutral_mass, peak.charge)
                # set the y-position to the highest peak in the isotopic
                # pattern
                pt = max(peak.envelope, key=lambda x: x.intensity)
                y = pt.intensity * 1.05
                # set the x-position to the weighted average m/z in the
                # isotopic pattern
                x = np.average(
                    [p.mz for p in peak.envelope],
                    weights=[p.intensity for p in peak.envelope])
                annotations.append(
                    ax.text(x, y, label, ha='center', **kwargs))
    else:
        for peak in subset:
            if peak.intensity > threshold:
                label = "%0.2f" % (peak.mz, )
                y = peak.intensity * 1.05
                x = peak.mz
                annotations.append(
                    ax.text(x, y, label, ha='center', **kwargs))
    return annotations, ax


def draw_spectrum_paths(graph, edge_color='red', peak_color='orange', alpha=0.8, fontsize=12, ax=None, **kwargs):
    """Draw all the paths given by `graph` over a rendered peak list.

    This function will draw all peaks which an edge connects in `peak_color`, and
    draws edges as lines connecting peaks in `edge_color`, with their annotations
    written across them.

    If a peak is not connected to an edge, it will not be drawn.

    Parameters
    ----------
    graph : :class:`~.SpectrumGraph` or :class:`list` of :class:`~.spectrum_graph.Path`
        The paths to annotate. If a :class:`~.SpectrumGraph` is given, the top 100 longest
        paths will be enumerated from it.
    edge_color : str, optional
        The color to draw the edge lines in (the default is 'red')
    peak_color : str, optional
        The color to draw the connected peaks in (the default is 'orange')
    alpha : float, optional
        The alpha channel value for peaks (the default is 0.8)
    fontsize : int, optional
        The font size to use for edge annotations (the default is 12)
    ax : :class:`matplotlib._axes.Axes`, optional
        The axes to draw the plot on (the default is None, which will cause a new figure with a
        single axes to be created)

    Returns
    -------
    :class:`matplotlib._axes.Axes`
    """
    if ax is None:
        _, ax = plt.subplots(1)
    try:
        paths = graph.longest_paths(limit=100)
    except AttributeError:
        paths = list(graph)

    for path in paths:
        for edge in path:
            for p1 in edge.start:
                for p2 in edge.end:
                    _draw_peak_pair(
                        (p1, p2), edge_color, peak_color, alpha, fontsize,
                        label=edge.annotation, ax=ax, **kwargs)


def _draw_peak_pair(pair, edge_color='red', peak_color='orange', alpha=0.8, fontsize=12, label=None, rotation=45,
                    ax=None, **kwargs):
    p1, p2 = pair
    ax.plot((p1.mz, p2.mz), (p1.intensity, p2.intensity),
            color=edge_color, alpha=alpha, **kwargs)
    kwargs.setdefault("clip_on", False)
    clip_on = kwargs['clip_on']
    draw_peaklist(pair, ax=ax, alpha=0.4, color=peak_color)
    if label:
        midx = (p1.mz + p2.mz) / 2
        # interpolate the midpoint's height
        midy = (p1.intensity * (p2.mz - midx) +
                p2.intensity * (midx - p1.mz)) / (p2.mz - p1.mz)

        # find the angle of the line connecting the two peaks
        xlo = min(p1.mz, p2.mz)
        xhi = max(p1.mz, p2.mz)
        adj = xhi - xlo
        ylo = min(p1.intensity, p2.intensity)
        yhi = max(p1.intensity, p2.intensity)
        opp = yhi - ylo
        hypot = np.hypot(adj, opp) # pylint: disable=assignment-from-no-return
        rotation = np.arccos(adj / hypot)  # pylint: disable=assignment-from-no-return

        if isinstance(label, (list, tuple)):
            label = '-'.join(map(str, label))
        else:
            label = str(label)
        ax.text(midx, midy, label, fontsize=fontsize,
                ha='center', va='bottom', rotation=rotation, clip_on=clip_on)
