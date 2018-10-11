# pragma: no cover
import math
import itertools

from ms_peak_picker.plot import draw_peaklist, draw_raw
from matplotlib import pyplot as plt, gridspec
has_plot = True


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


def annotate_scan(scan, products, nperrow=4, ax=None):
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
                peak = max(scan.peak_set.between(lower + 1.2, upper - 1.2), key=lambda x: x.intensity)
                local_intensity = peak.intensity
            except ValueError:
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


def annotate_scan_single(scan, product_scan, ax=None, standalone=True):
    if ax is None:
        fig, ax = plt.subplots(1)

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

    try:
        peak = max(scan.peak_set.between(lower + 1.2, upper - 1.2), key=lambda x: x.intensity)
        local_intensity = peak.intensity
    except ValueError:
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
    from .peak_set import DeconvolutedPeakSet

    if ax is None:
        fig, ax = plt.subplots(1)
    if color_cycle is None:
        color_cycle = _default_color_cycle()
    color_cycle = itertools.cycle(color_cycle)
    if isinstance(scan, DeconvolutedPeakSet):
        peaks = scan
    else:
        peaks = getattr(scan, "deconvoluted_peak_set", [])
    for peak in peaks:
        color = next(color_cycle)
        draw_peaklist(peak.envelope, ax=ax, color=color, alpha=0.75, **kwargs)
        ax.scatter(*zip(*peak.envelope), color=color, alpha=0.75)
    return ax
