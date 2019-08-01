import click
from six import string_types as basestring

import os
import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import (FigureCanvasAgg,)

from array import array

import ms_deisotope
from ms_deisotope import plot


@click.group('draw')
def draw():
    '''Draw images related to mass spectrometry data!
    '''

def _make_figure():
    figure = Figure(dpi=120)
    canvas = FigureCanvasAgg(figure)
    axis = figure.add_subplot(111)
    return figure, axis


@draw.command("tic")
@click.argument("path", type=click.Path(exists=True))
@click.option("-o", "--output-path", type=click.Path(writable=True),
              required=False, help="Path to write image to")
@click.option("-s", "--start-time", type=float, required=False, help='Time to start from')
@click.option("-e", "--end-time", type=float, required=False, help='Time to end at')
def draw_tic(path, output_path=None, start_time=None, end_time=None):
    """Draw the Total Ion Chromatogram (TIC), the total signal at each time point.
    """
    if output_path is None:
        output_path = path + '.tic.png'
    if start_time is None:
        start_time = 0
    if end_time is None:
        end_time = float('inf')

    figure, axis = _make_figure()

    reader = ms_deisotope.MSFileLoader(path)
    reader.start_from_scan(rt=start_time, grouped=False)

    time = array('d')
    intensity = array('d')

    bar = click.progressbar(reader, item_show_func=lambda x: str(
        x.id) if x is not None else '')
    with bar:
        for scan in bar:
            if scan.ms_level != 1:
                continue
            time.append(scan.scan_time)
            intensity.append(scan.arrays.intensity.sum())
    click.echo("Total Ion Current: %e" % np.sum(intensity))

    axis.plot(time, intensity)
    axis.set_xlabel("Scan Time (Min)", fontsize=16)
    axis.set_ylabel("Relative Intensity", fontsize=16)
    ylim = axis.get_ylim()
    axis.set_ylim(-10, ylim[1])
    axis.set_xlim(time[0] - 2, time[-1] + 2)
    figure.text(0.15, 0.8, "%0.3e" % np.sum(intensity), ha='left')
    figure.savefig(output_path, bbox_inches='tight', dpi=120)


