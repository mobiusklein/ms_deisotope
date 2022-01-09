"""The main entry point for the `ms-deisotope` command line program.
"""
import os

import click

import ms_peak_picker

from ms_peak_picker.scan_filter import parse as parse_filter

import ms_deisotope
from ms_deisotope import MSFileLoader
from ms_deisotope.data_source import RandomAccessScanSource

from ms_deisotope.tools.utils import processes_option, AveragineParamType, is_debug_mode, register_debug_hook
from ms_deisotope.tools.deisotoper import workflow, output



def configure_iterator(loader, start_time, end_time):
    """Configure `loader` to run between the time ranges specified, or
    if the `loader` does not support random access, the bounds are set to
    the start and end of the iterator.

    This function will also set the iteration mode of `loader` to "grouped".

    Parameters
    ----------
    loader : :class:`ScanDataSource`
        The source to load scans from
    start_time : float
        The start time to load scans from
    end_time : float
        The end time to stop loading scans at

    Returns
    -------
    str: start_scan_id
        The id value of the scan to start from, will be :cosnt:`None` if there is no random access.
    float: start_scan_time
        The actual time associated with the start scan
    str: end_scan_id
        The id value of the scan to stop at, will be :cosnt:`None` if there is no random access.
    float: end_scan_time
        The actual time associated with the end scan
    """
    if isinstance(loader, RandomAccessScanSource):
        last_scan = loader[len(loader) - 1]
        last_time = last_scan.scan_time

        start_scan = loader.get_scan_by_time(start_time)
        if loader.has_ms1_scans():
            start_scan = loader._locate_ms1_scan(start_scan)

        if end_time > last_time:
            end_time = last_time
        end_scan = loader.get_scan_by_time(end_time)
        if loader.has_ms1_scans():
            end_scan = loader._locate_ms1_scan(end_scan)

        start_scan_id = start_scan.id
        end_scan_id = end_scan.id

        start_scan_time = start_scan.scan_time
        end_scan_time = end_scan.scan_time

        loader.reset()
        loader.start_from_scan(start_scan_id, grouped=True)
    else:
        click.secho("The file format provided does not support random"
                    " access, start and end points will be ignored", fg='yellow')
        start_scan_time = 0
        start_scan_id = None

        end_scan_time = float('inf')
        end_scan_id = None
        loader.make_iterator(grouped=True)

    return start_scan_id, start_scan_time, end_scan_id, end_scan_time


def check_if_profile(loader):
    """Check if the spectra in `loader` are stored in profile mode,
    established by querying the next scan bunch.

    Parameters
    ----------
    loader : :class:`ScanIterator`
        The source to load scans from

    Returns
    -------
    :class:`bool`:
        Whether or not the spectra are profile mode.
    """
    first_bunch = next(loader)
    if first_bunch.precursor is not None:
        is_profile = (first_bunch.precursor.is_profile)
    elif first_bunch.products:
        is_profile = (first_bunch.products[0].is_profile)

    if is_profile:
        click.secho("Spectra are profile")
    else:
        click.secho("Spectra are centroided", fg='yellow')
    return is_profile


def determine_output_format(output_file):
    base, ext = os.path.splitext(output_file)
    if ext.lower() == '.mzml':
        return output.ThreadedMzMLScanStorageHandler
    elif ext.lower() == '.mgf':
        return output.ThreadedMGFScanStorageHandler
    else:
        click.secho("Could not infer output file format. Assuming mzML.", fg='yellow')
        return output.ThreadedMzMLScanStorageHandler


@click.command("deisotope",
               short_help=(
                   "Convert raw mass spectra data into deisotoped neutral mass peak lists written to mzML."
                   " Can accept mzML, mzXML, MGF with either profile or centroided scans."),
               context_settings=dict(help_option_names=['-h', '--help']))
@click.argument("ms-file", type=click.Path(exists=True))
@click.argument("outfile-path", type=click.Path(writable=True))
@click.option("-a", "--averagine", default=["peptide"],
              type=AveragineParamType(),
              help=('Averagine model to use for MS1 scans. '
                    'Either a name or formula. May specify multiple times.'),
              multiple=True)
@click.option("-an", "--msn-averagine", default=["peptide"],
              type=AveragineParamType(),
              help=('Averagine model to use for MS^n scans. '
                    'Either a name or formula. May specify multiple times.'),
              multiple=True)
@click.option("-s", "--start-time", type=float, default=0.0,
              help='Scan time to begin processing at in minutes')
@click.option("-e", "--end-time", type=float, default=float('inf'),
              help='Scan time to stop processing at in minutes')
@click.option("-c", "--maximum-charge", type=int, default=8,
              help=('Highest absolute charge state to consider'))
@click.option("-n", "--name", default=None,
              help="Name for the sample run to be stored. Defaults to the base name of the input data file")
@click.option("-t", "--score-threshold", type=float, default=workflow.SampleConsumer.MS1_SCORE_THRESHOLD,
              help="Minimum score to accept an isotopic pattern fit in an MS1 scan")
@click.option("-tn", "--msn-score-threshold", type=float, default=workflow.SampleConsumer.MSN_SCORE_THRESHOLD,
              help="Minimum score to accept an isotopic pattern fit in an MS^n scan")
@click.option("-m", "--missed-peaks", type=int, default=3,
              help="Number of missing peaks to permit before an isotopic fit is discarded")
@click.option("-mn", "--msn-missed-peaks", type=int, default=1,
              help="Number of missing peaks to permit before an isotopic fit is discarded in an MSn scan")
@processes_option
@click.option("-b", "--background-reduction", type=float, default=0., help=(
    "Background reduction factor. Larger values more aggresively remove low abundance"
    " signal in MS1 scans."))
@click.option("-bn", "--msn-background-reduction", type=float, default=0., help=(
    "Background reduction factor. Larger values more aggresively remove low abundance"
    " signal in MS^n scans."))
@click.option("-r", '--transform', multiple=True, type=parse_filter, help=(
    "Scan transformations to apply to MS1 scans. May specify more than once."))
@click.option("-rn", '--msn-transform', multiple=True, type=parse_filter, help=(
    "Scan transformations to apply to MS^n scans. May specify more than once."))
@click.option("-v", "--extract-only-tandem-envelopes", is_flag=True, default=False,
              help='Only work on regions that will be chosen for MS/MS')
@click.option("--verbose", is_flag=True, help="Log additional diagnostic information for each scan.")
@click.option("-g", "--ms1-averaging", default=0, type=int, help=(
    "The number of MS1 scans before and after the current MS1 "
    "scan to average when picking peaks."))
@click.option("--ignore-msn", is_flag=True, default=False, help="Ignore MS^n scans")
@click.option("-i", "--isotopic-strictness", default=2.0, type=float)
@click.option("-in", "--msn-isotopic-strictness", default=0.0, type=float)
@click.option("-snr", "--signal-to-noise-threshold", default=1.0, type=float, help=(
    "Signal-to-noise ratio threshold to apply when filtering peaks"))
@click.option("-mo", "--mass-offset", default=0.0, type=float, help=("Shift peak masses by the given amount"))
@click.option("-D", "--default-precursor-ion-selection-window", default=1.5, type=float,
              help="The isolation window width to assume when it is not specified.")
def deisotope(ms_file, outfile_path, averagine=None, start_time=None, end_time=None, maximum_charge=None,
              name=None, msn_averagine=None, score_threshold=35., msn_score_threshold=10., missed_peaks=1,
              msn_missed_peaks=1, background_reduction=0., msn_background_reduction=0.,
              transform=None, msn_transform=None, processes=4, extract_only_tandem_envelopes=False,
              ignore_msn=False, isotopic_strictness=2.0, ms1_averaging=0,
              msn_isotopic_strictness=0.0, signal_to_noise_threshold=1.0, mass_offset=0.0,
              default_precursor_ion_selection_window=1.5, deconvolute=True, verbose=False):
    '''Convert raw mass spectra data into deisotoped neutral mass peak lists written to mzML.
    '''
    if transform is None:
        transform = []
    if msn_transform is None:
        msn_transform = []

    if (ignore_msn and extract_only_tandem_envelopes):
        click.secho(
            "Cannot use both --ignore-msn and --extract-only-tandem-envelopes",
            fg='red')
        raise click.Abort("Cannot use both --ignore-msn and --extract-only-tandem-envelopes")

    cache_handler_type = determine_output_format(outfile_path)
    click.echo("Preprocessing %s" % ms_file)
    minimum_charge = 1 if maximum_charge > 0 else -1
    charge_range = (minimum_charge, maximum_charge)

    loader = MSFileLoader(ms_file)
    (start_scan_id, start_scan_time,
     end_scan_id, end_scan_time) = configure_iterator(loader, start_time, end_time)

    is_profile = check_if_profile(loader)

    if name is None:
        name = os.path.splitext(os.path.basename(ms_file))[0]

    if os.path.exists(outfile_path) and not os.access(outfile_path, os.W_OK):
        click.secho("Can't write to output file path", fg='red')
        raise click.Abort()

    click.secho("Initializing %s" % name, fg='green')
    click.echo("from %s (%0.2f) to %s (%0.2f)" % (
        start_scan_id, start_scan_time, end_scan_id, end_scan_time))
    if deconvolute:
        click.echo("charge range: %s" % (charge_range,))

    if is_profile:
        ms1_peak_picking_args = {
            "transforms": [
            ] + list(transform),
            "signal_to_noise_threshold": signal_to_noise_threshold
        }
        if background_reduction:
            ms1_peak_picking_args['transforms'].append(
                ms_peak_picker.scan_filter.FTICRBaselineRemoval(
                    scale=background_reduction, window_length=2))
            ms1_peak_picking_args['transforms'].append(ms_peak_picker.scan_filter.SavitskyGolayFilter())
    else:
        ms1_peak_picking_args = {
            "transforms": [
                ms_peak_picker.scan_filter.FTICRBaselineRemoval(
                    scale=background_reduction, window_length=2),
            ] + list(transform)
        }

    if msn_background_reduction > 0.0:
        msn_peak_picking_args = {
            "transforms": [
                ms_peak_picker.scan_filter.FTICRBaselineRemoval(
                    scale=msn_background_reduction, window_length=2),
            ] + list(msn_transform)
        }
    else:
        msn_peak_picking_args = {
            "transforms": [
            ] + list(msn_transform)
        }

    if mass_offset != 0.0:
        ms1_peak_picking_args['transforms'].append(
            ms_peak_picker.scan_filter.RecalibrateMass(offset=mass_offset))
        msn_peak_picking_args['transforms'].append(
            ms_peak_picker.scan_filter.RecalibrateMass(offset=mass_offset))

    if deconvolute:
        if len(averagine) == 1:
            averagine = averagine[0]
            ms1_deconvoluter_type = ms_deisotope.deconvolution.AveraginePeakDependenceGraphDeconvoluter
        else:
            ms1_deconvoluter_type = ms_deisotope.deconvolution.MultiAveraginePeakDependenceGraphDeconvoluter

        ms1_deconvolution_args = {
            "scorer": ms_deisotope.scoring.PenalizedMSDeconVFitter(score_threshold, isotopic_strictness),
            "max_missed_peaks": missed_peaks,
            "averagine": averagine,
            "truncate_after": workflow.SampleConsumer.MS1_ISOTOPIC_PATTERN_WIDTH,
            "ignore_below": workflow.SampleConsumer.MS1_IGNORE_BELOW,
            "deconvoluter_type": ms1_deconvoluter_type,
            "use_quick_charge": True
        }

        if msn_isotopic_strictness >= 1:
            msn_isotopic_scorer = ms_deisotope.scoring.PenalizedMSDeconVFitter(
                msn_score_threshold, msn_isotopic_strictness)
        else:
            msn_isotopic_scorer = ms_deisotope.scoring.MSDeconVFitter(msn_score_threshold)

        if len(msn_averagine) == 1:
            msn_averagine = msn_averagine[0]
            msn_deconvoluter_type = ms_deisotope.deconvolution.AveraginePeakDependenceGraphDeconvoluter
        else:
            msn_deconvoluter_type = ms_deisotope.deconvolution.MultiAveraginePeakDependenceGraphDeconvoluter

        msn_deconvolution_args = {
            "scorer": msn_isotopic_scorer,
            "averagine": msn_averagine,
            "max_missed_peaks": msn_missed_peaks,
            "truncate_after": workflow.SampleConsumer.MSN_ISOTOPIC_PATTERN_WIDTH,
            "ignore_below": workflow.SampleConsumer.MSN_IGNORE_BELOW,
            "use_quick_charge": True,
            "deconvoluter_type": msn_deconvoluter_type,
        }
    else:
        ms1_deconvolution_args = None
        msn_deconvolution_args = None

    consumer = workflow.SampleConsumer(
        ms_file,
        ms1_peak_picking_args=ms1_peak_picking_args,
        ms1_deconvolution_args=ms1_deconvolution_args,
        msn_peak_picking_args=msn_peak_picking_args,
        msn_deconvolution_args=msn_deconvolution_args,
        storage_path=outfile_path, sample_name=name,
        start_scan_id=start_scan_id, storage_type=cache_handler_type,
        end_scan_id=end_scan_id, n_processes=processes,
        extract_only_tandem_envelopes=extract_only_tandem_envelopes,
        ignore_tandem_scans=ignore_msn,
        ms1_averaging=ms1_averaging,
        default_precursor_ion_selection_window=default_precursor_ion_selection_window,
        deconvolute=deconvolute,
        verbose=verbose,
        end_scan_time=end_scan_time,
        start_scan_time=start_scan_time)
    consumer.start()


if is_debug_mode():
    register_debug_hook()

if __name__ == '__main__':
    deisotope.main()
