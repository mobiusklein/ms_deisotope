'''A collection of little command line utilities for inspecting mass
spectrum data.
'''
import io
import os
import math
import csv
import sys
from collections import Counter

import click
import six

import numpy as np

import ms_deisotope

from ms_deisotope.feature_map import quick_index
from ms_deisotope.feature_map import scan_interval_tree

from ms_deisotope.clustering.scan_clustering import (
    iterative_clustering, ScanClusterWriter, ScanClusterReader, _DynamicallyLoadingResolver)

from ms_deisotope.qc.isolation import isolation_window_valid, is_isolation_window_empty

from ms_deisotope.data_source import (
    _compression, ScanProxyContext, MSFileLoader)
from ms_deisotope.data_source.scan import RandomAccessScanSource
from ms_deisotope.data_source.metadata.file_information import SourceFile

from ms_deisotope.output import ProcessedMzMLDeserializer

from ms_deisotope.tools import conversion, draw, maintenance
from ms_deisotope.tools.utils import processes_option, is_debug_mode, progress, register_debug_hook, spinner


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    '''A collection of utilities for inspecting and manipulating
    mass spectrometry data.
    '''


@cli.command("describe", short_help=("Produce a minimal textual description"
                                     " of a mass spectrometry data file"))
@click.argument('path', type=click.Path(exists=True))
@click.option("-d", "--diagnostics", is_flag=True,
              help="Run more diagnostics, greatly increasing runtime but producing additional informatoin")
def describe(path, diagnostics=False):
    '''Produces a minimal textual description of a mass spectrometry data file.
    '''
    click.echo("Describing \"%s\"" % (path,))
    try:
        sf = SourceFile.from_path(path)
    except IOError as err:
        raise click.ClickException("Could not open file \"%s\":\n%s" % (path, err))

    if sf.file_format is None:
        raise click.ClickException("\"%s\" doesn't appear to be a mass spectrometry data file" % (path, ))
    click.echo("File Format: %s" % (sf.file_format, ))
    click.echo("ID Format: %s" % (sf.id_format, ))
    reader = MSFileLoader(path)
    if isinstance(reader, RandomAccessScanSource):
        click.echo("Format Supports Random Access: True")
        first_scan = reader[0]
        last_scan = reader[-1]
        click.echo("First Scan: %s at %0.3f minutes" % (first_scan.id, first_scan.scan_time))
        click.echo("Last Scan: %s at %0.3f minutes" % (last_scan.id, last_scan.scan_time))
    else:
        click.echo("Format Supports Random Access: False")
    try:
        finfo = reader.file_description()
        click.echo("Contents:")
        for key in finfo.contents:
            click.echo("    %s" % (key, ))
    except AttributeError:
        pass
    index_file_name = quick_index.ExtendedScanIndex.index_file_name(path)
    # Extra introspection if the extended index is available
    if os.path.exists(index_file_name):
        with open(index_file_name, 'rt') as fh:
            index = quick_index.ExtendedScanIndex.deserialize(fh)
        ms1_scans = len(index.ms1_ids)
        msn_scans = len(index.msn_ids)
        click.echo("MS1 Scans: %d" % (ms1_scans, ))
        click.echo("MSn Scans: %d" % (msn_scans, ))
        n_defaulted = 0
        n_orphan = 0

        charges = Counter()
        first_msn = float('inf')
        last_msn = 0
        for scan_info in index.msn_ids.values():
            n_defaulted += scan_info.get('defaulted', False)
            n_orphan += scan_info.get('orphan', False)
            charges[scan_info['charge']] += 1
            rt = scan_info['scan_time']
            if rt < first_msn:
                first_msn = rt
            if rt > last_msn:
                last_msn = rt
        click.echo("First MSn Scan: %0.3f minutes" % (first_msn,))
        click.echo("Last MSn Scan: %0.3f minutes" % (last_msn,))
        for charge, count in sorted(charges.items()):
            if not isinstance(charge, int):
                continue
            click.echo("Precursors with Charge State %d: %d" % (charge, count))
        if n_defaulted > 0:
            click.echo("Defaulted MSn Scans: %d" % (n_defaulted,))
        if n_orphan > 0:
            click.echo("Orphan MSn Scans: %d" % (n_orphan,))
    if diagnostics:
        reader.reset()
        scan_ids_with_invalid_isolation_window = []
        scan_ids_with_empty_isolation_window = []
        activation_counter = Counter()
        n_ms1 = 0
        n_msn = 0
        for precursor, products in reader:
            if precursor is not None:
                n_ms1 += 1
            for product in products:
                n_msn += 1
                if not isolation_window_valid(product):
                    scan_ids_with_invalid_isolation_window.append((precursor.id, product.id))
                if is_isolation_window_empty(product):
                    scan_ids_with_empty_isolation_window.append(
                        (precursor.id, product.id))
                activation_counter[product.activation] += 1

        click.echo("MS1 Spectra: %d" % (n_ms1, ))
        click.echo("MSn Spectra: %d" % (n_msn, ))
        click.echo("Invalid Isolation Windows: %d" % (
            len(scan_ids_with_invalid_isolation_window), ))
        click.echo("Empty Isolation Windows: %d" % (
            len(scan_ids_with_empty_isolation_window), ))

        click.echo("Activation Methods")
        for activation, count in activation_counter.items():
            if not activation.is_multiple_dissociation():
                click.echo("\t{method}:{energy} = {count}".format(
                    method=activation.method, energy=activation.energy, count=count))
            else:
                click.echo("\t{method}:{energy} = {count}".format(
                    method=','.join(map(str, activation.method)),
                    energy=','.join(map(str, activation.energies)), count=count))


@cli.command("byte-index", short_help='Build an external byte offset index for a mass spectrometry data file')
@click.argument('paths', type=click.Path(exists=True), nargs=-1)
def byte_index(paths):
    '''Build an external byte offset index for a mass spectrometry data file, saving time when
    opening the file with indexing enabled.

    Supported Formats: mzML, mzXML
    '''
    for path in paths:
        click.echo("Indexing %s" % (path, ))
        reader = ms_deisotope.MSFileLoader(path, use_index=False)
        try:
            fn = reader.prebuild_byte_offset_file
        except AttributeError:
            click.echo("\"%s\" does not support pre-indexing byte offsets" % (path,))
            return
        fn(path)


@cli.command("metadata-index", short_help='Build an external scan metadata index for a mass spectrometry data file')
@click.argument('paths', type=click.Path(exists=True), nargs=-1)
@click.option("-D", '--deconvoluted', is_flag=True, help='Whether or not to assume this file has been '
              'processed specifically by ms_deisotope\'s deconvolution algorithm')
@processes_option
def metadata_index(paths, processes=4, deconvoluted=False):
    '''Build an external scan metadata index for a mass spectrometry data file

    This extended index is saved in a separate JSON file that can be loaded with
    :class:`~.ExtendedScanIndex`. It includes the scan time of all scans, the precursor
    mass of MSn scans, as well as the relationships between precursor and product ion
    scans, as well as other details. See :class:`~.ExtendedScanIndex` for more information
    '''
    for path in paths:
        click.echo("Indexing %s" % (path, ))
        if deconvoluted:
            reader = ProcessedMzMLDeserializer(path, use_extended_index=False)
        else:
            reader = MSFileLoader(path)
        try:
            fn = reader.prebuild_byte_offset_file
            if not reader.source._check_has_byte_offset_file():
                fn(path)
        except AttributeError:
            pass
        if processes > 1:
            progbar = click.progressbar(label='Building Index', length=100)
            acc = [0]

            def update_bar(x):
                '''Progress Bar update callback for :func:`~.quick_index.index`
                '''
                x = int(x * 100)
                x -= acc[0]  # pylint: disable=cell-var-from-loop
                progbar.update(x)  # pylint: disable=cell-var-from-loop
                acc[0] += x  # pylint: disable=cell-var-from-loop

            with progbar:
                update_bar(0.0)
                index, _ = quick_index.index(
                    reader, processes, progress_indicator=update_bar)
        else:
            index = quick_index.ExtendedScanIndex()
            reader.reset()
            try:
                n = len(reader)
                progbar = click.progressbar(label='Building Index', length=n)
            except TypeError:
                progbar = spinner(title="Building Index")
            with progbar:
                for bunch in reader.make_iterator(grouped=True):
                    i = 0
                    i += bunch.precursor is not None
                    i += len(bunch.products)
                    index.add_scan_bunch(bunch)
                    progbar.update(i)

        name = path
        index_file_name = index.index_file_name(name)
        with open(index_file_name, 'w') as fh:
            index.serialize(fh)


def _partial_ms_file_iterator(reader, start, end):
    for scan_bunch in reader.start_from_scan(index=start):
        if scan_bunch.precursor.index > end:
            break
        yield scan_bunch


class _MSMSIntervalTask(object):
    def __init__(self, time_radius, mz_lower, mz_higher):
        self.time_radius = time_radius
        self.mz_lower = mz_lower
        self.mz_higher = mz_higher

    def __call__(self, payload):
        reader, start, end = payload
        iterator = _partial_ms_file_iterator(reader, start, end)
        intervals = scan_interval_tree.extract_intervals(
            iterator, time_radius=self.time_radius,
            mz_lower=self.mz_lower, mz_higher=self.mz_higher)
        return intervals


@cli.command("msms-intervals", short_help=(
    'Build an interval tree over precursor isolation events in time and m/z space'))
@click.argument('paths', type=click.Path(exists=True), nargs=-1)
@click.option("-o", "--output", type=click.Path(writable=True, file_okay=True, dir_okay=False), required=False)
@processes_option
def msms_intervals(paths, processes=4, time_radius=5, mz_lower=2., mz_higher=3., output=None):
    '''Construct an interval tree spanning time and m/z domains where MSn spectra were acquired
    in the LC-MS map. The interval tree is serialized to JSON.
    '''
    interval_extraction = _MSMSIntervalTask(time_radius, mz_lower, mz_higher)
    interval_set = []
    total_work_items = len(paths) * processes * 4

    def _run():
        for path in paths:
            reader = MSFileLoader(path)
            chunk_out_of_order = quick_index.run_task_in_chunks(
                reader, processes, processes * 4, task=interval_extraction)
            for chunk in chunk_out_of_order:
                interval_set.extend(chunk)
                yield 0
    work_iterator = _run()
    with click.progressbar(work_iterator, length=total_work_items, label='Extracting Intervals') as g:
        for _ in g:
            pass
    tree = scan_interval_tree.ScanIntervalTree(scan_interval_tree.make_rt_tree(interval_set))
    if output is not None:
        with open(output, 'wt') as fh:
            tree.serialize(fh)
    else:
        stream = click.get_text_stream('stdout')
        tree.serialize(stream)
        stream.flush()


def _ensure_metadata_index(path):
    reader = ms_deisotope.MSFileLoader(path)
    name = path
    index_file_name = quick_index.ExtendedScanIndex.index_file_name(name)
    if not os.path.exists(index_file_name):
        click.secho("Building Index For %s" % (path, ), fg='yellow', err=True)
        index = quick_index.ExtendedScanIndex()
        reader.reset()
        try:
            n = len(reader)
            progbar = click.progressbar(label='Building Index', length=n)
        except TypeError:
            progbar = spinner(title="Building Index")
        with progbar:
            for bunch in reader.make_iterator(grouped=True):
                i = 0
                i += bunch.precursor is not None
                i += len(bunch.products)
                index.add_scan_bunch(bunch)
                progbar.update(i)
        reader.reset()
        with open(index_file_name, 'wt') as fh:
            index.serialize(fh)
    else:
        with open(index_file_name, 'rt') as fh:
            index = quick_index.ExtendedScanIndex.deserialize(fh)
    return reader, index


@cli.command("charge-states", short_help='Count the different precursor charge states in a mass spectrometry data file')
@click.argument("path", type=click.Path(exists=True))
def charge_states(path):
    '''Count the different precursor charge states in a mass spectrometry data file.

    This command will construct a metadata index if it is not found.
    '''
    _, index = _ensure_metadata_index(path)

    charges = Counter()
    for _, msn_info in index.msn_ids.items():
        charges[msn_info.charge] += 1

    def sort_key(k):
        if isinstance(k, int):
            return abs(k)
        else:
            return 0

    for charge in sorted(charges, key=sort_key):
        click.echo("%s: %d" % (charge, charges[charge]))


def _binsearch(array, x):
    n = len(array)
    lo = 0
    hi = n

    while hi != lo:
        mid = (hi + lo) / 2
        y = array[mid][0]
        err = y - x
        if hi - lo == 1:
            return mid
        elif err > 0:
            hi = mid
        else:
            lo = mid
    return


@cli.command("precursor-clustering", short_help='Cluster precursor masses in a mass spectrometry data file')
@click.argument("path", type=click.Path(exists=True))
def precursor_clustering(path, grouping_error=2e-5):
    '''Cluster precursor masses in a mass spectrometry data file.

    This command will construct a metadata index if it is not found.
    '''
    _, index = _ensure_metadata_index(path)
    points = []
    for _, msn_info in index.msn_ids.items():
        points.append((msn_info.neutral_mass, msn_info.intensity))
    points.sort(key=lambda x: x[1], reverse=1)
    centroids = []
    if not points:
        click.secho("No MS/MS detected", fg='yellow', err=True)
        return

    for point in points:
        if not centroids:
            centroids.append((point[0], [point]))
            continue
        i = _binsearch(centroids, point[0])
        centroid = centroids[i]
        err = (centroid[0] - point[0]) / point[0]
        if abs(err) < grouping_error:
            centroid[1].append(point)
        else:
            if err < 0:
                i += 1
            centroids.insert(i, (point[0], [point]))
    acc = 0
    nt = 0
    for centroid, obs in centroids:
        n = len(obs)
        if n == 1:
            continue
        mean = sum(p[0] for p in obs) / n
        var = sum([(p[0] - mean) ** 2 for p in obs]) / (n - 1)
        acc += (var * (n - 1))
        nt += (n - 1)
    for centroid, obs in centroids:
        click.echo("%f: %d" % (centroid, len(obs)))
    click.echo("MS/MS Precursor Mass Std. Dev.: %f Da" % (math.sqrt(acc / nt),))


@cli.command("spectrum-clustering", short_help=("Cluster MSn spectra in a mass spectrometry data file using"
                                                " cosine similarity"))
@click.argument('paths', type=click.Path(exists=True), nargs=-1)
@click.option("-m", "--precursor-error-tolerance", type=float, default=1e-5)
@click.option("-t", "--similarity-threshold", "similarity_thresholds", multiple=True, type=float)
@click.option("-o", "--output", "output_path", type=click.Path(writable=True, file_okay=True, dir_okay=False),
              required=False)
@click.option("-c", "--cache-size", type=int, default=2**10, help=(
    "The number of scans to cache in memory when not using --in-memory. If you are clustering multiple "
    "files, this number will be apply to each file separately."))
@click.option("-M", "--in-memory", is_flag=True, default=False, help=(
    "Whether to load the entire dataset into memory for better performance"))
@click.option("-D", "--deconvoluted", is_flag=True, default=False, help=(
    "Whether to assume the spectra are deconvoluted or not"))
def spectrum_clustering(paths, precursor_error_tolerance=1e-5, similarity_thresholds=None, output_path=None,
                        in_memory=False, deconvoluted=False, cache_size=2**10):
    '''Cluster spectra by precursor mass and cosine similarity.

    Spectrum clusters are written out to a text file recording
    cluster precursor mass, within-cluster similarity, and the
    source file and scan ID for each cluster member.
    '''
    if not similarity_thresholds:
        similarity_thresholds = [0.1, 0.4, 0.7]
    else:
        similarity_thresholds = sorted(similarity_thresholds)
    if output_path is None:
        output_path = "-"
    msn_scans = []
    n_spectra = 0

    with click.progressbar(paths, label="Indexing", item_show_func=lambda x: str(x) if x else '') as progbar:
        key_seqs = []
        for path in progbar:
            if deconvoluted:
                reader = ProcessedMzMLDeserializer(path)
                index = reader.extended_index
            else:
                reader, index = _ensure_metadata_index(path)
            key_seqs.append((reader, index))
            n_spectra += len(index.msn_ids)

    with click.progressbar(label="Loading Spectra", length=n_spectra,
                           item_show_func=lambda x: str(x) if x else '') as progbar:
        for reader, index in key_seqs:
            if not in_memory:
                if not reader.has_fast_random_access:
                    click.secho(
                        "%s does not have fast random access, scan fetching may be slow!" % (
                            reader, ), fg='yellow')
                proxy_context = ScanProxyContext(reader, cache_size=cache_size)
                pinfo_map = {
                    pinfo.product_scan_id: pinfo for pinfo in
                    index.get_precursor_information()
                }
                for i in index.msn_ids:
                    progbar.current_item = i
                    progbar.update(1)
                    scan = proxy_context(i)
                    scan.precursor_information = pinfo_map[i]
                    msn_scans.append(scan)
            else:
                if reader.has_fast_random_access:
                    # We have fast random access so we can just loop over the index and pull out
                    # the MSn scans directly without completely traversing the file.
                    for i in index.msn_ids:
                        progbar.current_item = i
                        progbar.update(1)
                        scan = reader.get_scan_by_id(i)
                        if scan.peak_set is None and not deconvoluted:
                            scan = scan.pick_peaks().pack(bind=True)
                        msn_scans.append(scan)
                else:
                    # If we don't  have fast random access, it's better just to loop over the file,
                    # and absorb the cost of parsing the MS1 scans
                    reader.reset()
                    reader.make_iterator(grouped=False)
                    for scan in reader:
                        if scan.ms_level != 1:
                            progbar.current_item = scan.id
                            progbar.update(1)
                            if scan.peak_set is None and not deconvoluted:
                                scan = scan.pick_peaks().pack(bind=True)
                            msn_scans.append(scan)
                # Dispose of the state that is no longer required.
                reader.reset()
                index.clear()


    click.echo("Begin Clustering", err=True)
    clusters = iterative_clustering(
        msn_scans, precursor_error_tolerance, similarity_thresholds)
    click.echo("Clusering Finished", err=True)
    by_size = Counter()
    for cluster in clusters:
        by_size[len(cluster)] += 1
    click.echo("Clusters: {:d}".format(len(clusters)), err=True)
    for key, value in sorted(by_size.items()):
        click.echo("Size {:d}: {:d}".format(key, value), err=True)
    with click.open_file(output_path, mode='w') as outfh:
        writer = ScanClusterWriter(outfh)
        for cluster in clusters:
            writer.save(cluster)


@cli.command('cluster-statistics')
@click.argument('path', type=click.Path(readable=True, dir_okay=False))
def cluster_evaluation(path):
    """Calculate some statistics describing the clustered mass spectra described in
    the specified cluster file.
    """
    with click.open_file(path) as fh:
        dirname = os.path.dirname(path)
        reader = ScanClusterReader(fh, _DynamicallyLoadingResolver(dirname))
        clusters = list(reader)
        by_size = Counter()
        for cluster in clusters:
            by_size[len(cluster)] += 1
        click.echo("Clusters: {:d}".format(len(clusters)))
        for key, value in sorted(by_size.items()):
            click.echo("Size {:d}: {:d}".format(key, value))



@cli.command('ms1-spectrum-diagnostics')
@click.argument('path', type=click.Path(exists=True, readable=True))
@click.option("-o", "--output-path", type=click.Path(writable=True))
def ms1_spectrum_diagnostics(path, output_path=None):
    '''Collect diagnostic information from MS1 spectra.
    '''
    reader = ms_deisotope.MSFileLoader(path)

    reader.make_iterator(grouped=True)

    ms1_metric_names = [
        'scan_id', 'scan_index', 'scan_time', 'duty_cycle', 'tic',
        'base_peak_mz', 'base_peak_intensity', 'data_point_count',
        'injection_time', 'n_ms2_scans'
    ]
    ms1_metrics = []
    products = None
    last_ms1 = None
    prog = progress(length=len(reader), label='Processing Scans',
                    file=sys.stderr, item_show_func=lambda x: x.id if x else '')
    with prog:
        for precursor, products in reader:
            ms1_time = precursor.scan_time
            if last_ms1 is not None:
                duty_cycle = ms1_time - last_ms1
                ms1_metrics[-1]['duty_cycle'] = duty_cycle
            last_ms1 = ms1_time
            bp = precursor.base_peak()
            acquisition_info = precursor.acquisition_information
            if acquisition_info:
                scan_event = acquisition_info[0]
                inj = scan_event.injection_time
            else:
                inj = np.nan
            ms1_record = {
                "scan_id": precursor.id,
                "scan_index": precursor.index,
                "scan_time": precursor.scan_time,
                "duty_cycle": np.nan,
                "tic": precursor.tic(),
                "base_peak_mz": bp.mz,
                "base_peak_intensity": bp.intensity,
                "data_point_count": precursor.arrays.mz.size,
                "injection_time": inj,
                "n_ms2_scans": len([p for p in products if p.ms_level == 2])
            }
            ms1_metrics.append(ms1_record)
            prog.current_item = precursor
            prog.update(1 + len(products))

    if last_ms1 is not None:
        if products:
            last_time = max([p.scan_time for p in products])
            duty_cycle = last_time - last_ms1
            ms1_metrics[-1]['duty_cycle'] = duty_cycle


    if output_path is None:
        outfh = click.open_file("-", mode='wb')
    else:
        outfh = io.open(output_path, mode='wb')
    if six.PY3:
        stream = io.TextIOWrapper(outfh, encoding='utf8', newline='')
    else:
        stream = outfh
    writer = csv.DictWriter(stream, fieldnames=ms1_metric_names)
    writer.writeheader()
    writer.writerows(ms1_metrics)
    stream.flush()

if _compression.has_idzip:
    @cli.command("idzip", short_help='Compress a file with idzip, a gzip-compatible format with random access support')
    @click.argument('path', type=click.Path(allow_dash=True, readable=True))
    @click.option("-o", "--output", type=click.Path(writable=True, file_okay=True, dir_okay=False),
                  required=False, help="The output stream to write to. Defaults to STDOUT if omitted.")
    def idzip_compression(path, output):
        '''Compress a file using  idzip, a gzip-compatible format with random access support.
        '''
        if output is None:
            output = '-'
        with click.open_file(output, mode='wb') as outfh:
            writer = _compression.GzipFile(fileobj=outfh, mode='wb')
            with click.open_file(path, 'rb') as infh:
                try:
                    infh_wrap = io.BufferedReader(infh)
                    header = infh_wrap.peek(2)
                    if _compression.starts_with_gz_magic(header):
                        click.echo("Detected gzip input file", err=True)
                        infh_wrap = _compression.GzipFile(fileobj=infh_wrap)
                except AttributeError:
                    infh_wrap = infh
                buffer_size = _compression.WRITE_BUFFER_SIZE
                chunk = infh_wrap.read(buffer_size)
                while chunk:
                    writer.write(chunk)
                    chunk = infh_wrap.read(buffer_size)
            writer.close()



def _mount_group(group):
    try:
        for name, command in group.commands.items():
            cli.add_command(command, name)
    except Exception as e:
        click.secho("%r occurred while loading additional tools from %r" % (e, group), err=True, fg='yellow')


_mount_group(conversion.ms_conversion)

cli.add_command(draw.draw)
cli.add_command(maintenance.maintenance)
main = cli.main

if is_debug_mode():
    register_debug_hook()

if __name__ == '__main__':
    click.secho("Running Debug Mode", fg='yellow')
    register_debug_hook()
    main()
