import os
import math

from collections import Counter

import ms_deisotope

import click

from ms_deisotope.feature_map import quick_index
from ms_deisotope.feature_map import scan_interval_tree
from ms_deisotope.data_source import _compression


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@cli.command("byte-index")
@click.argument("path", type=click.Path(exists=True))
def byte_index(path):
    reader = ms_deisotope.MSFileLoader(path, use_index=False)
    try:
        fn = reader.prebuild_byte_offset_file
    except AttributeError:
        click.echo("\"%s\" does not support pre-indexing byte offsets" % (path,))
        return
    fn(path)


@cli.command("metadata-index")
@click.argument("path", type=click.Path(exists=True))
def metadata_index(path, processes=4):
    reader = ms_deisotope.MSFileLoader(path)
    index, interval_tree = quick_index.index(reader, processes)
    name = path
    index_file_name = index.index_file_name(name)
    with open(index_file_name, 'w') as fh:
        index.serialize(fh)


def partial_ms_file_iterator(reader, start, end):
    for scan_bunch in reader.start_from_scan(index=start):
        if scan_bunch.precursor.index > end:
            break
        yield scan_bunch


class MSMSIntervalTask(object):
    def __init__(self, time_radius, mz_lower, mz_higher):
        self.time_radius = time_radius
        self.mz_lower = mz_lower
        self.mz_higher = mz_higher

    def __call__(self, payload):
        reader, start, end = payload
        iterator = partial_ms_file_iterator(reader, start, end)
        intervals = scan_interval_tree.extract_intervals(
            iterator, time_radius=self.time_radius,
            mz_lower=self.mz_lower, mz_higher=self.mz_higher)
        return intervals


@cli.command("msms-intervals")
@click.argument('paths', type=click.Path(exists=True), nargs=-1)
@click.option("-o", "--output", type=click.Path(writable=True, file_okay=True, dir_okay=False), required=False)
def msms_intervals(paths, processes=4, time_radius=5, mz_lower=2., mz_higher=3., output=None):
    interval_extraction = MSMSIntervalTask(time_radius, mz_lower, mz_higher)
    interval_set = []
    total_work_items = len(paths) * processes * 4

    def run():
        for path in paths:
            reader = ms_deisotope.MSFileLoader(path)
            chunk_out_of_order = quick_index.run_task_in_chunks(
                reader, processes, processes * 4, task=interval_extraction)
            for chunk in chunk_out_of_order:
                interval_set.extend(chunk)
                yield 0
    work_iterator = run()
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
        click.secho("Building Index", fg='yellow', err=True)
        index = quick_index.ExtendedScanIndex()
        reader.reset()
        for bunch in reader:
            index.add_scan_bunch(bunch)
        reader.reset()
        with open(index_file_name, 'w') as fh:
            index.serialize(fh)
    else:
        with open(index_file_name, 'rt') as fh:
            index = quick_index.ExtendedScanIndex.deserialize(fh)
    return reader, index


@cli.command("charge-states")
@click.argument("path", type=click.Path(exists=True))
def charge_states(path):
    reader, index = _ensure_metadata_index(path)

    charges = Counter()
    for msn_id, msn_info in index.msn_ids.items():
        charges[msn_info.charge] += 1
    for charge in sorted(charges, key=abs):
        click.echo("%d: %d" % (charge, charges[charge]))


def binsearch(array, x):
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


@cli.command("precursor-clustering")
@click.argument("path", type=click.Path(exists=True))
def precursor_clustering(path, grouping_error=2e-5):
    reader, index = _ensure_metadata_index(path)
    points = []
    for msn_id, msn_info in index.msn_ids.items():
        points.append((msn_info.neutral_mass, msn_info.intensity))
    points.sort(key=lambda x: x[1], reverse=1)
    centroids = []
    if len(points) == 0:
        click.secho("No MS/MS detected", fg='yellow', err=True)
        return

    for point in points:
        if len(centroids) == 0:
            centroids.append((point[0], [point]))
            continue
        i = binsearch(centroids, point[0])
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


if _compression.has_idzip:

    @cli.command("idzip")
    @click.argument('path', type=str)
    @click.option("-o", "--output", type=click.Path(writable=True, file_okay=True, dir_okay=False), required=False)
    def idzip_compression(path, output):
        if output is None:
            output = '-'
        with click.open_file(output, mode='wb') as outfh:
            writer = _compression.GzipFile(fileobj=outfh, mode='wb')
            with click.open_file(path, 'rb') as infh:
                buffer_size = 2 ** 28
                chunk = infh.read(buffer_size)
                while chunk:
                    writer.write(chunk)
                    chunk = infh.read(buffer_size)
            writer.close()


main = cli.main

if __name__ == '__main__':
    main()
