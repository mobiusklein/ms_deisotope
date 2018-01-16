import os
from functools import partial

import ms_deisotope

import click

from ms_deisotope.feature_map import quick_index
from ms_deisotope.feature_map import scan_interval_tree


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
@click.argument("path")
def metadata_index(path, type=click.Path(exists=True), processes=4):
    reader = ms_deisotope.MSFileLoader(path)
    index, interval_tree = quick_index.index(reader, processes)
    name = os.path.splitext(path)[0]
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


main = cli.main
