import multiprocessing
from .scan_index import ExtendedScanIndex
from .scan_interval_tree import (
    ScanIntervalTree, extract_intervals, make_rt_tree)


def indexing_generator(reader, start, end, index):
    for scan_bunch in reader.start_from_scan(index=start):
        if scan_bunch.precursor.index > end:
            break
        index.add_scan_bunch(scan_bunch)
        yield scan_bunch


def index_chunk(reader, start, end):
    index = ExtendedScanIndex()
    generator = indexing_generator(reader, start, end, index)
    intervals = extract_intervals(generator)
    return (start, end, index, intervals)


def partition_work(n_items, n_workers, start_index=0):
    chunk_size = int(n_items / n_workers)
    n_items += start_index
    intervals = []
    start = start_index
    intervals.append([start, start + chunk_size])
    start += chunk_size

    while start + chunk_size < (n_items):
        end = start + chunk_size
        intervals.append([start, end])
        start = end

    last = intervals[-1]
    if last[1] > n_items:
        last[1] = n_items
    else:
        last[1] += (n_items - last[1])

    return intervals


class _Indexer(object):
    def __call__(self, payload):
        reader, start, end = payload
        return index_chunk(reader, start, end)


def run_task_in_chunks(reader, n_processes=4, n_chunks=None, scan_interval=None):
    if scan_interval is None:
        start_scan = 0
        end_scan = len(reader.index)
    else:
        start_scan, end_scan = scan_interval
    if n_chunks is None:
        n_chunks = n_processes
    n_items = end_scan - start_scan
    pool = multiprocessing.Pool(n_processes)
    scan_ranges = partition_work(n_items, n_chunks, start_scan)
    feeder = ((reader, scan_range[0], scan_range[1]) for scan_range in scan_ranges)
    return list(pool.imap_unordered(_Indexer(), feeder))


def merge_indices(indices):
    index = indices[0]
    for ind in indices:
        index = index.merge(ind)
    return index


def make_interval_tree(intervals):
    concat = []
    for i in intervals:
        concat.extend(i)
    return ScanIntervalTree(make_rt_tree(concat), None)


def index(reader, n_processes=4, scan_interval=None):
    chunks = run_task_in_chunks(
        reader, n_processes, scan_interval=scan_interval)
    indices = [chunk[2] for chunk in chunks]
    intervals = [chunk[3] for chunk in chunks]
    index = merge_indices(indices)
    interval_tree = make_interval_tree(intervals)
    return index, interval_tree
