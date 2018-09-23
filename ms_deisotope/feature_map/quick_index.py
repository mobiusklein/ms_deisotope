import multiprocessing
from .scan_index import ExtendedScanIndex
from .scan_interval_tree import (
    ScanIntervalTree, extract_intervals, make_rt_tree)


def indexing_iterator(reader, start, end, index):
    try:
        iterator = reader.start_from_scan(index=start, grouped=True)
    except AttributeError:
        if start != 0:
            raise
        iterator = reader
    for scan_bunch in iterator:
        try:
            ix = scan_bunch.precursor.index
        except AttributeError:
            ix = scan_bunch.products[0].index
        if ix > end:
            break
        index.add_scan_bunch(scan_bunch)
        yield scan_bunch


def index_chunk(reader, start, end):
    index = ExtendedScanIndex()
    iterator = indexing_iterator(reader, start, end, index)
    intervals = extract_intervals(iterator)
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
        try:
            result = index_chunk(reader, start, end)
        except Exception as e:
            print(reader, start, end, e)
            import traceback
            traceback.print_exc()
            raise e
        return result


def run_task_in_chunks(reader, n_processes=4, n_chunks=None, scan_interval=None, task=None):
    if task is None or not callable(task):
        raise ValueError("The task must be callable!")
    if scan_interval is None:
        start_scan = 0
        end_scan = len(reader.index)
    else:
        start_scan, end_scan = scan_interval
        if start_scan is None:
            start_scan = 0
        if end_scan is None:
            end_scan = len(reader.index)
    if n_chunks is None:
        n_chunks = n_processes
    n_items = end_scan - start_scan
    pool = multiprocessing.Pool(n_processes)
    scan_ranges = partition_work(n_items, n_chunks, start_scan)
    feeder = ((reader, scan_range[0], scan_range[1]) for scan_range in scan_ranges)
    result = list(pool.imap_unordered(task, feeder))
    pool.close()
    pool.terminate()
    pool.join()
    return result


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
    task = _Indexer()
    # indexing a :class:`ScanIterator` without random access, have to go in sequence
    if not hasattr(reader, 'start_from_scan'):
        reader.make_iterator(grouped=True)
        chunks = [task((reader, 0, len(reader)))]
    else:
        chunks = run_task_in_chunks(
            reader, n_processes, scan_interval=scan_interval, task=task)
    indices = [chunk[2] for chunk in chunks]
    intervals = [chunk[3] for chunk in chunks]
    index = merge_indices(indices)
    interval_tree = make_interval_tree(intervals)
    return index, interval_tree


def multi_index(readers, n_processes=4, scan_interval=None):
    index_collection = []
    interval_collection = []
    for reader in readers:
        chunks = run_task_in_chunks(
            reader, n_processes, scan_interval=scan_interval, task=_Indexer())
        indices = [chunk[2] for chunk in chunks]
        intervals = [chunk[3] for chunk in chunks]
        index_collection.append(indices)
        interval_collection.append(intervals)
