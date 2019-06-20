import multiprocessing

from .scan_index import ExtendedScanIndex
from .scan_interval_tree import (
    ScanIntervalTree, extract_intervals, make_rt_tree)


n_cores = multiprocessing.cpu_count()


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


class _TaskWrapper(object):
    '''A simple wrapper for a callable to capture the index range for a chunk created by
    :func:`run_task_in_chunks`, so that the start index of the chunk is known.
    '''
    def __init__(self, task):
        self.task = task

    def __call__(self, payload):
        _reader, start, _end = payload
        out = self.task(payload)
        return start, out


def run_task_in_chunks(reader, n_processes=None, n_chunks=None, scan_interval=None, task=None, progress_indicator=None):
    """Run a :class:`~.Callable` `task` over a :class:`~.ScanIterator` in chunks across multiple processes.

    This function breaks apart a :class:`~.ScanIterator`'s scans over `scan_interval`,
    or the whole sequence if not provided.

    Parameters
    ----------
    reader : :class:`~.ScanIterator`
        The set of :class:`~.Scan` objects to operate on
    n_processes : int, optional
        The number of worker processes to use (the default is 4 or the number of cores available,
        whichever is lower)
    n_chunks : int, optional
        The number of chunks to break the scan range into (the default is equal to `n_processes`)
    scan_interval : :class:`tuple` of (:class:`int`, :class:`int`), optional
        The start and stop scan index to apply the task over. If omitted, the entire scan range will
        be used. If either entry is :const:`None`, then the index will be assumed to be the first or
        last scan respectively.
    task : :class:`~.Callable`
        The callable object which will be executed on each chunk in a sub-process. It must take one
        argument, a :class:`tuple` of (`reader`, `start index`, `stop index`), and it must return a pickle-able
        object.
    progress_indicator : :class:`~.Callable`, optional
        A callable object which will be used to report progress as chunks finish processing. It must take
        one argument, a :class:`float` which represents the fraction of all work completed.

    Returns
    -------
    :class:`list`:
        The result of `task` on each chunk of `reader` in index sorted order.

    """
    if n_processes is None:
        n_processes = min(n_cores, 4)
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
    result = []
    for i, block in enumerate(pool.imap_unordered(_TaskWrapper(task), feeder), 1):
        result.append(block)
        if progress_indicator is not None:
            progress_indicator(i / float(n_chunks))
    pool.close()
    pool.terminate()
    pool.join()
    result.sort(key=lambda x: x[0])
    result = [x[1] for x in result]
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


def index(reader, n_processes=4, scan_interval=None, progress_indicator=None):
    task = _Indexer()
    # indexing a :class:`ScanIterator` without random access, have to go in sequence
    if not hasattr(reader, 'start_from_scan'):
        reader.make_iterator(grouped=True)
        chunks = [task((reader, 0, len(reader)))]
    else:
        chunks = run_task_in_chunks(
            reader, n_processes, scan_interval=scan_interval, task=task,
            progress_indicator=progress_indicator)
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
