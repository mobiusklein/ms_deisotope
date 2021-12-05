'''Implements the infrastructure to spread indexing tasks over a single
:class:`~.RandomAccessScanSource` across multiple processes for speed.

The primary function in this module is :func:`index`, which will perform
this dispatch.
'''
import multiprocessing
import logging

import dill

from .scan_index import ExtendedScanIndex
from .scan_interval_tree import (
    ScanIntervalTree, extract_intervals, make_rt_tree)

logger = logging.getLogger(__name__)
n_cores = multiprocessing.cpu_count()


def indexing_iterator(reader, start, end, index):
    """A helper function which will iterate over an interval of a :class:`~.RandomAccessScanSource`
    while feeding each yielded :class:`~.ScanBunch` into a provided :class:`~.ExtendedScanIndex`.

    [description]

    Parameters
    ----------
    reader : :class:`~.RandomAccessScanSource` or :class:`~.ScanIterator`
        The scan data source to loop over.
    start : int
        The starting index
    end : int
        The stopping index
    index : :class:`~.ExtendedScanIndex`
        The scan index to fill.

    Yields
    ------
    :class:`~.ScanBunch`
    """
    assert end >= start, "End cannot precede Start"
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
    """The task function for :func:`quick_index`, which will build an
    :class:`~.ExtendedIndex` and :class:`ScanIntervalTree` from an index
    range over a :class:`~.ScanIterator`

    Parameters
    ----------
    reader : :class:`~.ScanIterator` or :class:`~.RandomAccessScanSource`
        The scan source to iterate over
    start : int
        The starting index
    end : int
        The stopping index

    Returns
    -------
    start: int
        The starting index for this chunk
    end: int
        The stopping index for this chunk
    index: :class:`~.ExtendedIndex`
        The constructed scan metadata index for this chunk
    intervals: :class:`~.ScanIntervalTree
        The constructed scan interval tree for this chunk
    """
    index = ExtendedScanIndex()
    iterator = indexing_iterator(reader, start, end, index)
    intervals = extract_intervals(iterator)
    return (start, end, index, intervals)


def partition_work(n_items, n_workers, start_index=0, reader=None):
    """Given an index range and a number of workers to work on them,
    break the index range into approximately evenly sized sub-intervals.

    This is a helper function for :func:`run_task_in_chunks` used to
    compute the chunks.

    Parameters
    ----------
    n_items : int
        The maximum value of the index range
    n_workers : int
        The number of workers to split the work between
    start_index : int, optional
        The starting value of the index range (the default is 0)

    Returns
    -------
    list:
        A list of (start, end) pairs defining the index ranges each worker will
        handle.
    """
    if n_workers == 1:
        return [[start_index, start_index + n_items]]
    chunk_size = int(n_items / n_workers)
    n_items += start_index
    intervals = []
    start = start_index
    end = start + chunk_size
    if reader is not None:
        if reader.has_fast_random_access:
            scan = reader.get_scan_by_index(start)
            if scan.ms_level > 1:
                start_scan = reader.find_previous_ms1(start)
                if start_scan is not None:
                    start = start_scan.index
            end = min(end, len(reader) - 1)
            scan = scan = reader.get_scan_by_index(end)
            if scan.ms_level > 1:
                end_scan = reader.find_next_ms1(end)
                if end_scan is not None:
                    end = end_scan.index
    intervals.append([start, end])
    start += chunk_size

    while start  < n_items:
        end = start + chunk_size
        if reader is not None:
            scan = reader.get_scan_by_index(start)
            if scan.ms_level > 1:
                start_scan = reader.find_previous_ms1(start)
                if start_scan is not None:
                    start = start_scan.index
            end = min(end, len(reader))
            if end != len(reader):
                scan = scan = reader.get_scan_by_index(end)
                if scan.ms_level > 1:
                    end_scan = reader.find_next_ms1(end)
                    if end_scan is not None:
                        end = end_scan.index
        intervals.append([start, end])
        start = end
    # Make sure that the last chunk actually covers the end of
    # the interval.
    last = intervals[-1]
    if last[1] > n_items:
        last[1] = n_items
    else:
        last[1] += (n_items - last[1])
    return intervals


class _Indexer(object):
    '''A pickle-able callable object which wraps :func:`index_chunk` for the call signature
    used by :func:`run_task_in_chunks`
    '''
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

    def __getstate__(self):
        state = {
            "task": dill.dumps(self.task)
        }
        return state

    def __setstate__(self, state):
        self.task = dill.loads(state['task'])

    def __call__(self, payload):
        _reader, start, _end = payload
        out = self.task(payload)
        return start, out


class _TaskPayload(object):
    """A wrapper for the input to the distributed task which transmits the :class:`~.RandomAccessScanSource`
    via :mod:`dill` to make pickling of any wrapped file objects possible.

    Mocks a subset of the :class:`~.Sequence` API to allow it to be treated like a :class:`tuple`

    Attributes
    ----------
    reader: :class:`~.RandomAccessScanSource`
        The scan data source to be shared with the worker.
    start: int
        The scan index to start processing from
    end: int
        The scan index to stop processing at.
    options: dict
        A dictionary of extra arguments that the task might use.
    """
    def __init__(self, reader, start, end, **kwargs):
        self.reader = reader
        self.start = start
        self.end = end
        self.options = kwargs

    def __iter__(self):
        yield self.reader
        yield self.start
        yield self.end

    def __getitem__(self, i):
        if i == 0:
            return self.reader
        elif i == 1:
            return self.start
        elif i == 2:
            return self.end
        else:
            raise IndexError(i)

    def __len__(self):
        return 3

    def __getstate__(self):
        state = {
            "reader": dill.dumps(self.reader, -1),
            "start": self.start,
            "end": self.end,
            "options": self.options
        }
        return state

    def __setstate__(self, state):
        self.reader = dill.loads(state['reader'])
        self.start = state['start']
        self.end = state['end']
        self.options = state['options']

    def __repr__(self):
        template = "{self.__class__.__name__}({self.reader}, {self.start}, {self.end}, {self.options})"
        return template.format(self=self)


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
    scan_ranges = partition_work(
        n_items, n_chunks, start_scan, reader=reader)
    feeder = (
        _TaskPayload(reader, scan_range[0], scan_range[1])
        for scan_range in scan_ranges)
    result = []
    for i, block in enumerate(pool.imap_unordered(_TaskWrapper(task), feeder), 1):
        result.append(block)
        if progress_indicator is not None:
            progress_indicator(i / float(n_chunks))
    pool.close()
    pool.join()
    result.sort(key=lambda x: x[0])
    result = [x[1] for x in result]
    return result


def _merge_indices(indices):
    index = indices[0]
    for ind in indices:
        index = index.merge(ind)
    return index


def _make_interval_tree(intervals):
    concat = []
    for i in intervals:
        concat.extend(i)
    return ScanIntervalTree(make_rt_tree(concat), None)


def index(reader, n_processes=4, scan_interval=None, progress_indicator=None):
    """Generate a :class:`~.ExtendedScanIndex` and :class:`~.ScanIntervalTree` for
    `reader` between `scan_interval` start and end points across `n_processes` worker
    processes.

    If a :class:`~.ScanIterator` is passed instead of a :class:`~.RandomAccessScanSource`,
    only a single process will be used.

    Parameters
    ----------
    reader : :class:`~.RandomAccessScanSource` or :class:`~.ScanIterator`
        The scan data source to index.
    n_processes : int, optional
        The number of worker processes to use (the default is 4 or however many CPUs are available, whichever is lower)
    scan_interval : tuple, optional
        The start and stop scan indices to operate on (the default is None, which will index the entire file)
    progress_indicator : :class:`Callable`, optional
        A callable object which will be used to report progress as chunks finish processing. It must take
        one argument, a :class:`float` which represents the fraction of all work completed.

    Returns
    -------
    :class:`~.ExtendedScanIndex`:
        The extended metadata index for this data file.
    :class:`~.ScanIntervalTree`:
        The scan interval tree for this data file.

    See Also
    --------
    run_task_in_chunks
    """
    task = _Indexer()
    # indexing a :class:`ScanIterator` without random access, have to go in sequence
    if not hasattr(reader, 'start_from_scan'):
        logger.info("A non-random access ScanIterator was passed, defaulting to a single worker.")
        reader.make_iterator(grouped=True)
        chunks = [task(_TaskPayload(reader, 0, len(reader)))]
        n_processes = 1
    else:
        chunks = run_task_in_chunks(
            reader, n_processes, scan_interval=scan_interval, task=task,
            progress_indicator=progress_indicator)
    indices = [chunk[2] for chunk in chunks]
    intervals = [chunk[3] for chunk in chunks]
    index = _merge_indices(indices)
    interval_tree = _make_interval_tree(intervals)
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
