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
