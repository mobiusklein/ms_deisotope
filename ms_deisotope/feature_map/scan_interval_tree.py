from collections import namedtuple
from ms_deisotope.peak_dependency_network import Interval, IntervalTreeNode


class BoundingBox(namedtuple("BoundingBox", ['mz', 'rt'])):
    def merge(self, other):
        min_mz = min(self[0].start, other[0].start)
        max_mz = max(self[0].end, other[0].end)
        min_rt = min(self[1].start, other[1].start)
        max_rt = max(self[1].end, other[1].end)
        return self.__class__(Interval(min_mz, max_mz), Interval(min_rt, max_rt))

    def overlaps(self, other):
        return self.mz.overlaps(other.mz) and self.rt.overlaps(other.rt)

    def __add__(self, other):
        return self.merge(other)


def extract_intervals(scan_iterator, time_radius=5., mz_lower=2., mz_higher=3.):
    intervals = []
    for scan, products in scan_iterator:
        for product in products:
            intervals.append(BoundingBox(
                Interval(max(0, product.precursor_information.mz - mz_lower),
                         product.precursor_information.mz + mz_higher),
                Interval(max(0, product.scan_time - time_radius),
                         product.scan_time + time_radius)))
    return intervals


def combine_intervals(a, b):
    min_mz = min(a[0].start, b[0].start)
    max_mz = max(a[0].end, b[0].end)
    min_rt = min(a[1].start, b[1].start)
    max_rt = max(a[1].end, b[1].end)
    return (Interval(min_mz, max_mz), Interval(min_rt, max_rt))


def overlaps_2d(a, b):
    return a[0].overlaps(b[0]) and a[1].overlaps(b[1])


def merge_interval_set(intervals, minimum_overlap_size=0.3):
    merged_intervals = []
    for interval in intervals:
        for i, candidate in enumerate(merged_intervals):
            if interval.overlaps(candidate) and interval.mz.overlap_size(candidate.mz) > (
                    (interval.mz.end - interval.mz.start) * minimum_overlap_size):
                merged_intervals[i] = interval.merge(candidate)
                break
        else:
            merged_intervals.append(interval)
    return merged_intervals


def make_rt_tree(intervals):
    temp = []
    for interval in intervals:
        mz, rt = interval
        rt.members.append(mz)
        temp.append(rt)
    tree = IntervalTreeNode.build(temp)
    return tree


class ScanIntervalTree(object):
    @classmethod
    def build(cls, scan_iterator, time_radius=3., mz_lower=2., mz_higher=3.):
        intervals = extract_intervals(
            scan_iterator, time_radius=time_radius,
            mz_lower=mz_lower, mz_higher=mz_higher)
        return cls(make_rt_tree(intervals), intervals)

    def __init__(self, rt_tree, original_intervals=None):
        self.rt_tree = rt_tree
        self.original_intervals = original_intervals

    def get_mz_intervals_for_rt(self, rt_point):
        if self.rt_tree is None:
            return [[0, float('inf')]]
        intervals = [
            i.members[0] for i in self.rt_tree.contains_point(rt_point)
        ]
        intervals = sorted([[i.start, i.end] for i in intervals])
        out = []
        if len(intervals) == 0:
            return []
        last = intervals[0]
        for interval in intervals[1:]:
            if interval[0] < last[1]:
                last[1] = interval[1]
            else:
                out.append(last)
                last = interval
        out.append(last)
        return out

    def __call__(self, scan):
        intervals = self.get_mz_intervals_for_rt(scan.scan_time)
        return intervals
