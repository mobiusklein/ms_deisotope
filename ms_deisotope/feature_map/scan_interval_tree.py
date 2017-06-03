from ms_deisotope.peak_dependency_network import Interval, IntervalTreeNode


def extract_intervals(scan_iterator, time_radius=5., mz_lower=2., mz_higher=3.):
    intervals = []
    for scan, products in scan_iterator:
        for product in products:
            intervals.append((
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
            if overlaps_2d(interval, candidate) and interval[0].overlap_size(
                    candidate[0]) > (
                    interval[0].end - interval[0].start * minimum_overlap_size):
                merged_intervals[i] = combine_intervals(candidate, interval)
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
