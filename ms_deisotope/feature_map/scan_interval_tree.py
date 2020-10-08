import json

from collections import namedtuple, deque

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
            if product.precursor_information is None:
                continue
            intervals.append(BoundingBox(
                Interval(max(0, product.precursor_information.mz - mz_lower),
                         product.precursor_information.mz + mz_higher),
                Interval(max(0, product.scan_time - time_radius),
                         product.scan_time + time_radius)))
    return intervals


def nest_2d_intervals(intervals):
    out = []
    for interval in intervals:
        mz, rt = interval
        rt.members.append(mz)
        out.append(rt)
    return out


def make_rt_tree(intervals):
    nested = nest_2d_intervals(intervals)
    tree = IntervalTreeNode.build(nested)
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

    def serialize(self, handle):
        work_queue = deque()
        id_counter = 0
        node_id_map = {}

        def serialize_contained(interval_list):
            out = []
            for interval in interval_list:
                item = {"start": interval.start, "end": interval.end}
                inner = [{"start": i.start, "end": i.end}
                         for i in interval.members]
                item['members'] = inner
                out.append(item)
            return out

        root = self.rt_tree
        node_dict = {
            "id": id_counter,
            "contained": serialize_contained(root.contained),
            "center": root.center,
            "start": root.start,
            "end": root.end,
            "left": None,
            "right": None
        }
        root_id = node_dict['id']
        node_id_map[node_dict['id']] = node_dict
        id_counter += 1
        work_queue.append((root.left, node_dict['id'], 'left'))
        work_queue.append((root.right, node_dict['id'], 'right'))

        while work_queue:
            node, parent_id, side = work_queue.popleft()
            if node is None:
                continue
            node_dict = {
                "id": id_counter,
                "contained": serialize_contained(node.contained),
                "center": node.center,
                "start": node.start,
                "end": node.end,
                "left": None,
                "right": None
            }
            node_id_map[node_dict['id']] = node_dict
            id_counter += 1
            parent = node_id_map[parent_id]
            parent[side] = node_dict
            work_queue.append((node.left, node_dict['id'], 'left'))
            work_queue.append((node.right, node_dict['id'], 'right'))

        json.dump(node_id_map[root_id], handle)

    @classmethod
    def load(cls, handle):
        data = json.load(handle)

        def load_contained(interval_list):
            out = []
            for item in interval_list:
                inner = [Interval(**i) for i in item['members']]
                out.append(Interval(start=item['start'], end=item['end'], members=inner))
            return out

        def load_node_single(node_dict):
            contained = node_dict.get("contained", [])
            contained = load_contained(contained)
            node = IntervalTreeNode(
                node_dict['center'], left=None, contained=contained,
                right=None)
            node.start = node_dict['start']
            node.end = node_dict['end']
            return node
        root = load_node_single(data)

        work_queue = deque()

        work_queue.append((data['left'], root, 'left'))
        work_queue.append((data['right'], root, 'right'))

        while work_queue:
            node_dict, parent, side = work_queue.popleft()
            if node_dict is None:
                continue
            node = load_node_single(node_dict)
            node.parent = parent
            node.level = parent.level + 1
            if side == 'left':
                parent.left = node
            elif side == 'right':
                parent.right = node
            else:
                raise ValueError(side)
            work_queue.append((node_dict['left'], node, 'left'))
            work_queue.append((node_dict['right'], node, 'right'))
        return cls(root)
