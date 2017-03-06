from ms_deisotope.utils import uid

from collections import deque

from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_GetSlice
from ms_peak_picker._c.peak_set cimport FittedPeak

cimport numpy as np
import numpy as np


cdef double INF = float("inf")


cdef class LCMSFeatureTreeList(object):
    def __init__(self, roots=None):
        if roots is None:
            roots = []
        self.roots = list(roots)
        self._node_id_hash = None

    cdef void _invalidate(self):
        self._node_id_hash = None

    cpdef tuple find_time(self, double retention_time):
        if len(self.roots) == 0:
            raise ValueError()
        cdef:
            int lo, hi
            double rt
            LCMSFeatureTreeNode node
        lo = 0
        hi = PyList_GET_SIZE(self.roots)
        while lo != hi:
            i = (lo + hi) / 2
            node = <LCMSFeatureTreeNode>PyList_GET_ITEM(self.roots, i)
            rt = node.retention_time
            if rt == retention_time:
                return node, i
            elif (hi - lo) == 1:
                return None, i
            elif rt < retention_time:
                lo = i
            elif rt > retention_time:
                hi = i

    def _build_node_id_hash(self):
        self._node_id_hash = set()
        for node in self.unspool():
            self._node_id_hash.add(node.node_id)

    @property
    def node_id_hash(self):
        if self._node_id_hash is None:
            self._build_node_id_hash()
        return self._node_id_hash

    def insert_node(self, node):
        self._invalidate()
        try:
            root, i = self.find_time(node.retention_time)
            if root is None:
                if i != 0:
                    self.roots.insert(i + 1, node)
                else:
                    slot = self.roots[i]
                    if slot.retention_time < node.retention_time:
                        i += 1
                    self.roots.insert(i, node)
            else:
                root.add(node)
            return i
        except ValueError:
            self.roots.append(node)
            return 0

    def insert(self, retention_time, peaks):
        node = LCMSFeatureTreeNode(retention_time, peaks)
        return self.insert_node(node)

    def extend(self, iterable):
        for peaks, retention_time in iterable:
            self.insert(retention_time, peaks)

    def __getitem__(self, i):
        return self.roots[i]

    def __len__(self):
        return len(self.roots)

    def __iter__(self):
        return iter(self.roots)

    def __reduce__(self):
        return self.__class__, (self.nodes,)

    def clone(self):
        return LCMSFeatureTreeList(node.clone() for node in self)

    def unspool(self):
        out_queue = []
        for root in self:
            stack = [root]
            while len(stack) != 0:
                node = stack.pop()
                out_queue.append(node)
                stack.extend(node.children)
        return out_queue

    def common_nodes(self, other):
        return len(self.node_id_hash & other.node_id_hash)

    def __repr__(self):
        return "LCMSFeatureTreeList(%d nodes, %0.2f-%0.2f)" % (
            len(self), self[0].retention_time, self[-1].retention_time)


cdef class LCMSFeatureTreeNode(object):
    def __init__(self, retention_time, members=None):
        if members is None:
            members = []
        self.retention_time = retention_time
        self.members = members
        self._most_abundant_member = None
        self._mz = 0
        self._recalculate()
        self.node_id = uid()

    def clone(self):
        cdef LCMSFeatureTreeNode node
        node = self.__class__(
            self.retention_time, list(self.members))
        node.node_id = self.node_id
        return node

    def _unspool_strip_children(self):
        node = self.__class__(
            self.retention_time, list(self.members))
        yield node

    def _recalculate(self):
        self._calculate_most_abundant_member()
        self._mz = self._most_abundant_member.mz

    def _calculate_most_abundant_member(self):
        if PyList_GET_SIZE(self.members) == 1:
            self._most_abundant_member = <FittedPeak>PyList_GET_ITEM(self.members, 0)
        else:
            if PyList_GET_SIZE(self.members) == 0:
                self._most_abundant_member = None
            else:
                self._most_abundant_member = max(
                    self.members, key=lambda x: x.intensity)
        return

    def __getstate__(self):
        return (self.retention_time, self.members, self.node_id)

    def __setstate__(self, state):
        self.retention_time, self.members, self.node_id = state
        self._recalculate()

    def __reduce__(self):
        return self.__class__, (self.retention_time, []), self.__getstate__()

    @property
    def mz(self):
        if self._mz is None:
            if self._most_abundant_member is not None:
                self._mz = self._most_abundant_member.mz
        return self._mz

    def add(self, node, recalculate=True):
        if self.node_id == node.node_id:
            raise ValueError("Duplicate Node %s" % node)
        self.members.extend(node.members)
        if recalculate:
            self._recalculate()

    cpdef double _total_intensity_members(self):
        total = 0.
        for peak in self.members:
            total += peak.intensity
        return total

    cpdef double max_intensity(self):
        return self._most_abundant_member.intensity

    cpdef double total_intensity(self):
        return self._total_intensity_members()

    cpdef bint _eq(self, LCMSFeatureTreeNode other):
        return self.members == other.members

    cpdef bint _ne(self, LCMSFeatureTreeNode other):
        return not (self._eq(other))

    def __richcmp__(self, other, int code):
        if code == 2:
            return self._eq(other)
        elif code == 3:
            return self._ne(other)
        else:
            return NotImplemented

    def __hash__(self):
        return hash(self.uid)

    @property
    def peaks(self):
        peaks = list(self.members)
        return peaks
    
    def __repr__(self):
        return "%s(%f, %0.4f %s|%d)" % (
            self.__class__.__name__,
            self.retention_time, self.peaks[0].mz, "",
            len(self.members))


cdef class FeatureBase(object):
    def __init__(self, nodes):
        self.nodes = LCMSFeatureTreeList(nodes)

    cdef double get_mz(self):
        return 0.0

    cdef double get_start_time(self):
        return 0.0

    cdef double get_end_time(self):
        return INF


cdef class LCMSFeature(FeatureBase):
    def __init__(self, nodes=None, adducts=None, used_as_adduct=None):
        if nodes is None:
            nodes = []
        if adducts is None:
            adducts = []
        if used_as_adduct is None:
            used_as_adduct = []

        # self.nodes = LCMSFeatureTreeList(nodes)
        FeatureBase.__init__(self, nodes)
        self._total_intensity = -1
        self._mz = -1
        self._last_mz = 0.0
        self._most_abundant_member = 0.0
        self._retention_times = None
        self._peaks = None
        self._scan_ids = None
        self._start_time = -1
        self._end_time = -1
        self.adducts = adducts
        self.used_as_adduct = used_as_adduct

    def invalidate(self):
        self._invalidate()

    def _invalidate(self):
        self._total_intensity = -1
        self._last_mz = self._mz if self._mz != -1 else 0.
        self._mz = -1
        self._retention_times = None
        self._peaks = None
        self._scan_ids = None
        self._start_time = -1
        self._end_time = -1
        self._last_most_abundant_member = self._most_abundant_member
        self._most_abundant_member = None

    @property
    def most_abundant_member(self):
        if self._most_abundant_member is None:
            self._most_abundant_member = max(
                node.max_intensity() for node in self.nodes)
        return self._most_abundant_member

    def retain_most_abundant_member(self):
        self._mz = self._last_mz
        self._most_abundant_member = self._last_most_abundant_member

    @property
    def total_signal(self):
        if self._total_intensity == -1:
            total = 0.
            for node in self.nodes:
                total += node.total_intensity()
            self._total_intensity = total
        return self._total_intensity

    @property
    def intensity(self):
        return self.total_signal

    @property
    def neutral_mass(self):
        return self.mz

    cdef double get_mz(self):
        if self._mz == -1:
            best_mz = 0
            maximum_intensity = -1
            for node in self.nodes:
                intensity = node.max_intensity()
                if intensity > maximum_intensity:
                    maximum_intensity = intensity
                    best_mz = node.mz
            assert (len(self) == 0 or best_mz != 0)
            self._last_mz = self._mz = best_mz
        return self._mz

    @property
    def mz(self):
        if self._mz == -1:
            best_mz = 0
            maximum_intensity = -1
            for node in self.nodes:
                intensity = node.max_intensity()
                if intensity > maximum_intensity:
                    maximum_intensity = intensity
                    best_mz = node.mz
            assert (len(self) == 0 or best_mz != 0)
            self._last_mz = self._mz = best_mz
        return self._mz

    @property
    def retention_times(self):
        if self._retention_times is None:
            self._retention_times = np.array([node.retention_time for node in self.nodes])
        return self._retention_times

    @property
    def scan_ids(self):
        if self._scan_ids is None:
            self._scan_ids = tuple(node.scan_id for node in self.nodes)
        return self._scan_ids

    @property
    def peaks(self):
        if self._peaks is None:
            self._peaks = tuple(node.peaks for node in self.nodes)
        return self._peaks

    @property
    def start_time(self):
        if self._start_time == -1:
            self._start_time = self.nodes[0].retention_time
        return self._start_time

    @property
    def end_time(self):
        if self._end_time == -1:
            self._end_time = self.nodes[-1].retention_time
        return self._end_time

    cdef double get_start_time(self):
        if self._start_time == -1:
            self._start_time = self.nodes[0].retention_time
        return self._start_time

    cdef double get_end_time(self):
        if self._end_time == -1:
            self._end_time = self.nodes[-1].retention_time
        return self._end_time

    cpdef bint overlaps_in_time(self, FeatureBase interval):
        cdef:
            double self_start_time, self_end_time
            double other_start_time, other_end_time
        self_start_time = self.get_start_time()
        self_end_time = self.get_end_time()
        other_start_time = interval.get_start_time()
        other_end_time = interval.get_end_time()

        cond = ((self_start_time <= other_start_time and self_end_time >= other_end_time) or (
            self_start_time >= other_start_time and self_end_time <= other_end_time) or (
            self_start_time >= other_start_time and self_end_time >= other_end_time and
            self_start_time <= other_end_time) or (
            self_start_time <= other_start_time and self_end_time >= other_start_time) or (
            self_start_time <= other_end_time and self_end_time >= other_end_time))
        return cond

    def as_arrays(self):
        rts = np.array(
            [node.retention_time for node in self.nodes], dtype=np.float64)
        signal = np.array([node.total_intensity()
                           for node in self.nodes], dtype=np.float64)
        return rts, signal

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        return "%s(%0.4f, %0.2f, %0.2f)" % (
            self.__class__.__name__, self.mz, self.start_time, self.end_time)

    def split_sparse(self, delta_rt=1.):
        chunks = []
        current_chunk = []
        last_rt = self.nodes[0].retention_time

        for node in self.nodes:
            if (node.retention_time - last_rt) > delta_rt:
                x = self.__class__(LCMSFeatureTreeList(current_chunk))
                x.used_as_adduct = list(self.used_as_adduct)

                chunks.append(x)
                current_chunk = []

            last_rt = node.retention_time
            current_chunk.append(node)

        x = self.__class__(LCMSFeatureTreeList(current_chunk))
        x.used_as_adduct = list(self.used_as_adduct)

        chunks.append(x)
        for chunk in chunks:
            chunk.created_at = self.created_at

        for member in chunks:
            for other in chunks:
                if member == other:
                    continue
                assert not member.overlaps_in_time(other)

        return chunks

    def truncate_before(self, time):
        _, i = self.nodes.find_time(time)
        if self.nodes[i].retention_time < time:
            i += 1
        self.nodes = LCMSFeatureTreeList(self.nodes[i:])
        self._invalidate()

    def truncate_after(self, time):
        _, i = self.nodes.find_time(time)
        if self.nodes[i].retention_time < time:
            i += 1
        self.nodes = LCMSFeatureTreeList(self.nodes[:i])
        self._invalidate()

    def clone(self, cls=None):
        if cls is None:
            cls = self.__class__
        c = cls(self.nodes.clone(), list(self.adducts), list(self.used_as_adduct))
        c.created_at = self.created_at
        return c

    def insert_node(self, node):
        self.nodes.insert_node(node)
        self._invalidate()

    def insert(self, peak, retention_time):
        self.nodes.insert(retention_time, [peak])
        self._invalidate()

    def merge(self, other):
        new = self.clone()
        for node in other.nodes:
            node = node.clone()
            new.insert_node(node)
        new.created_at = "merge"
        return new

    @property
    def apex_time(self):
        rt, intensity = self.as_arrays()
        return rt[np.argmax(intensity)]

    def __iter__(self):
        return iter(self.nodes)

    cpdef bint _eq(self, other):
        if other is None:
            return False
        if len(self) != len(other):
            return False
        else:
            if abs(self.start_time - other.start_time) > 1e-4:
                return False
            elif abs(self.end_time - other.end_time) > 1e-4:
                return False
            else:
                for a, b in zip(self, other):
                    if a != b:
                        return False
                return True

    cpdef bint _ne(self, other):
        return not self._eq(other)

    def __richcmp__(self, other, int code):
        if code == 2:
            return self._eq(other)
        elif code == 3:
            return self._ne(other)
        else:
            return NotImplemented

    def __hash__(self):
        return hash((self.mz, self.start_time, self.end_time))

    def __getitem__(self, i):
        return self.nodes[i]


cdef class EmptyFeature(FeatureBase):
    def __init__(self, mz):
        FeatureBase.__init__(self, [])
        self._mz = mz
        self._start_time = 0
        self._end_time = INF

    @property
    def mz(self):
        return self._mz

    cdef double get_mz(self):
        return self._mz

    def __len__(self):
        return 0

    @property
    def start_time(self):
        return self._start_time

    @property
    def end_time(self):
        return self._end_time


def binsearch(array, value):
    lo = 0
    hi = len(array)
    while hi != lo:
        mid = (hi + lo) / 2
        point = array[mid]
        if value == point:
            return mid
        elif hi - lo == 1:
            return mid
        elif point > value:
            hi = mid
        else:
            lo = mid


cdef class FeatureSetIterator(object):
    @staticmethod
    def extract_time_points(features, start_time, end_time):
        time_points = set()
        for feature in features:
            rts = (feature.retention_times)
            time_points.update(rts)
        time_points = sorted(time_points)
        si = binsearch(time_points, start_time)
        ei = binsearch(time_points, end_time)
        time_points = time_points[si:ei + 1]
        return time_points

    def __init__(self, features, start_time=None, end_time=None):
        self.features = list(features)
        self.real_features = [f for f in features if f is not None and not isinstance(f, EmptyFeature)]
        self.start_time = max(
            [f.start_time for f in self.real_features]
        ) if start_time is None else start_time
        self.end_time = min(
            [f.end_time for f in self.real_features]
        ) if end_time is None else end_time
        self.time_point_queue = deque([self.start_time])
        self.last_time_seen = -1

    cpdef double get_next_time(self):
        cdef:
            set options
            list options_uniq
            object next_time
            double node_time
            tuple find_node
            object index
            size_t i, n
            LCMSFeature feature
            LCMSFeatureTreeNode node

        if self.time_point_queue:
            next_time = self.time_point_queue.popleft()
            return next_time
        else:
            options = set()
            n = len(self.real_features)
            for i in range(n):
                feature = <LCMSFeature>PyList_GET_ITEM(self.real_features, i)
                find_node = feature.nodes.find_time(self.last_time_seen)
                index = find_node[1]
                try:
                    node = <LCMSFeatureTreeNode>feature[index + 1]
                    node_time = node.retention_time
                    if node_time >= self.end_time:
                        continue
                    options.add(node_time)
                except IndexError:
                    continue
            if len(options) == 0:
                return -1
            options_uniq = sorted(options)
            next_time = options_uniq[0]
            self.time_point_queue.extend(
                <object>PyList_GetSlice(options_uniq, 1, PyList_GET_SIZE(options_uniq)))
            return next_time

    @property
    def current_time(self):
        return self.last_time_seen

    cpdef list get_peaks_for_time(self, double time):
        cdef:
            list peaks
            FeatureBase feature
            size_t i, n
            object ix
            tuple out
            LCMSFeatureTreeNode node
        peaks = []
        n = PyList_GET_SIZE(self.features)
        for i in range(n):
            feature = <FeatureBase>PyList_GET_ITEM(self.features, i)
            if feature is not None:
                try:
                    out = feature.nodes.find_time(time)
                    node = <LCMSFeatureTreeNode>out[0]
                except ValueError:
                    if not isinstance(feature, EmptyFeature):
                        raise
                    node = None
                if node is not None:
                    peaks.append(node.members[0])
                else:
                    peaks.append(None)
            else:
                peaks.append(None)
        return peaks

    def __next__(self):
        cdef:
            double time
            list peaks
        time = self.get_next_time()
        if time < 0:
            raise StopIteration()
        self.last_time_seen = time
        peaks = self.get_peaks_for_time(time)
        return peaks

    def __iter__(self):
        return self
