import random
import numpy as np

from collections import OrderedDict

from ms_deisotope.utils import uid


class FeatureBase(object):
    def __init__(self, nodes):
        self.nodes = LCMSFeatureTreeList(nodes)

    def find_time(self, time_point):
        return self.nodes.find_time(time_point)


class LCMSFeature(FeatureBase):
    created_at = "new"

    def __init__(self, nodes=None, adducts=None, used_as_adduct=None, feature_id=None):
        if nodes is None:
            nodes = []
        if adducts is None:
            adducts = []
        if used_as_adduct is None:
            used_as_adduct = []
        if feature_id is None:
            feature_id = uid()

        FeatureBase.__init__(self, nodes)
        self._total_intensity = None
        self._mz = None
        self._last_mz = 0.0
        self._times = None
        self._peaks = None
        self._start_time = None
        self._end_time = None
        self.adducts = adducts
        self.used_as_adduct = used_as_adduct
        self.feature_id = feature_id
        self._peak_averager = RunningWeightedAverage()
        if len(self) > 0:
            self._feed_peak_averager()

    def _feed_peak_averager(self):
        for node in self:
            self._peak_averager.update(node.members)

    def invalidate(self, reaverage=False):
        self._invalidate(reaverage)

    def _invalidate(self, reaverage=False):
        self._total_intensity = None
        self._last_mz = self._mz if self._mz is not None else 0.
        self._mz = None
        self._times = None
        self._peaks = None
        self._start_time = None
        self._end_time = None

        if reaverage:
            self._peak_averager = RunningWeightedAverage()
            if len(self) > 0:
                self._feed_peak_averager()

    @property
    def total_signal(self):
        if self._total_intensity is None:
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

    @property
    def mz(self):
        if self._mz is None:
            best_mz = self._peak_averager.current_mean
            self._last_mz = self._mz = best_mz
        return self._mz

    @property
    def times(self):
        if self._times is None:
            self._times = np.array(
                [node.time for node in self.nodes])
        return self._times

    @property
    def peaks(self):
        if self._peaks is None:
            self._peaks = tuple(node.peaks for node in self.nodes)
        return self._peaks

    @property
    def start_time(self):
        if self._start_time is None:
            self._start_time = self.nodes[0].time
        return self._start_time

    @property
    def end_time(self):
        if self._end_time is None:
            self._end_time = self.nodes[-1].time
        return self._end_time

    def overlaps_in_time(self, interval):
        cond = ((self.start_time <= interval.start_time and self.end_time >= interval.end_time) or (
            self.start_time >= interval.start_time and self.end_time <= interval.end_time) or (
            self.start_time >= interval.start_time and self.end_time >= interval.end_time and
            self.start_time <= interval.end_time) or (
            self.start_time <= interval.start_time and self.end_time >= interval.start_time) or (
            self.start_time <= interval.end_time and self.end_time >= interval.end_time))
        return cond

    def spans_in_time(self, time):
        return self.start_time <= time <= self.end_time

    def as_arrays(self):
        rts = np.array(
            [node.time for node in self.nodes], dtype=np.float64)
        signal = np.array([node.total_intensity()
                           for node in self.nodes], dtype=np.float64)
        return rts, signal

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        return "%s(%0.4f, %0.2f, %0.2f)" % (
            self.__class__.__name__, self.mz, self.start_time, self.end_time)

    def _copy_chunk(self, nodes, *args, **kwargs):
        x = self.__class__(LCMSFeatureTreeList(nodes))
        x.used_as_adduct = list(self.used_as_adduct)
        return x

    def split_at(self, time):
        _, i = self.nodes.find_time(time)
        if self.nodes[i].time < time:
            i += 1
        return LCMSFeature(self.nodes[:i]), LCMSFeature(self.nodes[i:])

    def split_sparse(self, delta_rt=1.):
        chunks = []
        current_chunk = []
        last_rt = self.nodes[0].time

        for node in self.nodes:
            if (node.time - last_rt) > delta_rt:
                x = self._copy_chunk(current_chunk)

                chunks.append(x)
                current_chunk = []

            last_rt = node.time
            current_chunk.append(node)

        x = self._copy_chunk(current_chunk)

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
        if self.nodes[i].time < time:
            i += 1
        self.nodes = LCMSFeatureTreeList(self.nodes[i:])
        self._invalidate()

    def truncate_after(self, time):
        _, i = self.nodes.find_time(time)
        if self.nodes[i].time < time:
            i += 1
        self.nodes = LCMSFeatureTreeList(self.nodes[:i])
        self._invalidate()

    def clone(self, deep=False, cls=None):
        if cls is None:
            cls = self.__class__
        c = cls(self.nodes.clone(deep=deep), list(
            self.adducts), list(self.used_as_adduct))
        c.feature_id = self.feature_id
        c.created_at = self.created_at
        return c

    def insert_node(self, node):
        self._peak_averager.update(node.members)
        self.nodes.insert_node(node)
        self._invalidate()

    def insert(self, peak, time):
        self._peak_averager.add(peak)
        self.nodes.insert(time, [peak])
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

    def __eq__(self, other):
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

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.mz, self.start_time, self.end_time))

    def __getitem__(self, i):
        return self.nodes[i]


class LCMSFeatureTreeList(object):

    def __init__(self, roots=None):
        if roots is None:
            roots = []
        self.roots = list(roots)
        self._node_id_hash = None

    def _invalidate(self):
        self._node_id_hash = None

    def find_time(self, time):
        if len(self.roots) == 0:
            raise ValueError()
        lo = 0
        hi = len(self.roots)
        while lo != hi:
            i = (lo + hi) // 2
            node = self.roots[i]
            if node.time == time:
                return node, i
            elif (hi - lo) == 1:
                return None, i
            elif node.time < time:
                lo = i
            elif node.time > time:
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
            root, i = self.find_time(node.time)
            if root is None:
                if i != 0:
                    self.roots.insert(i + 1, node)
                else:
                    slot = self.roots[i]
                    if slot.time < node.time:
                        i += 1
                    self.roots.insert(i, node)
            else:
                root.add(node)
            return i
        except ValueError:
            self.roots.append(node)
            return 0

    def insert(self, time, peaks):
        node = LCMSFeatureTreeNode(time, peaks)
        return self.insert_node(node)

    def extend(self, iterable):
        for peaks, time in iterable:
            self.insert(time, peaks)

    def __getitem__(self, i):
        return self.roots[i]

    def __len__(self):
        return len(self.roots)

    def __iter__(self):
        return iter(self.roots)

    def clone(self, deep=False):
        return LCMSFeatureTreeList(node.clone(deep=deep) for node in self)

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
            len(self), self[0].time, self[-1].time)


class LCMSFeatureTreeNode(object):

    __slots__ = ["time", "members",
                 "_most_abundant_member", "_mz", "node_id"]

    def __init__(self, time=None, members=None):
        if members is None:
            members = []
        self.time = time
        self.members = members
        self._most_abundant_member = None
        self._mz = 0
        self._recalculate()
        self.node_id = uid()

    def clone(self, deep=False):
        if deep:
            peaks = [peak.clone() for peak in self.members]
        else:
            peaks = list(self.members)
        node = self.__class__(
            self.time, peaks)
        node.node_id = self.node_id
        return node

    def _unspool_strip_children(self):
        node = self.__class__(
            self.time, list(self.members))
        yield node

    def _calculate_most_abundant_member(self):
        if len(self.members) == 1:
            self._most_abundant_member = self.members[0]
        else:
            if len(self.members) == 0:
                self._most_abundant_member = None
            else:
                self._most_abundant_member = max(
                    self.members, key=lambda x: x.intensity)

    def _recalculate(self):
        self._calculate_most_abundant_member()
        self._mz = self._most_abundant_member.mz

    def __getstate__(self):
        return (self.time, self.members, self.node_id)

    def __setstate__(self, state):
        self.time, self.members, self.node_id = state
        self._recalculate()

    def __reduce__(self):
        return self.__class__, (self.time, []), self.__getstate__()

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

    def _total_intensity_members(self):
        total = 0.
        for peak in self.members:
            total += peak.intensity
        return total

    def max_intensity(self):
        return self._most_abundant_member.intensity

    def total_intensity(self):
        return self._total_intensity_members()

    def __eq__(self, other):
        return self.members == other.members and abs(
            self.time - other.time) < 1e-4

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash(self.uid)

    @property
    def peaks(self):
        peaks = list(self.members)
        return peaks

    def __repr__(self):
        return "%s(%f, %0.4f %s|%d)" % (
            self.__class__.__name__,
            self.time, self.peaks[0].mz, "",
            len(self.members))


class EmptyFeature(FeatureBase):

    def __init__(self, mz):
        FeatureBase.__init__(self, [])
        self.mz = mz
        self.start_time = 0
        self.end_time = float('inf')

    def __len__(self):
        return 0

    @property
    def times(self):
        return np.array([])

    def __iter__(self):
        return iter([])

    def __repr__(self):
        return "EmptyFeature(%f)" % (self.mz,)


class FeatureSetIterator(object):

    def __init__(self, features, start_time=None, end_time=None):
        self.features = features
        self.real_features = [
            f for f in features if f is not None and not isinstance(f, EmptyFeature)]
        self.start_time = max(
            [f.start_time for f in self.real_features]
        ) if start_time is None else start_time
        self.end_time = min(
            [f.end_time for f in self.real_features]
        ) if end_time is None else end_time

        self.index_list = [0 for f in self.features]
        self.init_indices()
        self.last_seen_time = -1

    def init_indices(self):
        for i, f in enumerate(self.features):
            if f is None:
                self.index_list[i] = 0
                continue
            try:
                node, ix = f.find_time(self.start_time)
            except ValueError:
                if not isinstance(f, EmptyFeature):
                    raise
                else:
                    ix = 0
            self.index_list[i] = ix

    def get_next_time(self):
        time = float('inf')
        for i, ix in enumerate(self.index_list):
            try:
                f = self.features[i]
                if f is None:
                    continue
                if ix < len(f):
                    ix_time = self.features[i][ix].time
                    if ix_time < time and ix_time < self.end_time:
                        time = ix_time
            except IndexError:
                continue
        if time == float('inf'):
            return None
        else:
            return time

    def has_more(self):
        j = 0
        for i, f in enumerate(self.features):
            if f is None:
                j += 1
                continue
            ix = self.index_list[i]
            done = ix >= len(f)
            if not done:
                at_end_time = f[ix].time >= self.end_time
            else:
                at_end_time = done
            if done or at_end_time:
                j += 1
        return j != len(self.features)

    def get_peaks_for_time(self, time):
        peaks = []
        i = 0
        for feature in self.features:
            if feature is not None:
                try:
                    node, ix = feature.nodes.find_time(time)
                    if ix == self.index_list[i] and self.index_list[i] == 0 and node is None:
                        pass
                    elif ix >= self.index_list[i]:
                        self.index_list[i] += 1
                except ValueError:
                    if not isinstance(feature, EmptyFeature):
                        raise
                    node = None
                if node is not None:
                    peaks.append(
                        node.members[0])
                else:
                    peaks.append(None)
            else:
                peaks.append(None)
            i += 1
        return peaks

    @property
    def current_time(self):
        return self.last_seen_time

    def get_next_value(self):
        time = self.get_next_time()
        if time is None:
            raise StopIteration()
        self.last_seen_time = time
        peaks = self.get_peaks_for_time(time)
        return peaks

    def __next__(self):
        return self.get_next_value()

    def next(self):
        return self.get_next_value()

    def __iter__(self):
        return self

    def __repr__(self):
        return "FeatureSetIterator(%d/%d, %0.3f-%0.3f)" % (
            len(self.real_features), len(self.features), self.start_time, self.end_time)


try:
    has_c = True
    _LCMSFeatureTreeList = LCMSFeatureTreeList
    _LCMSFeatureTreeNode = LCMSFeatureTreeNode
    _LCMSFeature = LCMSFeature
    _FeatureSetIterator = FeatureSetIterator
    _EmptyFeature = EmptyFeature
    _FeatureBase = FeatureBase
    from ms_deisotope._c.feature_map.lcms_feature import (
        LCMSFeatureTreeList,
        LCMSFeatureTreeNode,
        LCMSFeature,
        FeatureSetIterator,
        EmptyFeature,
        FeatureBase
    )
except ImportError:
    has_c = False


class NodeFeatureSetIterator(FeatureSetIterator):

    def get_peaks_for_time(self, time):
        peaks = []
        for feature in self.features:
            if feature is not None:
                try:
                    node, ix = feature.nodes.find_time(time)
                except ValueError:
                    if not isinstance(feature, EmptyFeature):
                        raise
                    node = None
                if node is not None:
                    peaks.append(node)
                else:
                    peaks.append(None)
            else:
                peaks.append(None)
        return peaks


def iterpeaks(feature):
    for node in feature:
        for peak in node.members:
            yield peak


class RunningWeightedAverage(object):

    def __init__(self, iterable=None):
        self.total_weight = 0
        self.current_mean = 0
        self.current_count = 0

        if iterable is not None:
            self.update(iterable)

    def add(self, peak):
        if peak.intensity == 0:
            if self.current_mean == 0 and self.total_weight == 0:
                self.current_mean = peak.mz
                self.total_weight = 1
            else:
                return
        agg = (self.total_weight * self.current_mean) + \
            (peak.mz * peak.intensity)
        self.total_weight += peak.intensity
        self.current_mean = agg / self.total_weight
        self.current_count += 1
        return self

    def update(self, iterable):
        for x in iterable:
            self.add(x)
        return self

    def __repr__(self):
        return "RunningWeightedAverage(%r, %d)" % (self.current_mean, self.current_count)


class SimpleLCMSFeature(object):
    def __init__(self, storage=None):
        if storage is None:
            storage = OrderedDict()
        self.storage = storage

    def as_arrays(self):
        return (
            np.array(self.storage.keys()),
            np.array(self.storage.values()))

    @property
    def total_signal(self):
        return sum(self.storage.values())

    @property
    def apex_time(self):
        rt, intensity = self.as_arrays()
        return rt[np.argmax(intensity)]

    @property
    def start_time(self):
        return tuple(self.storage.keys())[0]

    @property
    def end_time(self):
        return tuple(self.storage.keys())[-1]


def smooth_feature(feature, level=1):
    from ms_deisotope.feature_map.profile_transform import smooth_leveled
    time, intensity = feature.as_arrays()
    smoothed = smooth_leveled(time, intensity, level)
    for i, node in enumerate(feature):
        peaks = node.members
        if len(peaks) > 1:
            peaks.sort(key=lambda x: x.intensity, reverse=True)
            peak = peaks[0]
            for other in peaks[1:]:
                peak.intensity += other.intensity
            node.members = [peak]
        else:
            peak = node.members[0]
        peak.intensity = smoothed[i]
    feature.invalidate(True)
    return feature
