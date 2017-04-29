# cython: embedsignature=True

cimport cython

from ms_deisotope.utils import uid as _uid

from collections import deque
from random import sample as _sample

from libc.stdlib cimport malloc, free
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_GetSlice, PyList_GetItem
from ms_peak_picker._c.peak_set cimport FittedPeak, PeakBase
from brainpy._c.isotopic_distribution cimport TheoreticalPeak

cimport numpy as np
import numpy as np


np.import_array()


cdef object _zeros = np.zeros

cdef double INF = float("inf")


cdef class LCMSFeatureTreeList(object):  
    def __init__(self, roots=None):
        if roots is None:
            roots = []
        self.roots = list(roots)
        self._node_id_hash = None

    @staticmethod
    cdef LCMSFeatureTreeList _create(list roots):
        cdef:
            LCMSFeatureTreeList inst
        inst = LCMSFeatureTreeList.__new__(LCMSFeatureTreeList)
        inst.roots = list(roots)
        inst._node_id_hash = None
        return inst

    cdef void _invalidate(self):
        self._node_id_hash = None

    cpdef tuple find_time(self, double time):
        cdef:
            size_t indexout
            LCMSFeatureTreeNode node
        node = self._find_time(time, &indexout)
        return node, indexout

    cdef LCMSFeatureTreeNode _find_time(self, double time, size_t* indexout):
        cdef:
            int lo, hi
            double rt
            LCMSFeatureTreeNode node
        hi = PyList_GET_SIZE(self.roots)
        if hi == 0:
            raise ValueError()
        lo = 0
        while lo != hi:
            i = (lo + hi) / 2
            node = <LCMSFeatureTreeNode>PyList_GET_ITEM(self.roots, i)
            rt = node.time
            if rt == time:
                indexout[0] = i
                return node
            elif (hi - lo) == 1:
                indexout[0] = i
                return None
            elif rt < time:
                lo = i
            elif rt > time:
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

    def insert_node(self, LCMSFeatureTreeNode node):
        self._invalidate()
        cdef:
            LCMSFeatureTreeNode root
            size_t i
        try:
            root = self._find_time(node.time, &i)
            if root is None:
                if i != 0:
                    self.roots.insert(i + 1, node)
                else:
                    slot = self.getitem(i)
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

    cdef LCMSFeatureTreeNode getitem(self, size_t i):
        return <LCMSFeatureTreeNode>PyList_GetItem(self.roots, i)

    def __len__(self):
        return PyList_GET_SIZE(self.roots)

    cdef size_t get_size(self):
        return PyList_GET_SIZE(self.roots)

    def __iter__(self):
        return iter(self.roots)

    def __reduce__(self):
        return self.__class__, (self.roots,)

    def clone(self, deep=False):
        return LCMSFeatureTreeList(node.clone(deep=deep) for node in self)

    def unspool(self):
        cdef:
            list out_queue, stack
            LCMSFeatureTreeNode node, root
            size_t i, n
        out_queue = []
        n = self.get_size()
        for i in range(n):
            root = self.getitem(i)
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


cdef object uid = _uid
cdef object sample = _sample


cdef class LCMSFeatureTreeNode(object):
    def __init__(self, time, members=None):
        if members is None:
            members = []
        self.time = time
        self.members = members
        self._most_abundant_member = None
        self._mz = 0
        self._recalculate()
        self.node_id = uid()

    def clone(self, deep=False):
        cdef LCMSFeatureTreeNode node
        cdef list peaks
        cdef size_t i, n
        cdef PeakBase peak
        if deep:
            peaks = []
            n = self.get_members_size()
            for i in range(n):
                peak = self.getitem(i)
                if peak is not None:
                    peak = peak.clone()
                peaks.append(peak)
        else:
            peaks = list(self.members)
        node = self.__class__(self.time, peaks)
        node.node_id = self.node_id
        return node

    def _unspool_strip_children(self):
        node = self.__class__(
            self.time, list(self.members))
        yield node

    cpdef _recalculate(self):
        self._calculate_most_abundant_member()
        self._mz = self._most_abundant_member.mz

    cdef PeakBase _find_most_abundant_member(self):
        cdef:
            size_t i, n
            PeakBase peak, best
            double best_intensity
        n = self.get_members_size()
        best = self.getitem(0)
        best_intensity = best.intensity
        for i in range(n):
            peak = self.getitem(i)
            if peak.intensity > best_intensity:
                best = peak
                best_intensity = best.intensity
        return best

    cpdef _calculate_most_abundant_member(self):
        cdef:
            size_t n
        n = self.get_members_size()
        if n == 1:
            self._most_abundant_member = self.getitem(0)
        else:
            if n == 0:
                self._most_abundant_member = None
            else:
                self._most_abundant_member = self._find_most_abundant_member()
        return

    def __getstate__(self):
        return (self.time, self.members, self.node_id)

    def __setstate__(self, state):
        self.time, self.members, self.node_id = state
        self._recalculate()

    def __reduce__(self):
        return self.__class__, (self.time, []), self.__getstate__()

    cdef inline size_t get_members_size(self):
        return PyList_GET_SIZE(self.members)

    @property
    def mz(self):
        if self._mz == 0:
            if self._most_abundant_member is not None:
                self._mz = self._most_abundant_member.mz
        return self._mz

    cdef double get_mz(self):
        if self._mz == 0:
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
        cdef:
            double total
            size_t n, i
            PeakBase peak
        total = 0.
        n = self.get_members_size()
        for i in range(n):
            peak = self.getitem(i)
            total += peak.intensity
        return total

    cpdef double max_intensity(self):
        return self._most_abundant_member.intensity

    cpdef double total_intensity(self):
        return self._total_intensity_members()

    cpdef bint _eq(self, LCMSFeatureTreeNode other):
        return self.members == other.members and abs(
            self.time - other.time) < 1e-4

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
            self.time, self.peaks[0].mz, "",
            len(self.members))

    cdef inline FittedPeak getpeak(self, size_t i):
        return <FittedPeak>PyList_GET_ITEM(self.members, i)

    cdef PeakBase getitem(self, size_t i):
        return <PeakBase>PyList_GET_ITEM(self.members, i)


@cython.freelist(10000000)
cdef class FeatureBase(object):
    def __init__(self, nodes):
        self.nodes = LCMSFeatureTreeList(nodes)

    cdef double get_mz(self):
        return 0.0

    cdef double get_start_time(self):
        return 0.0

    cdef double get_end_time(self):
        return INF

    cpdef tuple find_time(self, double time):
        cdef:
            size_t indexout
            LCMSFeatureTreeNode node
        node = self._find_time(time, &indexout)
        return node, indexout

    cdef inline LCMSFeatureTreeNode _find_time(self, double time, size_t* indexout):
        return self.nodes._find_time(time, indexout)

    cdef inline size_t get_size(self):
        return self.nodes.get_size()

    def __len__(self):
        return self.get_size()

    def __iter__(self):
        return iter(self.nodes)


cdef class LCMSFeature(FeatureBase):
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
        self._total_intensity = -1
        self._mz = -1
        self._last_mz = 0.0
        self._times = None
        self._peaks = None
        self._start_time = -1
        self._end_time = -1
        self.adducts = adducts
        self.used_as_adduct = used_as_adduct
        self.feature_id = feature_id
        self._peak_averager = RunningWeightedAverage._create(None)
        if self.get_size() > 0:
            self._feed_peak_averager()

    cdef void _feed_peak_averager(self):
        cdef:
            size_t i, n
            LCMSFeatureTreeNode node
        n = self.get_size()
        for i in range(n):
            node = self.getitem(i)
            self._peak_averager._update(node.members)

    def invalidate(self, reaverage=False):
        self._invalidate(reaverage)

    cpdef _invalidate(self, bint reaverage=False):
        self._total_intensity = -1
        self._last_mz = self._mz if self._mz != -1 else 0.
        self._mz = -1
        self._times = None
        self._peaks = None
        self._start_time = -1
        self._end_time = -1
        if reaverage:
            self._peak_averager = RunningWeightedAverage._create(None)
            self._feed_peak_averager()

    @property
    def total_signal(self):
        cdef:
            double total
            size_t i, n
            LCMSFeatureTreeNode node
        if self._total_intensity == -1:
            total = 0.
            n = self.get_size()
            for i in range(n):
                node = self.getitem(i)
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
        cdef:
            size_t i, n
            double best_mz, temp
            LCMSFeatureTreeNode node
        if self._mz == -1:
            temp = best_mz = self._peak_averager.current_mean
            if best_mz < 0:
                raise ValueError("best_mz < 0")
            self._last_mz = self._mz = best_mz
        return self._mz

    @property
    def mz(self):
        return self.get_mz()

    @property
    def times(self):
        cdef:
            size_t i, n
            LCMSFeatureTreeNode node
            np.ndarray[double, ndim=1, mode='c'] acc
        if self._times is None:
            n = self.get_size()
            acc = _zeros(n)
            for i in range(n):
                node = self.getitem(i)
                acc[i] = node.time
            self._times = acc
        return self._times

    @property
    def peaks(self):
        if self._peaks is None:
            self._peaks = tuple(node.peaks for node in self.nodes)
        return self._peaks

    @property
    def start_time(self):
        if self._start_time == -1:
            self._start_time = self.nodes[0].time
        return self._start_time

    @property
    def end_time(self):
        if self._end_time == -1:
            self._end_time = self.nodes[-1].time
        return self._end_time

    cdef double get_start_time(self):
        if self._start_time == -1:
            self._start_time = self.getitem(0).time
        return self._start_time

    cdef double get_end_time(self):
        cdef:
            size_t n
        if self._end_time == -1:
            n = self.get_size()
            self._end_time = self.getitem(n - 1).time
        return self._end_time

    cpdef bint overlaps_in_time(self, FeatureBase interval):
        cdef:
            double self_start_time, self_end_time
            double other_start_time, other_end_time
            bint cond

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
            [node.time for node in self.nodes], dtype=np.float64)
        signal = np.array([node.total_intensity()
                           for node in self.nodes], dtype=np.float64)
        return rts, signal

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
        c = cls(self.nodes.clone(deep=deep), list(self.adducts), list(self.used_as_adduct))
        c.feature_id = self.feature_id
        c.created_at = self.created_at
        return c

    cpdef insert_node(self, LCMSFeatureTreeNode node):
        self._peak_averager._update(node.members)
        self.nodes.insert_node(node)
        self._invalidate()

    cpdef insert(self, PeakBase peak, double time):
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

    def __reduce__(self):
        return self.__class__, (self.nodes, self.adducts, self.used_as_adduct, self.feature_id)

    def __hash__(self):
        return hash((self.mz, self.start_time, self.end_time))

    def __getitem__(self, i):
        return self.nodes[i]

    cdef LCMSFeatureTreeNode getitem(self, size_t i):
        return self.nodes.getitem(i)


cdef class EmptyFeature(FeatureBase):
    def __init__(self, mz):
        FeatureBase.__init__(self, [])
        self._mz = mz
        self._start_time = 0
        self._end_time = INF

    @staticmethod
    cdef EmptyFeature _create(double mz):
        cdef EmptyFeature inst = EmptyFeature.__new__(EmptyFeature)
        inst._mz = mz
        inst.nodes = LCMSFeatureTreeList._create([])
        inst._start_time = 0
        inst._end_time = INF
        return inst

    @property
    def mz(self):
        return self._mz

    cdef double get_mz(self):
        return self._mz

    def __len__(self):
        return 0

    def __getitem__(self, i):
        raise IndexError()

    @property
    def start_time(self):
        return self._start_time

    @property
    def end_time(self):
        return self._end_time


cdef class FeatureSetIterator(object):
    def __init__(self, features, start_time=None, end_time=None):
        cdef:
            size_t n, i

        self.features = list(features)
        n = PyList_GET_SIZE(self.features)
        self.real_features = []
        self.index_list = <size_t*>malloc(sizeof(size_t) * n)
        self.start_time = 0
        self.end_time = INF

        if start_time is not None:
            self.start_time = start_time
        if end_time is not None:
            self.end_time = end_time

        for i in range(n):
            feature = <FeatureBase>PyList_GET_ITEM(self.features, i)
            if feature is None or isinstance(feature, EmptyFeature):
                continue
            self.real_features.append(feature)
            time = feature.get_start_time()
            if (self.start_time < time) and (start_time is None):
                self.start_time = time
            time = feature.get_end_time()
            if (self.end_time > time) and (end_time is None):
                self.end_time = time
            self.index_list[i] = 0

        self.init_indices()
        self.last_time_seen = -1
    
    cdef void _initialize(self, list features):
        cdef:
            size_t i, n
        self.features = features
        n = PyList_GET_SIZE(self.features)
        self.real_features = []
        self.index_list = <size_t*>malloc(sizeof(size_t) * n)

        self.start_time = 0
        self.end_time = INF

        for i in range(n):
            feature = <FeatureBase>PyList_GET_ITEM(self.features, i)
            if feature is None or isinstance(feature, EmptyFeature):
                continue
            self.real_features.append(feature)
            time = feature.get_start_time()
            if (self.start_time < time):
                self.start_time = time
            time = feature.get_end_time()
            if (self.end_time > time):
                self.end_time = time
            self.index_list[i] = 0

        self.init_indices()
        self.last_time_seen = -1

    @staticmethod
    cdef FeatureSetIterator _create(list features):
        cdef FeatureSetIterator self
        self = FeatureSetIterator.__new__(FeatureSetIterator)
        self._initialize(features)
        return self

    @staticmethod
    cdef FeatureSetIterator _create_with_threshold(list features, list theoretical_distribution, double detection_threshold):
        cdef:
            FeatureSetIterator self
            TheoreticalPeak p
            size_t i
            bint passed_threshold
            double start_time, end_time
        self = FeatureSetIterator.__new__(FeatureSetIterator)
        self._initialize(features)

        start_time = -INF
        end_time = INF
        for i in range(self.get_size()):
            feature = self.getitem(i)
            if feature is None or isinstance(feature, EmptyFeature):
                continue
            p = <TheoreticalPeak>PyList_GET_ITEM(theoretical_distribution, i)
            
            passed_threshold = p.intensity > detection_threshold

            time = feature.get_start_time()
            if (start_time < time) and passed_threshold:
                start_time = time
            time = feature.get_end_time()
            if (end_time > time) and passed_threshold:
                end_time = time
        if end_time != INF:
            self.start_time = start_time
            self.end_time = end_time

        self.init_indices()
        self.last_time_seen = -1
        return self

    def __dealloc__(self):
        free(self.index_list)

    cpdef init_indices(self):
        cdef:
            size_t i, ix, n
            FeatureBase f
            LCMSFeatureTreeNode node
        n = self.get_size()
        for i in range(n):
            f = self.getitem(i)
            if f is None:
                self.index_list[i] = 0
                continue
            if f.get_size() > 0:
                node = f._find_time(self.start_time, &ix)
            else:
                ix = 0
            self.index_list[i] = ix
    
    cpdef double get_next_time(self):
        cdef:
            double time, ix_time
            size_t i, n, ix
            FeatureBase f
        time = INF
        n = self.get_size()
        for i in range(n):
            f = self.getitem(i)
            if f is None:
                continue
            ix = self.index_list[i]
            if ix < f.get_size():
                ix_time = f.nodes.getitem(ix).time
                if ix_time < time and ix_time < self.end_time:
                    time = ix_time
        if time == INF:
            return -1
        else:
            return time
    
    cpdef bint has_more(self):
        cdef:
            size_t j, i, n, ix
            FeatureBase f
            bint done, at_end_time
        j = 0
        n = self.get_size()
        for i in range(n):
            f = self.getitem(i)
            if f is None:
                j += 1
                continue
            ix = self.index_list[i]
            done = ix >= f.get_size()
            if not done:
                at_end_time = f.nodes.getitem(ix).time >= self.end_time
            else:
                at_end_time = done
            if done or at_end_time:
                j += 1
        return j != n
    
    cdef inline size_t get_size(self):
        return PyList_GET_SIZE(self.features)

    cdef inline FeatureBase getitem(self, size_t i):
        return <FeatureBase>PyList_GET_ITEM(self.features, i)

    def get_index_list(self):
        out = []
        for i in range(self.get_size()):
            out.append(self.index_list[i])
        return out

    cpdef list get_peaks_for_time(self, double time):
        cdef:
            list peaks
            LCMSFeatureTreeNode node
            FeatureBase feature
            size_t i, n, ix
        peaks = []
        i = 0
        n = self.get_size()
        for i in range(n):
            feature = self.getitem(i)
            if feature is not None:
                if feature.get_size() > 0:
                    node = feature._find_time(time, &ix)
                    if ix == self.index_list[i] and self.index_list[i] == 0 and node is None:
                        pass
                    elif ix >= self.index_list[i]:
                        self.index_list[i] += 1
                else:
                    if not isinstance(feature, EmptyFeature):
                        raise ValueError("Empty Feature is not of type EmptyFeature")
                    node = None
                if node is not None:
                    peaks.append(
                        node.getpeak(0))
                else:
                    peaks.append(None)
            else:
                peaks.append(None)
        return peaks

    @property
    def current_time(self):
        return self.last_time_seen
    
    cpdef double get_current_time(self):
        return self.last_time_seen

    cpdef list get_next_value(self):
        cdef:
            double time
            list peaks
        time = self.get_next_time()
        if time < 0:
            return None
        self.last_time_seen = time
        peaks = self.get_peaks_for_time(time)
        return peaks

    def __next__(self):
        cdef:
            list peaks
        peaks = self.get_next_value()
        if peaks is None:
            raise StopIteration()
        return peaks
    
    def __iter__(self):
        return self

    def __repr__(self):
        return "FeatureSetIterator(%d/%d, %0.3f-%0.3f)" % (
            len(self.real_features), len(self.features), self.start_time, self.end_time)


cdef class RunningWeightedAverage(object):

    def __cinit__(self, *args, **kwargs):
        self.accumulator = []
        self.current_mean = 0
        self.total_weight = 0
        self.current_count = 0

    def __init__(self, iterable=None):
        if iterable is not None:
            self.update(iterable)

    @staticmethod
    cdef RunningWeightedAverage _create(list peaks):
        cdef:
            RunningWeightedAverage inst

        inst = RunningWeightedAverage.__new__(RunningWeightedAverage)
        if peaks is None:
            return inst
        else:
            inst._update(peaks)
            return inst

    @cython.nonecheck(False)
    @cython.cdivision(True)
    cpdef add(self, PeakBase peak):
        if peak.intensity == 0:
            if self.current_mean == 0 and self.total_weight == 0:
                self.current_mean = peak.mz
                self.total_weight = 1
            else:
                return

        self.accumulator.append(peak)
        agg = (self.total_weight * self.current_mean) + \
            (peak.mz * peak.intensity)
        self.total_weight += peak.intensity
        self.current_count += 1

        if self.total_weight != 0:
            self.current_mean = agg / (self.total_weight)
        else:
            print("NaN produced in add()")
        return self

    cpdef double recompute(self):
        cdef:
            size_t i
            PeakBase peak
            double weight, total
        weight = 0
        total = 0
        for i in range(self.current_count):
            peak = <PeakBase>PyList_GET_ITEM(self.accumulator, i)
            weight += peak.intensity
            total += peak.intensity * peak.mz
        return total / weight

    cpdef RunningWeightedAverage update(self, iterable):
        for x in iterable:
            self.add(x)
        return self

    cdef void _update(self, list peaks):
        cdef:
            size_t i, n
            PeakBase p
        n = PyList_GET_SIZE(peaks)
        for i in range(n):
            p = <PeakBase>PyList_GET_ITEM(peaks, i)
            self.add(p)

    def __repr__(self):
        return "RunningWeightedAverage(%r, %d)" % (self.current_mean, self.current_count)

    @cython.nonecheck(False)
    @cython.cdivision(True)
    cpdef double _bootstrap(self, size_t n=150, size_t k=40):
        cdef:
            size_t i
            size_t tn, ti
            list traces, sampled_peaks
            RunningWeightedAverage inst
            double total_weight, accumulation, center
        traces = []
        if self.current_count < k:
            k = self.current_count
        for i in range(n):
            sampled_peaks = sample(self.accumulator, k)
            inst = RunningWeightedAverage._create(sampled_peaks)
            traces.append(inst)
        tn = i + 1
        total_weight = 0
        accumulation = 0
        for ti in range(tn):
            inst = traces[ti]
            accumulation += inst.total_weight * inst.current_mean
            total_weight += inst.total_weight
        center = accumulation / total_weight
        return center

    cpdef RunningWeightedAverage bootstrap(self, size_t n=150, size_t k=40):
        self.current_mean = self._bootstrap(n, k)
        return self
