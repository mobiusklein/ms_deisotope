# cython: embedsignature=True

cimport cython

from ms_deisotope.utils import uid as _uid

from collections import deque
from random import sample as _sample

from libc.stdlib cimport malloc, free

from cpython cimport PyErr_SetString
from cpython.list cimport PyList_GET_SIZE, PyList_GET_ITEM, PyList_GetSlice, PyList_GetItem
from ms_peak_picker._c.peak_set cimport FittedPeak, PeakBase
from brainpy._c.isotopic_distribution cimport TheoreticalPeak

from ms_deisotope._c.peak_set cimport DeconvolutedPeak

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

    cpdef LCMSFeatureTreeNodeBase _make_node(self, double time, list peaks):
        return LCMSFeatureTreeNode(time, peaks)

    cdef void _invalidate(self):
        self._node_id_hash = None

    cpdef tuple find_time(self, double time):
        cdef:
            size_t indexout
            LCMSFeatureTreeNodeBase node
        node = self._find_time(time, &indexout)
        return node, indexout

    cdef LCMSFeatureTreeNodeBase _find_time(self, double time, size_t* indexout):
        cdef:
            int lo, hi
            double rt
            LCMSFeatureTreeNodeBase node
        hi = PyList_GET_SIZE(self.roots)
        if hi == 0:
            raise ValueError()
        lo = 0
        while lo != hi:
            i = (lo + hi) // 2
            node = <LCMSFeatureTreeNodeBase>PyList_GET_ITEM(self.roots, i)
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

    def insert_node(self, LCMSFeatureTreeNodeBase node):
        self._invalidate()
        cdef:
            LCMSFeatureTreeNodeBase root
            size_t i
        i = 0
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
        node = self._make_node(time, peaks)
        return self.insert_node(node)

    def extend(self, iterable):
        for peaks, time in iterable:
            self.insert(time, peaks)

    def __getitem__(self, i):
        return self.roots[i]

    cdef LCMSFeatureTreeNodeBase getitem(self, size_t i):
        return <LCMSFeatureTreeNodeBase>PyList_GetItem(self.roots, i)

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
            LCMSFeatureTreeNodeBase node, root
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


cdef class LCMSFeatureTreeNodeBase(object):

    def clone(self, deep=False):
        cdef LCMSFeatureTreeNodeBase node
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

    def __getstate__(self):
        return (self.time, self.members, self.node_id)

    def __setstate__(self, state):
        self.time, self.members, self.node_id = state
        self._recalculate()

    def __reduce__(self):
        return self.__class__, (self.time, []), self.__getstate__()

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

    cpdef _recalculate(self):
        self._calculate_most_abundant_member()

    cdef inline size_t get_members_size(self):
        return PyList_GET_SIZE(self.members)

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

    cdef PeakBase getitem(self, size_t i):
        return <PeakBase>PyList_GET_ITEM(self.members, i)

    @property
    def peaks(self):
        peaks = list(self.members)
        return peaks

    def __hash__(self):
        return hash(self.uid)


cdef class LCMSFeatureTreeNode(LCMSFeatureTreeNodeBase):
    def __init__(self, time, members=None):
        if members is None:
            members = []
        self.time = time
        self.members = members
        self._most_abundant_member = None
        self._mz = 0
        self._recalculate()
        self.node_id = uid()

    cpdef _recalculate(self):
        self._calculate_most_abundant_member()
        if self._most_abundant_member is not None:
            self._mz = self._most_abundant_member.mz

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

    def __repr__(self):
        return "%s(%f, %0.4f %s|%d)" % (
            self.__class__.__name__,
            self.time, self.peaks[0].mz, "",
            len(self.members))

    cdef inline FittedPeak getpeak(self, size_t i):
        return <FittedPeak>PyList_GET_ITEM(self.members, i)


@cython.freelist(10000000)
cdef class FeatureBase(object):
    def __init__(self, nodes):
        self.nodes = LCMSFeatureTreeList(nodes)

    cdef double get_mz(self):
        return 0.0

    cdef double get_neutral_mass(self):
        return self.get_mz()

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
        self._initialize_averager()

    cpdef _initialize_averager(self):
        self._peak_averager = RunningWeightedAverage._create(None)

    cpdef _update_from_averager(self):
        best_mz = self._peak_averager.current_mean
        if best_mz < 0:
            raise ValueError("best_mz < 0")
        self._last_mz = self._mz = best_mz

    cpdef _reaverage_and_update(self):
        self._peak_averager.reset()
        self._peak_averager.feed_from_feature(self)
        self._update_from_averager()

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
            self._reaverage_and_update()

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
        if self._mz == -1:
            self._reaverage_and_update()
        return self._mz

    @property
    def mz(self):
        return self.get_mz()

    @property
    def times(self):
        cdef:
            size_t i, n
            LCMSFeatureTreeNode node
            np.npy_intp knd
            np.ndarray[double, ndim=1, mode='c'] acc
        if self._times is None:
            n = self.get_size()
            knd = n
            acc = np.PyArray_ZEROS(1, &knd, np.NPY_FLOAT64, 0)
            # acc = _zeros(n)
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
        """Test whether or not this feature overlaps with `interval`.

        Parameters
        ----------
        interval : FeatureBase
            The feature to test spanning

        Returns
        -------
        bool:
            Whether this feature spans `interval` or the `interval` spans this feature,
            partially or one is fully contained in the other.
        """
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

    cpdef bint spans_in_time(self, double time):
        """Test whether or not `time` is between :attr:`start_time`
        and :attr:`end_time`.

        Parameters
        ----------
        time : float
            The time to test for

        Returns
        -------
        bool:
            Whether or not the time is contained
        """
        return self.get_start_time() <= time <= self.get_end_time()

    def as_arrays(self, dtype=np.float64):
        """Convert this object into a pair of time and intensity arrays, summing
        all peaks at the same time point.

        Returns
        -------
        times : np.ndarray
            The time array for this feature.
        intensities : np.ndarray
            The intensity array for this feature.
        """
        rts = np.array(
            [node.time for node in self.nodes], dtype=dtype)
        signal = np.array([node.total_intensity()
                           for node in self.nodes], dtype=dtype)
        return rts, signal

    def __repr__(self):
        return "%s(%0.4f, %0.2f, %0.2f)" % (
            self.__class__.__name__, self.mz, self.start_time, self.end_time)

    def _copy_chunk(self, nodes, *args, **kwargs):
        x = self.__class__(LCMSFeatureTreeList(nodes))
        x.used_as_adduct = list(self.used_as_adduct)
        return x

    def split_at(self, time):
        """Split this feature at a specific time point into two pieces, before
        and after that point in time.

        Parameters
        ----------
        time : float
            The time point to split at

        Returns
        -------
        before : LCMSFeature
            The component up to the provided time
        after : LCMSFeature
            The component after the provided time
        """
        _, i = self.nodes.find_time(time)
        if self.nodes[i].time < time:
            i += 1
        return LCMSFeature(self.nodes[:i]), LCMSFeature(self.nodes[i:])

    def between_times(self, start, end):
        cdef:
            size_t i, j
            LCMSFeatureTreeNode node
        i = 0
        j = 0
        node = self._find_time(start, &i)
        if node is not None and node.time < start:
            i += 1
        node = self._find_time(end, &j)
        return self._copy_chunk(self.nodes[i:j + 1])

    cpdef list split_sparse(self, double delta_rt=1.):
        """Split this feature at any point where the distance from the
        time between one observed peak and the next exceeds `delta_rt`.

        Parameters
        ----------
        delta_rt : float
            The maximum time between observed peaks to tolerate before
            breaking the feature.

        Returns
        -------
        list of :class:`LCMSFeature`
            The feature components broken at the gaps.
        """
        cdef:
            list chunks, current_chunk
            double last_rt
            LCMSFeatureTreeNode node
            size_t i, n

        chunks = []
        current_chunk = []
        n = self.get_size()
        if n == 0:
            return chunks

        last_rt = self.getitem(0).time
        for i in range(n):
            node = self.getitem(i)
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

    cpdef LCMSFeature clone(self, deep=False, cls=None):
        if cls is None:
            cls = self.__class__
        c = cls(self.nodes.clone(deep=deep), list(self.adducts), list(self.used_as_adduct))
        c.feature_id = self.feature_id
        c.created_at = self.created_at
        return c

    cpdef insert_node(self, LCMSFeatureTreeNode node):
        """Insert a fully formed :class:`LCMSFeatureTreeNode` into this feature
        at the appropriate locaton in the time sequence.

        The node's :attr:`members` attribute will update the peak averager and
        triggers state invalidation.

        Parameters
        ----------
        node : LCMSFeatureTreeNode
            The node to add.
        """
        self._peak_averager._update(node.members)
        self.nodes.insert_node(node)
        self._invalidate()

    cpdef insert(self, PeakBase peak, double time):
        """Insert a peak at a known time into this feature.

        The peak updates the peak averager and triggers state invalidation.

        Parameters
        ----------
        peak : PeakBase
            The a peak with known m/z and intensity
        time : float
            The time point at which this peak was observed
        """
        self._peak_averager.add(peak)
        self.nodes.insert(time, [peak])
        self._invalidate()

    def merge(self, other):
        """Create a copy of this feature and copies all nodes from `other` into it,
        combining information from both features in a new instance.

        Parameters
        ----------
        other : FeatureBase
            Another feature whose nodes (which are distinct from this feature's)
            will be copied into a new copy of this feature.

        Returns
        -------
        LCMSFeature
        """
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

    cpdef double max_intensity(self):
        cdef:
            size_t i, n
            double max_intensity
        max_intensity = 0.0
        n = self.get_size()
        for i in range(n):
            max_intensity = max(self.getitem(i).max_intensity(), max_intensity)
        return max_intensity


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
        if self.index_list == NULL:
            raise MemoryError("Failed to allocate index list for LCMS Feature Iterator")
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

    cdef int _initialize(self, list features) except 1:
        cdef:
            size_t i, n
        self.features = features
        n = PyList_GET_SIZE(self.features)
        self.real_features = []
        self.index_list = <size_t*>malloc(sizeof(size_t) * n)
        if self.index_list == NULL:
            PyErr_SetString(MemoryError, "Failed to allocate index list for LCMS Feature Iterator")
            return 1
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
        return 0

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
            ix = 0
            if f.get_size() > 0:
                node = f._find_time(self.start_time, &ix)
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
                    ix = 0
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

    cpdef _initialize(self):
        self.current_mean = 0
        self.total_weight = 0
        self.current_count = 0

    def __init__(self, iterable=None):
        self._initialize()
        if iterable is not None:
            self.update(iterable)

    cpdef reset(self):
        self._initialize()

    @staticmethod
    cdef RunningWeightedAverage _create(list peaks):
        cdef:
            RunningWeightedAverage inst

        inst = RunningWeightedAverage.__new__(RunningWeightedAverage)
        inst._initialize()
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

        agg = (self.total_weight * self.current_mean) + \
            (peak.mz * peak.intensity)
        self.total_weight += peak.intensity
        self.current_count += 1

        if self.total_weight != 0:
            self.current_mean = agg / (self.total_weight)
        else:
            print("NaN produced in add()")
        return self

    cpdef int feed_from_feature(self, LCMSFeature feature):
        cdef:
            size_t i, j, n, m
            LCMSFeatureTreeNode node
            PeakBase peak

        n = feature.get_size()
        for i in range(n):
            node = feature.getitem(i)
            m = node.get_members_size()
            for j in range(m):
                peak = node.getitem(j)
                self.add(peak)
        return n

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
        return "%s(%r, %d)" % (self.__class__.__name__, self.current_mean, self.current_count)


cdef class RunningWeightedAverageNeutralMass(RunningWeightedAverage):

    @staticmethod
    cdef RunningWeightedAverageNeutralMass _create(list peaks):
        cdef:
            RunningWeightedAverageNeutralMass inst

        inst = RunningWeightedAverageNeutralMass.__new__(RunningWeightedAverageNeutralMass)
        inst._initialize()
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
                self.current_mean = (<DeconvolutedPeak?>peak).neutral_mass
                self.total_weight = 1
            else:
                return

        agg = (self.total_weight * self.current_mean) + \
            ((<DeconvolutedPeak?>peak).neutral_mass * peak.intensity)
        self.total_weight += peak.intensity
        self.current_count += 1

        if self.total_weight != 0:
            self.current_mean = agg / (self.total_weight)
        else:
            print("NaN produced in add()")
        return self
