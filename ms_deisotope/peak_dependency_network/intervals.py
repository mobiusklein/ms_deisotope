from collections import deque


class SpanningMixin(object):
    """Provides methods for checking whether an entity
    which has a defined start and end point over a single
    dimension contains or overlaps with another entity in
    that same dimension.
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __contains__(self, i):
        """Tests for point inclusion, `start <= i <= end`

        Parameters
        ----------
        i : Number
            The point to be tested

        Returns
        -------
        bool
        """
        return self.contains(i)

    def contains(self, i):
        """Tests for point inclusion, `start <= i <= end`

        Parameters
        ----------
        i : Number
            The point to be tested

        Returns
        -------
        bool
        """
        return self.start <= i <= self.end

    def overlaps(self, interval):
        """Tests whether another spanning entity with
        a defined start and end point overlaps with
        this spanning entity

        Parameters
        ----------
        interval : SpanningMixin

        Returns
        -------
        bool
        """
        cond = ((self.start <= interval.start and self.end >= interval.start) or (
            self.start >= interval.start and self.end <= interval.end) or (
            self.start >= interval.start and self.end >= interval.end and self.start <= interval.end) or (
            self.start <= interval.start and self.end >= interval.start) or (
            self.start <= interval.end and self.end >= interval.end))
        return cond

    def overlap_size(self, interval):
        if interval.start in self and interval.end in self:
            return interval.end - interval.start
        elif self.start in interval and self.end in interval:
            return self.end - self.start
        if self.end >= interval.start:
            return self.end - interval.start
        elif self.start >= interval.end:
            return self.start - interval.end

    def is_contained_in_interval(self, interval):
        """Tests whether this spanning entity is
        completely contained inside the other entity

        Parameters
        ----------
        interval : SpanningMixin

        Returns
        -------
        bool
        """
        return self.start >= interval.start and self.end <= interval.end

    def contains_interval(self, interval):
        """Tests whether the other spanning entity is
        completely contained inside this entity

        Parameters
        ----------
        interval : SpanningMixin

        Returns
        -------
        bool
        """
        return self.start <= interval.start and self.end >= interval.end


class Interval(SpanningMixin):
    """A generic wrapper around data associated
    with a particular interval over a single dimension

    Attributes
    ----------
    data : dict
        A holder of arbitrary extra data not associated
        with any single member contained in this interval
    end : Number
        The end point of the interval described
    members : list
        A list of arbitrary objects which are associated with this
        interval.
    start : Number
        The start point of the interval described
    """
    def __init__(self, start, end, members=None, **kwargs):
        if members is None:
            members = []
        self.start = start
        self.end = end
        self.members = members
        self.data = kwargs

    def __repr__(self):
        return "%s(start=%r, end=%r, data=%r)" % (self.__class__.__name__, self.start, self.end, self.data)

    def __getitem__(self, k):
        return self.members[k]

    def __len__(self):
        return len(self.members)

    def __iter__(self):
        return iter(self.members)


class IntervalTreeNode(object):
    """A single node in an Interval Tree. The
    root node of an interval tree can be treated
    as the tree itself.

    Attributes
    ----------
    center : Number
        The center point of this node's collection
    contained : list
        The list of Interval-like objects which
        were aggregated under this node. This collection
        of intervals all span this node's center.
    end : Number
        The end point of this node's collection or it's
        rightmost child's end point, whichever is larger
    left : IntervalTreeNode
        The left child node of this node. May be `None`
        if this node is a leaf node. Contains all
        intervals whose end point is less than this
        node's center
    level : int
        Depth in the tree
    parent : IntervalTree
        This node's parent node. May be None if this
        node is the root node.
    right : IntervalTree
        The right child node of this node. May be `None`
        if this node is a leaf node. Contains all
        intervals whose start point is greater than this
        node's center
    start : Number
        The start point of this node's collection or it's
        leftmost child's start point, whichever is smaller
    """
    def __init__(self, center, left, contained, right, level=0, parent=None):
        self.center = center
        self.left = left
        self.contained = contained
        self.right = right
        self.level = level
        self.parent = parent

        start = float('inf')
        end = -float('inf')
        i = 0

        for interval in self.contained:
            i += 1
            if interval.start < start:
                start = interval.start
            if interval.end > end:
                end = interval.end

        if i > 0:
            self.start = start
            self.end = end
        else:
            self.start = self.end = center

    def contains_point(self, x):
        """Returns the list of contained intervals for all
        nodes which contain the point `x`.

        Parameters
        ----------
        x : Number
            The query point

        Returns
        -------
        list
            A list of objects which span `x` from all spanning
            nodes' `contained` list.
        """
        inner = [
            i for i in self.contained
            if i.start <= x <= i.end
        ]
        if self.left is not None and self.left.start <= x <= self.left.end:
            return self.left.contains_point(x) + inner
        if self.right is not None and self.right.start <= x <= self.right.end:
            return self.right.contains_point(x) + inner
        else:
            return inner

    def _overlaps_interval(self, start, end):
        result = []
        query = Interval(start, end)
        for i in self.contained:
            if query.overlaps(i):
                result.append(i)
        return result

    def overlaps(self, start, end):
        """Returns the list of all contained intervals which
        overlap with the interval described by `start` and `end`

        Parameters
        ----------
        start : Number
            The start of the query interval
        end : Number
            The end of the query interval

        Returns
        -------
        list
            A list of all objects which overlap the argument interval
            from all spanning nodes' `contained` list.
        """
        result = []
        if start > self.end:
            return result
        if end < self.start:
            return result
        elif start <= self.start:
            if end < self.start:
                return result
            else:
                if self.left is not None:
                    result.extend(self.left.overlaps(start, end))
                result.extend(self._overlaps_interval(start, end))
                if self.right is not None and end >= self.right.start:
                    result.extend(self.right.overlaps(start, end))
        elif start > self.start:
            if self.left is not None and self.left.end >= start:
                result.extend(self.left.overlaps(start, end))
            result.extend(self._overlaps_interval(start, end))
            if self.right is not None and end >= self.right.start:
                result.extend(self.right.overlaps(start, end))
        elif end > self.start:
            if self.left is not None:
                result.extend(self.left.overlaps(start, end))
            result.extend(self._overlaps_interval(start, end))
            if self.right is not None and end >= self.right.start:
                result.extend(self.right.overlaps(start, end))
        return result

    def node_contains(self, point):
        return self.start <= point <= self.end

    def __repr__(self):
        return "IntervalTreeNode(level=%d, center=%0.4f, start=%0.4f, end=%0.4f)" % (
            self.level, self.center, self.start, self.end)

    def __eq__(self, other):
        if other is None:
            return False
        else:
            result = self.contained == other.contained
            if result:
                result = self.left == other.left
                if result:
                    result = self.right == other.right
            return result

    def __hash__(self):
        return hash((self.start, self.center, self.right, self.level))

    def __diagnostic_eq__(self, other):  # pragma: no cover
        if other is None:
            return False
        else:
            result = self.contained == other.contained
            if result:
                try:
                    result = self.left.__diagnostic_eq__(other.left)
                except AttributeError:
                    result = self.left == other.left
                if result:
                    try:
                        result = self.right.__diagnostic_eq__(other.right)
                    except AttributeError:
                        result = self.right == other.right
                    if not result:
                        print(self.right, "r!=", other.right)
                else:
                    print(self.left, "l!=", other.left)
            else:
                print(self, "!=", other)
            return result

    def _update_bounds(self):
        changed = False
        if self.left is not None and self.left.start < self.start:
            changed = True
            self.start = self.left.start
        if self.right is not None and self.right.end > self.end:
            changed = True
            self.end = self.right.end
        if self.parent is not None:
            self.parent._update_bounds()
        return changed

    def insert(self, interval):
        insert_in_self = False
        if self.node_contains(interval.start):
            if self.left is not None and self.left.node_contains(interval.end):
                return self.left.insert(interval)
            elif (self.right is not None and self.right.node_contains(interval.end) and
                  self.right.node_contains(interval.start)):
                return self.right.insert(interval)
            else:
                insert_in_self = True
        elif self.node_contains(interval.end):
            if self.right is not None and self.right.node_contains(interval.start):
                return self.right.insert(interval)
            else:
                insert_in_self = True
        if not insert_in_self and self.parent is None:
            insert_in_self = True
        if insert_in_self:
            self.contained.append(interval)
            changed = False
            if interval.start < self.start:
                self.start = interval.start
                changed = True
            if interval.end > self.end:
                self.end = interval.end
                changed = True
            if changed:
                self._update_bounds()

    def flatten(self):
        # perform an infix traversal of the tree and collect
        # the elements
        items = []
        stack = deque([])
        current = self
        done = False
        while not done:
            if current is not None:
                stack.append(current)
                current = current.left
            else:
                if stack:
                    current = stack.pop()
                    items.extend(current.contained)
                    current = current.right
                else:
                    done = True
        return items

    def balance(self):
        items = self.flatten()
        tree = self.build(items)
        self.left = tree.left
        self.right = tree.right
        self.contained = tree.contained
        self.start = tree.start
        self.end = tree.end
        self.center = tree.center
        self.left.parent = self
        self.right.parent = self
        self._update_bounds()


def iterative_build_interval_tree(cls, intervals):
    """Builds an IntervalTreeNode hierarchy from `intervals`. This
    iterative method is preferable to avoid recursion limits.

    Parameters
    ----------
    intervals : Iterable of Intervals
        The set of Interval-like objects to use to construct
        the Interval Tree

    Returns
    -------
    IntervalTreeNode
        The root of the constructed Interval Tree
    """
    stack = []
    root = cls(0, None, [], None, -1)
    if not intervals:
        return root

    stack.append((root, intervals, "left"))

    while stack:
        parent, members, side = stack.pop()
        if len(members) > 0:
            centers = [
                float(i.start + i.end) / 2. for i in members
            ]
            left = []
            right = []
            contained = []
            center = sum(centers) / (len(centers))
            if len(members) > 20:
                for i in members:
                    if abs(i.start - center) < 1e-6 and abs(i.end - center) < 1e-6:
                        contained.append(i)
                    elif center > i.end:
                        left.append(i)
                    elif center < i.start:
                        right.append(i)
                    else:
                        contained.append(i)
            else:
                contained = members[:]

            if len(right) == len(members) or len(left) == len(members):
                contained = members[:]
                left = []
                right = []
            node = cls(center, left=None, contained=contained, right=None, level=parent.level + 1, parent=parent)
            if side == 'left':
                parent.left = node
                parent.start = min(node.start, parent.start)
                up = parent.parent
                while up is not None:
                    if up.start > node.start:
                        up.start = node.start
                    up = up.parent
            elif side == 'right':
                parent.right = node
                parent.end = max(node.end, parent.end)
                up = parent.parent
                while up is not None:
                    if up.end < node.end:
                        up.end = node.end
                    up = up.parent
            else:
                raise ValueError(side)
            stack.append((node, left, "left"))
            stack.append((node, right, "right"))
    return root.left


IntervalTreeNode.build = classmethod(iterative_build_interval_tree)


def recursive_build_interval_tree(cls, intervals, level=0):  # pragma: no cover
    if len(intervals) > 0:
        centers = [
            float(i.start + i.end) / 2. for i in intervals
        ]
        left = []
        right = []
        contained = []
        center = sum(centers) / (len(centers) + 1.)
        if len(intervals) > 20:
            for i in intervals:
                if center >= i.end:
                    left.append(i)
                elif center < i.start:
                    right.append(i)
                else:
                    contained.append(i)
            left = cls.recursive_build_interval_tree(left, level + 1)
            right = cls.recursive_build_interval_tree(right, level + 1)
        else:
            contained = intervals
            left = None
            right = None

        inst = cls(center, left, contained, right, level)
        if left is not None:
            inst.left.parent = inst
            inst.start = min(inst.start, left.start)
        if right is not None:
            inst.right.parent = inst
            inst.end = (max(inst.end, right.end))
        return inst
    else:
        return None


IntervalTreeNode.recursive_build_interval_tree = classmethod(
    recursive_build_interval_tree)

try:
    has_c = True
    _SpanningMixin = SpanningMixin
    _Interval = Interval
    _IntervalTreeNode = IntervalTreeNode
    from ms_deisotope._c.peak_dependency_network.intervals import (
        SpanningMixin, Interval, IntervalTreeNode)
except ImportError:
    has_c = False
