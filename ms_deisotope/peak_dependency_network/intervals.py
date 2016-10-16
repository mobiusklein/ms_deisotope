

class SpanningMixin(object):
    def __contains__(self, i):
        return self.start <= i <= self.end

    def overlaps(self, interval):
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
        return self.start >= interval.start and self.end <= interval.end

    def contains_interval(self, interval):
        return self.start <= interval.start and self.end >= interval.end


class Interval(SpanningMixin):
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

    def node_contains_point(self, x):
        if x < self.start:
            if self.left is not None:
                return self.left.node_contains_point(x)
            return None
        elif x > self.end:
            if self.right is not None:
                return self.right.node_contains_point(x)
            return None
        else:
            return self

    def contains_point(self, x):
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
        for i in self.contained:
            if ((start < i.start and i.end <= end) or
                    (i.start <= start and i.end >= end)):
                result.append(i)
        return result

    def overlaps(self, start, end):
        result = []
        if (start < self.start and end >= self.end):
            if self.left is not None:
                result.extend(self.left.overlaps(start, end))
            result.extend(self._overlaps_interval(start, end))
            return result
        elif start >= self.start and end >= self.end:
            if self.right is not None:
                result.extend(self.right.overlaps(start, end))
            result.extend(self._overlaps_interval(start, end))
            return result
        elif start > self.start and end <= self.end:
            return self._overlaps_interval(start, end)
        elif self.start > end:
            if self.left is not None:
                return self.left.overlaps(start, end)
            return result
        elif self.end < start:
            if self.right is not None:
                return self.right.overlaps(start, end)
            return result
        else:
            if self.left is not None:
                result.extend(self.left.overlaps(start, end))
            result.extend(self._overlaps_interval(start, end))
            if self.right is not None:
                result.extend(self.right.overlaps(start, end))
            return result

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

    def __diagnostic_eq__(self, other):
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
                        print self.right, "r!=", other.right
                else:
                    print self.left, "l!=", other.left
            else:
                print self, "!=", other
            return result


def iterative_build_interval_tree(cls, intervals):
    stack = []
    root = cls(0, None, [], None, -1)

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


def recursive_build_interval_tree(cls, intervals, level=0):
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
