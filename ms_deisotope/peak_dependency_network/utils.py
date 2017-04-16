from collections import deque


def bisect_at_score(scored_dependencies, splitting_value):
    lo = 0
    hi = len(scored_dependencies)

    while hi != lo:
        mid = (hi + lo) / 2
        mid_fit = scored_dependencies[mid]
        if hi == lo + 1:
            return mid
        elif abs(mid_fit.score - splitting_value) < 1e-3:
            return mid
        elif mid_fit.score > splitting_value:
            lo = mid
        elif mid_fit.score < splitting_value:
            hi = mid


def bisect_at_start(clusters, mz):
    lo = 0
    hi = len(clusters)
    while hi != lo:
        mid = (hi + lo) / 2
        mid_cluster = clusters[mid]
        if hi == lo + 1:
            return mid
        elif abs(mid_cluster.start - mz) < 1e-3:
            return mid
        elif mid_cluster.start < mz:
            hi = mid
        elif mid_cluster.start > mz:
            lo = mid


class GeneratorQueue(object):  # pragma: no cover
    sentinel = object()

    def __init__(self, queue=None):
        if queue is None:
            queue = deque()
        self.queue = queue
        self.active = None
        self._peek = self.sentinel

    def extend(self, iterable):
        self.queue.appendleft(iter(iterable))

    def load_next_generator(self):
        self.active = self.queue.popleft()

    def get_next_value(self):
        if self._peek is not self.sentinel:
            value = self._peek
            self._peek = self.sentinel
        else:
            try:
                value = next(self.active)
            except (StopIteration, TypeError):
                try:
                    self.load_next_generator()
                    value = self.get_next_value()
                except IndexError:
                    print("Terminating with size", len(self.queue))
                    raise StopIteration()
        return value

    def __next__(self):
        value = self.get_next_value()
        return value

    def next(self):
        return self.__next__()

    def __iter__(self):
        while True:
            yield self.next()
        print("Terminating with size", len(self.queue))

    def __len__(self):
        return len(self.queue)

    @property
    def peek(self):
        if self._peek is self.sentinel:
            try:
                self._peek = self.get_next_value()
            except StopIteration:
                self._peek = self.sentinel
        return self._peek

    def hasnext(self):
        return self.peek is not self.sentinel
