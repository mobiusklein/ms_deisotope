from .utils import GeneratorQueue
from .intervals import SpanningMixin


def ident(x):
    return x


class FitNode(SpanningMixin):
    def __init__(self, fit, index=None):
        self.index = index
        self.fit = fit
        self.edges = set()
        self.overlap_edges = set()
        # Placeholder Peaks have a peak_count of -1
        self.peak_indices = {p.peak_count for p in fit.experimental if p.peak_count >= 0}
        self._hash = None
        self.score = fit.score
        self.start = fit.experimental[0].mz
        self.end = fit.experimental[-1].mz

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.fit)
        return self._hash

    def __eq__(self, other):
        return self.fit is other.fit

    def __gt__(self, other):
        return self.fit > other.fit

    def __lt__(self, other):
        return self.fit < other.fit

    def __ne__(self, other):
        return self.fit is not other.fit

    def isdisjoint(self, other):
        return self.peak_indices.isdisjoint(other.peak_indices)

    def visit(self, other):
        if self.isdisjoint(other):
            self.edges.add(other)
            other.edges.add(self)
        else:
            self.overlap_edges.add(other)
            other.overlap_edges.add(self)

    def __repr__(self):
        return "FitNode(%r)" % self.fit


def tupleslices(iterable, shard_size=4):
    i = 0
    t = len(iterable)
    while i < t:
        yield iterable[i:i + shard_size]
        i += shard_size


def span_overlap(a, b):
    return a.__contains__(b.start) or a.__contains__(b.end - 1) or\
        b.__contains__(a.start) or b.__contains__(a.end - 1)


def peak_overlap(a, b):
    return len(a.peak_indices & b.peak_indices) > 0


def layout_layers(envelopes, overlap_fn=peak_overlap, maximize=True):
    '''
    Produce a non-overlapping stacked layout of individual envelopes.
    '''
    layers = [[]]
    envelopes.sort(key=lambda x: x.score, reverse=maximize)
    for envelope in envelopes:
        placed = False
        for layer in layers:
            collision = False
            for member in layer:
                if overlap_fn(envelope, member):
                    collision = True
                    break
            if not collision:
                layer.append(envelope)
                placed = True
                break
        if not placed:
            layers.append([envelope])
    return layers


class GreedySubgraphSelection(object):
    def __init__(self, subgraph, maximize=True, overlap_fn=peak_overlap):
        self.intervals = list(subgraph)
        self._layers = [[]]
        self.maximize = maximize
        self._build_layers(overlap_fn=overlap_fn)

    def _build_layers(self, overlap_fn):
        self._layers = layout_layers(self.intervals, overlap_fn=overlap_fn, maximize=self.maximize)

    def _select_best_subset(self):
        return self._layers[0]

    def select(self):
        return self._select_best_subset()

    @classmethod
    def solve(cls, nodes, maximize=True, overlap_fn=peak_overlap):
        solver = cls(nodes, maximize=maximize, overlap_fn=overlap_fn)
        solution = solver.select()
        return solution


class ExhaustiveDisjointSolutionSelection(object):  # pragma: no cover
    shard_size = 7

    def __init__(self, active_set, remaining_components, solutions=None, parent=None, score=None):
        if solutions is None:
            solutions = dict()
        if score is None:
            score = sum(f.score for f in active_set)
        self.active_set = active_set
        self.remaining_components = remaining_components
        self.solutions = solutions
        self.accepted = []
        self.rejected = []
        self.possible_next = []
        self.score = score

    def is_valid_addition(self, current, addition):
        is_valid = True
        for member in current:
            if not member.isdisjoint(addition):
                is_valid = False
                break
        return is_valid

    def check_has_solution(self, active_set, addition):
        for chunk in tupleslices(active_set, self.shard_size):
            try:
                check = self.solutions[chunk, addition.index]
            except KeyError:
                check = self.solutions[chunk, addition.index] = self.is_valid_addition(chunk, addition)
            if not check:
                return False
        return True

    def find_acceptable_peaks(self):
        accepted = []
        rejected = []
        for comp in self.remaining_components:
            if self.check_has_solution(self.active_set, comp):
                accepted.append(comp)
            else:
                rejected.append(comp)
        self.accepted = accepted
        self.rejected = rejected
        self.possible_next = self.remaining_components - set(rejected)

    def solve_level(self):
        if len(self.accepted) == 0:
            self.find_acceptable_peaks()
        accepted = self.accepted
        # next_steps = []
        current_step = set(self.active_set)
        extensions = self.possible_next
        for acc_ in accepted:
            acc = {acc_}
            next_set = tuple(sorted(current_step | acc))
            next_components = (extensions - acc) - acc_.overlaps
            new_score = self.score + acc_.score
            # total_possible_score = new_score + sum(f.score for f in next_components)
            # if total_possible_score > self.get_best_score() - 1:
            x = self.__class__(next_set, next_components, solutions=self.solutions,
                               parent=self, score=new_score)
            yield x

    def total_possible_score(self):
        if len(self.accepted) == 0:
            self.find_acceptable_peaks()
        if self.accepted:
            accepted_score = max(f.score for f in self.accepted)
        else:
            accepted_score = 0
        return self.score + sum(f.score for f in self.possible_next) + accepted_score

    @classmethod
    def solve(cls, nodes):
        initial_set = tuple()
        best_score = 0
        best_solution = None
        work_queue = GeneratorQueue()
        work_queue.extend(iter([cls(initial_set, set(nodes))]))
        i = 0
        # t = time.time()
        while work_queue.hasnext():
            i += 1
            solution = work_queue.next()
            solution.find_acceptable_peaks()

            if i % 300000 == 0 or solution.score > best_score and best_score != 0:
                # print i, best_score, len(best_solution.active_set), len(
                #     best_solution.remaining_components), len(
                #     best_solution.solutions), len(work_queue), best_solution.total_possible_score()
                # print '----', time.time() - t
                pass
            if solution.total_possible_score() < best_score:
                continue

            if solution.score > best_score:
                best_score = solution.score
                best_solution = solution

            work_queue.extend(solution.solve_level())
        return best_solution


class ConnectedSubgraph(object):
    def __init__(self, fits, maximize=True):
        self.nodes = tuple(map(FitNode, fits))
        self.maximize = maximize
        i = 0
        for node in self.nodes:
            node.index = i
            i += 1
        self.populate_edges()

    def __getitem__(self, i):
        return self.nodes[i]

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def populate_edges(self):
        for i in range(len(self.nodes)):
            node = self.nodes[i]
            for j in range(len(self.nodes)):
                if i == j:
                    continue
                other = self.nodes[j]
                node.visit(other)

    def find_heaviest_path(self, method="greedy"):
        if len(self) == 1:
            return set(self.nodes)
        if method == "greedy":
            solution = GreedySubgraphSelection.solve(tuple(self), maximize=self.maximize)
            return solution
        else:
            raise NotImplementedError(method)


try:
    has_c = True
    _ConnectedSubgraph = ConnectedSubgraph
    _FitNode = FitNode

    from ms_deisotope._c.peak_dependency_network.subgraph import ConnectedSubgraph, FitNode
except ImportError:
    has_c = False
