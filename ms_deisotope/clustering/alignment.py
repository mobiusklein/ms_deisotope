from collections import defaultdict, namedtuple

from .similarity_methods import SpectrumAlignment

Neighbor = namedtuple("Neighbor", ('score', 'shift', 'scan_id', 'weight'))


class SpectrumAlignmentGraph(object):
    '''A graph of spectra that can find similar pairs of shifted spectra from
    the population, constructing a support network for each spectrum.

    Attributes
    ----------
    scans : list of :class:`~.ScanBase`
        The spectra to be related.
    threshold : float
        The minimum alignment similarity to accept for supporters
    error_tolerance : float
        The ppm error tolerance to use when matching peaks, defaults to 2e-5
    min_delta : float
        The minimum precursor mass delta recognized. Below this number and the
        precursor mass delta is assumed to be the same and is skipped.
    match_charge : bool
        Whether or not to require charge states to match to construct an edge.
    '''

    def __init__(self, scans, threshold=0.5, error_tolerance=2e-5, min_delta=0.01, match_charge=True):
        self.scans = list(scans)
        self.threshold = threshold
        self.error_tolerance = error_tolerance
        self.min_delta = min_delta
        self.match_charge = match_charge

        self.edges = defaultdict(dict)
        self.supporters = dict()

        self.build_edges()
        self.trim()

    def build_edges(self):
        edges = defaultdict(dict)
        self.scans.sort(key=lambda x: x.precursor_information.neutral_mass)
        for i, scan1 in enumerate(self.scans):
            for scan2 in self.scans[i + 1:]:
                if scan1 is scan2:
                    continue
                if self.match_charge and scan1.precursor_information.charge != scan2.precursor_information.charge:
                    continue
                delta = scan2.precursor_information.neutral_mass - \
                    scan1.precursor_information.neutral_mass
                if abs(delta) < self.min_delta:
                    continue
                aln = SpectrumAlignment(
                    scan1.peaks(), scan2.peaks(),
                    error_tolerance=self.error_tolerance,
                    shift=scan2.precursor_information.neutral_mass - scan1.precursor_information.neutral_mass)
                weight = sum([pp.score for pp in aln.peak_pairs])
                edges[scan1.id][scan2.id] = (aln.score, aln.shift, weight)
                edges[scan2.id][scan1.id] = (aln.score, -aln.shift, weight)
        self.edges = edges

    def trim(self):
        supporters = {}
        for scan1_id, neighbors in self.edges.items():
            buckets = {}
            for scan2_id, (score, shift, weight) in neighbors.items():
                if score < self.threshold:
                    continue
                key = int(shift * 100.0)
                if key in buckets and buckets[key][0] < score:
                    buckets[key] = Neighbor(score, shift, scan2_id, weight)
                else:
                    buckets[key] = Neighbor(score, shift, scan2_id, weight)
            supporters[scan1_id] = list(buckets.values())
        self.supporters = supporters

    def average_shifts(self):
        seen = set()

        shifts = []

        for a, dst in self.edges.items():
            for b, ed in dst.items():
                key = frozenset((a, b))
                if key in seen:
                    continue
                seen.add(key)
                (score, shift, weight) = ed
                if score < self.threshold:
                    continue
                shift = abs(shift)
                shifts.append((shift, weight, score))

        shifts.sort(key=lambda x: x[1], reverse=True)

        bins = []
        for shift in shifts:
            for sbin in bins:
                if sbin.test(shift[0], 0.1):
                    sbin.add(shift)
                    break
            else:
                bins.append(ShiftBin(shift[0], [shift]))

        bins.sort(key=lambda x: x.average())
        return bins

    def to_json(self):
        container = {
            "edges": dict(self.edges),
            "supporters": self.supporters,
            "scan_ids": [s.id for s in self.scans],
            "parameters": {
                "threshold": self.threshold,
                "error_tolerance": self.error_tolerance,
                "min_delta": self.min_delta,
                "match_charge": self.match_charge
            }
        }
        return container


class ShiftBin(object):
    def __init__(self, mass, observations=None):
        if observations is None:
            observations = []
        self.mass = mass
        self.observations = observations

    def add(self, observation):
        self.observations.append(observation)

    def average(self):
        total = 0.0
        norm = 0.0
        for mass, weight in self.observations:
            total += mass * weight
            norm += weight
        if norm == 0:
            return 0.0
        return total / norm

    def test(self, mass, max_delta=0.01):
        return abs(self.mass - mass) < max_delta

    def __repr__(self):
        return "{self.__class__.__name__}({self.mass}, {self.observations})"
