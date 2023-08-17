import logging
from collections import defaultdict
from operator import attrgetter
from functools import total_ordering

from typing import Any, List, Dict, DefaultDict, Tuple, Optional, TYPE_CHECKING, Union, NamedTuple
from concurrent import futures

from ms_deisotope.data_source.scan import Scan, ProcessedScan, PrecursorInformation

# from .similarity_methods import SpectrumAlignment
from ms_deisotope._c.alignment import SpectrumAlignment, _AlignableSpectrum


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class Neighbor(NamedTuple):
    score: float
    shift: float
    scan_id: str
    weight: float


class Edge(NamedTuple):
    score: float
    shift: float
    weight: float


class Shift(NamedTuple):
    shift: float
    weight: float
    score: float


@total_ordering
class ShiftBin(object):
    mass: float
    observations: List[Shift]

    def __init__(self, mass, observations=None):
        if observations is None:
            observations = []
        self.mass = mass
        self.observations = observations

    def __eq__(self, other):
        return abs(self.mass - other.mass) < 1e-3

    def __lt__(self, other):
        return self.mass < other.mass

    def __len__(self):
        return len(self.observations)

    def __getitem__(self, i):
        return self.observations[i]

    def add(self, observation: Shift):
        self.observations.append(observation)

    def average(self) -> float:
        total = 0.0
        norm = 0.0
        for mass, weight, _score in self.observations:
            total += mass * weight
            norm += weight
        if norm == 0:
            return 0.0
        return total / norm

    def test(self, mass: float, max_delta: float = 0.01) -> bool:
        return abs(self.mass - mass) < max_delta

    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass}, {self.observations})"


GET_NEUTRAL_MASS = attrgetter("neutral_mass")
GET_EXTRACTED_NEUTRAL_MASS = attrgetter("extracted_neutral_mass")


class SpectrumAlignmentGraph(object):
    """
    A graph of spectra that can find similar pairs of shifted spectra from
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
    """

    scans: List[Union[Scan, ProcessedScan]]
    _scans: List[_AlignableSpectrum]
    threshold: float
    error_tolerance: float
    min_delta: float
    max_delta: float
    min_peaks_matched: int
    edges: DefaultDict[str, Dict[str, Edge]]
    supporters: Dict[str, List[Neighbor]]

    def __init__(self, scans, threshold=0.5, error_tolerance=2e-5, min_delta=0.01, max_delta=1e3, min_peaks_matched=6, match_charge=True,
                 edges=None, supporters=None):
        self.scans = list(scans)
        self._scans = []
        self.threshold = threshold
        self.error_tolerance = error_tolerance
        self.min_delta = min_delta
        self.max_delta = max_delta
        self.min_peaks_matched = min_peaks_matched
        self.match_charge = match_charge

        self.edges = defaultdict(dict, edges or {})
        self.supporters = dict(supporters or {})

        if not self.edges:
            self.build_edges_concurrent()
            self.trim()

    def _build_edges_for(self, i: int, edges: DefaultDict[Any, Dict[Any, Edge]]):
        n_edges = 0
        scan1 = self.scans[i]
        _scan1 = self._scans[i]
        pinfo1_mass = _scan1.precursor_mass
        for j in range(i + 1, len(self.scans)):
            scan2 = self.scans[j]
            if scan1 is scan2:
                continue
            _scan2 = self._scans[j]

            if self.match_charge and _scan1.precursor_charge != _scan2.precursor_charge:
                continue

            delta = _scan2.precursor_mass - pinfo1_mass
            if abs(delta) < self.min_delta:
                continue
            elif delta > self.max_delta:
                break

            if scan1.id in edges and scan2.id in edges[scan1.id]:
                continue

            aln = SpectrumAlignment(
                _scan1,
                _scan2,
                error_tolerance=self.error_tolerance,
                shift=delta,
            )
            if aln.score < self.threshold:
                continue
            if len(aln) < self.min_peaks_matched:
                continue
            weight = aln.score
            edges[scan1.id][scan2.id] = Edge(aln.score, aln.shift, weight)
            edges[scan2.id][scan1.id] = Edge(aln.score, -aln.shift, weight)
            n_edges += 1
        return i, n_edges

    def build_edges_concurrent(self):
        edges = defaultdict(dict)
        use_extracted = (self.scans and self.scans[0].precursor_information.extracted_neutral_mass != 0)
        sort_predicate = (lambda x: x.precursor_information.extracted_neutral_mass) if use_extracted \
                            else (lambda x: x.precursor_information.neutral_mass)
        self.scans.sort(key=sort_predicate)
        self._scans = list(map(_AlignableSpectrum, self.scans))
        promises: List[futures.Future] = []
        executor = futures.ThreadPoolExecutor(6)
        for i in range(len(self.scans)):
            promises.append(executor.submit(self._build_edges_for, i, edges))
        n_edges = 0
        for promise in promises:
            i, edges_added = promise.result()
            n_edges += edges_added
            if i % 1000 == 0 and i:
                logger.info("Processing scan %d, %d edges", i, n_edges)
        self.edges = edges

    def build_edges(self):
        edges = defaultdict(dict)
        use_extracted = (self.scans and self.scans[0].precursor_information.extracted_neutral_mass != 0)
        sort_predicate = (lambda x: x.precursor_information.extracted_neutral_mass) if use_extracted \
                            else (lambda x: x.precursor_information.neutral_mass)
        self.scans.sort(key=sort_predicate)
        self._scans = list(map(_AlignableSpectrum, self.scans))
        n_edges = 0
        for i, scan1 in enumerate(self.scans):
            if i % 1000 == 0 and i:
                logger.info("Processing scan %d, %d edges", i, n_edges)
            _scan1 = self._scans[i]
            pinfo1_mass = _scan1.precursor_mass
            for j in range(i + 1, len(self.scans)):
                scan2 = self.scans[j]
                if scan1 is scan2:
                    continue
                _scan2 = self._scans[j]

                if self.match_charge and _scan1.precursor_charge != _scan2.precursor_charge:
                    continue

                delta = _scan2.precursor_mass - pinfo1_mass
                if abs(delta) < self.min_delta:
                    continue
                elif delta > self.max_delta:
                    break

                if scan1.id in edges and scan2.id in edges[scan1.id]:
                    continue

                aln = SpectrumAlignment(
                    _scan1,
                    _scan2,
                    error_tolerance=self.error_tolerance,
                    shift=delta,
                )
                if aln.score < self.threshold:
                    continue
                if len(aln) < self.min_peaks_matched:
                    continue
                weight = aln.score
                edges[scan1.id][scan2.id] = Edge(aln.score, aln.shift, weight)
                edges[scan2.id][scan1.id] = Edge(aln.score, -aln.shift, weight)
                n_edges += 1
        self.edges = edges

    def iteredges(self):
        for s1, s2_edges in self.edges.items():
            for s2, edge in s2_edges.items():
                yield s1, s2, edge

    def trim(self):
        supporters = {}
        last_scan1 = None
        buckets = None
        for scan1_id, scan2_id, edge in self.iteredges():
            if scan1_id != last_scan1:
                if buckets is not None:
                    supporters[last_scan1] = buckets
                buckets = {}
                last_scan1 = scan1_id
                if edge.score < self.threshold:
                    continue
                key = int(edge.shift * 100.0)
                if key in buckets and buckets[key][0] < edge.score:
                    buckets[key] = Neighbor(edge.score, edge.shift, scan2_id, edge.weight)
                else:
                    buckets[key] = Neighbor(
                        edge.score, edge.shift, scan2_id, edge.weight)
        if buckets is not None:
            supporters[last_scan1] = buckets
        self.supporters = supporters

    def collect_shifts(self) -> List[Shift]:
        seen = set()

        shifts: List[Shift] = []
        for a, b, ed in self.iteredges():
            key = frozenset((a, b))
            if key in seen:
                continue
            seen.add(key)
            (score, shift, weight) = ed
            if score < self.threshold:
                continue
            shift = abs(shift)
            shifts.append(Shift(shift, weight, score))

        shifts.sort(key=lambda x: x[1], reverse=True)
        return shifts

    def average_shifts(self, error_tolerance: float=0.1) -> List[ShiftBin]:
        shifts = self.collect_shifts()

        bins = []
        for shift in shifts:
            for sbin in bins:
                if sbin.test(shift[0], error_tolerance):
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
