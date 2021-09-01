from collections import defaultdict, namedtuple

from .similarity_methods import SpectrumAlignment

Neighbor = namedtuple("Neighbor", ('score', 'shift', 'scan_id'))


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
        for scan1 in self.scans:
            for scan2 in self.scans:
                if scan1 is scan2:
                    continue
                if self.match_charge and scan1.precursor_information.charge != scan2.precursor_information.charge:
                    continue
                if scan1.id in edges and scan2.id in edges[scan1.id]:
                    continue
                delta = scan2.precursor_information.neutral_mass - \
                    scan1.precursor_information.neutral_mass
                if abs(delta) < self.min_delta:
                    continue
                aln = SpectrumAlignment(
                    scan1.deconvoluted_peak_set, scan2.deconvoluted_peak_set, error_tolerance=self.error_tolerance,
                    shift=scan2.precursor_information.neutral_mass - scan1.precursor_information.neutral_mass)
                edges[scan1.id][scan2.id] = (aln.score, aln.shift)
                edges[scan2.id][scan1.id] = (aln.score, -aln.shift)
        self.edges = edges

    def trim(self):
        supporters = {}
        for scan1_id, neighbors in self.edges.items():
            buckets = {}
            for scan2_id, (score, shift) in neighbors.items():
                if score < self.threshold:
                    continue
                key = int(shift * 100.0)
                if key in buckets and buckets[key][0] < score:
                    buckets[key] = Neighbor(score, shift, scan2_id)
                else:
                    buckets[key] = Neighbor(score, shift, scan2_id)
            supporters[scan1_id] = list(buckets.values())
        self.supporters = supporters
