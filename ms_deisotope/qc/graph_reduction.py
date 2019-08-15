from ms_deisotope.peak_set import DeconvolutedPeakSet
from ms_deisotope.spectrum_graph import (PathFinder, amino_acids, MassWrapper)


class GraphBasedDenoiser(object):
    def __init__(self, components=None, precursor_error_tolerance=2e-5, product_error_tolerance=2e-5):
        if components is None:
            components = amino_acids
        self.precursor_error_tolerance = precursor_error_tolerance
        self.product_error_tolerance = product_error_tolerance
        self.path_finder = PathFinder(components, product_error_tolerance)
        self.components = self.path_finder.components

    def initialize_graph(self, peak_set, precursor_mass):
        graph = self.path_finder.find_edges(peak_set)
        p = MassWrapper("Precursor", precursor_mass)
        self.path_finder.find_complements(peak_set, graph, p)
        return graph

    def walk_graph(self, graph):
        keep = set()
        for e in graph.transitions:
            keep.add(e.start.index)
            keep.add(e.end.index)
        return keep

    def select_peaks(self, peak_set, indices):
        keep = []
        for i in indices:
            keep.append(peak_set[i].clone())
        keep = DeconvolutedPeakSet(keep)
        keep.reindex()
        return keep

    def filter(self, scan):
        precursor_mass = scan.precursor_information.extracted_neutral_mass
        if precursor_mass == 0:
            precursor_mass = scan.precursor_information.neutral_mass
        peak_set = scan.deconvoluted_peak_set
        graph = self.initialize_graph(peak_set, precursor_mass)
        indices = self.walk_graph(graph)
        peak_set = self.select_peaks(peak_set, indices)
        return peak_set

    def __call__(self, scan):
        return self.filter(scan)
