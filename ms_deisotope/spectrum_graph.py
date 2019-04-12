from ms_deisotope._c.spectrum_graph import (
    PathFinder,
    MassWrapper,
    PeakGroupNode,
    PeakNode,
    NodeBase,
    Path,
    SpectrumGraph)

amino_acids = [
    MassWrapper('G', 57.02146372057),
    MassWrapper('A', 71.03711378471),
    MassWrapper('S', 87.03202840427),
    MassWrapper('P', 97.05276384884999),
    MassWrapper('V', 99.06841391299),
    MassWrapper('T', 101.04767846841),
    MassWrapper('C', 103.00918478471),
    MassWrapper('J', 113.08406397713),
    MassWrapper('N', 114.04292744114),
    MassWrapper('D', 115.02694302383),
    MassWrapper('Q', 128.05857750528),
    MassWrapper('K', 128.094963014),
    MassWrapper('E', 129.04259308797),
    MassWrapper('M', 131.04048491299),
    MassWrapper('H', 137.05891185845),
    MassWrapper('F', 147.06841391299),
    MassWrapper('R', 156.1011110236),
    MassWrapper('Y', 163.06332853255),
    MassWrapper('W', 186.07931294986),
]

def find_paths(peaks, components=None, error_tolerance=1e-5, merge=False):
    if components is None:
        components = amino_acids
    sequencer = PathFinder(components, error_tolerance)
    paths = sequencer.paths(peaks, merge=merge)
    return paths
