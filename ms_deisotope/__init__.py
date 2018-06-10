try:
    import dill
except ImportError:
    pass

from .version import version

from .averagine import (
    Averagine, peptide, glycan, glycopeptide, heparin,
    heparan_sulfate, permethylated_glycan,
    mass_charge_ratio, neutral_mass,
    calculate_mass, isotopic_shift)
from .deconvolution import (
    AveragineDeconvoluter, CompositionListDeconvoluter,
    AveraginePeakDependenceGraphDeconvoluter,
    CompositionListPeakDependenceGraphDeconvoluter,
    deconvolute_peaks)
from .scoring import MSDeconVFitter, PenalizedMSDeconVFitter, DistinctPatternFitter, IsotopicFitRecord
from .peak_set import DeconvolutedPeak, DeconvolutedPeakSet, DeconvolutedPeakSolution
from .processor import ScanProcessor
from .data_source import MzMLLoader, MzXMLLoader, MSFileLoader


def get_include():
    import os
    return os.path.join(
        os.path.dirname(__file__),
        "_c")


__all__ = [
    "Averagine", 'peptide', 'glycan', 'glycopeptide', 'heparin', "heparan_sulfate",
    "permethylated_glycan",
    "mass_charge_ratio", "neutral_mass", "isotopic_shift", "calculate_mass",
    "AveragineDeconvoluter", "CompositionListDeconvoluter",
    "AveraginePeakDependenceGraphDeconvoluter", "CompositionListPeakDependenceGraphDeconvoluter",
    "MSDeconVFitter", "PenalizedMSDeconVFitter", "DistinctPatternFitter", "IsotopicFitRecord",
    "DeconvolutedPeak", "DeconvolutedPeakSet", "DeconvolutedPeakSolution",
    "MzMLLoader", "MzXMLLoader", "MSFileLoader", "ScanProcessor",
    "deconvolute_peaks", 'version'
]
