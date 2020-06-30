'''
ms_deisotope
------------
A deisotoping and charge state deconvolution library for high resolution mass spectra
with integrated support for reading and writing common mass spectrometry data formats.
'''
try:
    import dill
except ImportError:
    pass

from .version import version

from .averagine import (
    Averagine, AveragineCache,
    peptide, glycan, glycopeptide, heparin,
    heparan_sulfate, permethylated_glycan,
    mass_charge_ratio, neutral_mass,
    calculate_mass, isotopic_shift)
from .deconvolution import (
    AveragineDeconvoluter, CompositionListDeconvoluter,
    AveraginePeakDependenceGraphDeconvoluter,
    CompositionListPeakDependenceGraphDeconvoluter,
    deconvolute_peaks)
from .scoring import MSDeconVFitter, PenalizedMSDeconVFitter, DistinctPatternFitter, IsotopicFitRecord
from .peak_set import (DeconvolutedPeak, DeconvolutedPeakSet, DeconvolutedPeakSolution, decharge)
from .processor import ScanProcessor, process
from .data_source import (MzMLLoader, MzXMLLoader, MSFileLoader, MGFLoader)


def get_include():
    """Get the include path for C extension headers.

    Returns
    -------
    :class:`str`
    """
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
    "deconvolute_peaks", 'version', "process", "decharge"
]


try:
    from ms_deisotope._c import deconvoluter_base as _cdeconvoluter_base
    has_c = True
    _has_c_error = None
except ImportError as _has_c_error:
    has_c = False


def check_c_extensions():
    if has_c:
        print("C extensions appear to have imported successfully.")
        return True
    else:
        print("Could not import deconvolution machinery: %r" % (_has_c_error, ))
    return False
