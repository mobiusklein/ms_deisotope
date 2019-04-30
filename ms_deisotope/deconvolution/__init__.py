# -*- coding: utf-8 -*-
'''This module defines a collection of isotopic envelope search strategies for
deisotoping and charge state deconvolution. Each strategy is implemented through
a subtype of :class:`~.DeconvoluterBase` and provide a common set of methods for
deconvolving a peak list.
'''
# from .deconvolution import *

from .utils import (
    prepare_peaklist,
    charge_range_,
    ChargeIterator,
    quick_charge,
    drop_placeholders_parallel,
    drop_placeholders,
    first_peak,
    from_fitted_peak,
    count_placeholders,
    mean)

from .base import (
    DeconvoluterBase)

from .exhaustive import (
    ExhaustivePeakSearchDeconvoluterBase,
    PeakDependenceGraphDeconvoluterBase)

from .averagine_based import (
    AveragineDeconvoluterBase,
    AveragineDeconvoluter,
    AveraginePeakDependenceGraphDeconvoluter,
    MultiAveragineDeconvoluterBase,
    MultiAveragineDeconvoluter,
    MultiAveraginePeakDependenceGraphDeconvoluter)

from .composition_list import (
    CompositionListDeconvoluterBase,
    CompositionListDeconvoluter,
    CompositionListPeakDependenceGraphDeconvoluter)

from .api import deconvolute_peaks
