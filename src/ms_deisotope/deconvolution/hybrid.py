from typing import List, Tuple
from ms_deisotope.constants import (ERROR_TOLERANCE, TRUNCATE_AFTER, IGNORE_BELOW)
from ms_deisotope.averagine import PROTON

from .utils import charge_range_, drop_placeholders
from .averagine_based import AveraginePeakDependenceGraphDeconvoluter
from .composition_list import CompositionListDeconvoluterBase, CompositionType

_APDGD = AveraginePeakDependenceGraphDeconvoluter


class HybridAveragineCompositionListPeakDependenceGraphDeconvoluter(_APDGD, CompositionListDeconvoluterBase):

    composition_list: List[CompositionType]

    def __init__(self, peaklist, composition_list=None, *args, **kwargs):
        if composition_list is None:
            composition_list = []
        super(HybridAveragineCompositionListPeakDependenceGraphDeconvoluter, self).__init__(
            peaklist, composition_list=composition_list, *args, **kwargs)

    def deconvolute_composition(self, composition: CompositionType,
                                error_tolerance: float = ERROR_TOLERANCE,
                                charge_range: Tuple[int, int] = (1, 8),
                                charge_carrier: float = PROTON,
                                truncate_after: float = TRUNCATE_AFTER,
                                ignore_below: float = IGNORE_BELOW,
                                mass_shift: float = None):
        for charge in charge_range_(*charge_range):
            fit = self.fit_composition_at_charge(
                composition, charge, error_tolerance, charge_carrier=charge_carrier,
                truncate_after=truncate_after, mass_shift=mass_shift, ignore_below=ignore_below)
            if fit is None:
                continue
            rep_eid = drop_placeholders(fit.experimental)
            if len(rep_eid) == 1 and fit.charge > 1:
                continue
            if not self.scorer.reject(fit):
                fit.data = (composition, charge)
                self.peak_dependency_network.add_fit_dependence(fit)
            if self.incremental_truncation is not None:
                for case in self.fit_incremental_truncation(fit, self.incremental_truncation):
                    if not self.scorer.reject(case):
                        self.peak_dependency_network.add_fit_dependence(case)

    def populate_graph(self, error_tolerance: float=ERROR_TOLERANCE,
                       charge_range: Tuple[int, int]=(1, 8),
                       left_search_limit: int=1,
                       right_search_limit: int=0,
                       charge_carrier: float=PROTON,
                       truncate_after: float=TRUNCATE_AFTER,
                       ignore_below: float=IGNORE_BELOW,
                       mass_shift: float=None):
        for composition in self.composition_list:
            self.deconvolute_composition(
                composition, error_tolerance, charge_range=charge_range,
                truncate_after=truncate_after,
                charge_carrier=charge_carrier,
                ignore_below=ignore_below,
                mass_shift=mass_shift)
        AveraginePeakDependenceGraphDeconvoluter.populate_graph(
            self, error_tolerance,
            charge_range=charge_range,
            left_search_limit=left_search_limit,
            right_search_limit=right_search_limit,
            charge_carrier=charge_carrier,
            truncate_after=truncate_after,
            ignore_below=ignore_below)
