import operator
import numpy as np

from .averagine import Averagine, MultipleAveragine, peptide, multiple, neutral_mass, isotopic_variants
from .peak_set import DeconvolutedPeak, DeconvolutedPeakSet
from .scoring import g_test_scaled
from .utils import range, Base

from ms_peak_picker import FittedPeak


class DeconvoluterBase(object):
    use_subtraction = False

    def match_theoretical_isotopic_distribution(self, theoretical_distribution, error_tolerance=2e-5):
        experimental_distribution = [self.peaklist.has_peak(
            p.mz, error_tolerance) for p in theoretical_distribution]
        experimental_distribution = [p if p is not None else FittedPeak(-1., 1.0, 0, 0, 0, 0)
                                     for p in experimental_distribution]
        return experimental_distribution

    def scale_theoretical_distribution(self, theoretical_distribution, experimental_distribution):
        total_abundance = sum(p.intensity for p in experimental_distribution)
        for peak in theoretical_distribution:
            peak.intensity *= total_abundance
        return theoretical_distribution

    def subtraction(self, isotopic_cluster, error_tolerance=2e-5):
        if self.use_subtraction:
            for peak in isotopic_cluster:
                match = self.peaklist.has_peak(peak.mz, error_tolerance)
                if match is not None:
                    match.intensity -= peak.intensity
                    if match.intensity < 0:
                        match.intensity = 1.


class MinimizeDecider(Base):
    def __init__(self, threshold):
        self.threshold = threshold

    def __call__(self, score):
        if score <= self.threshold:
            return True


def charge_range_(lo, hi, step=None):
    sign = -1 if lo < 0 else 1
    abs_lo, abs_hi = abs(lo), abs(hi)
    upper = max(abs_lo, abs_hi)
    lower = min(abs_lo, abs_hi)

    for c in range(upper, lower - 1, -1):
        yield c * sign


class AveragineDeconvoluter(DeconvoluterBase):
    def __init__(self, peaklist, averagine=None, scorer=g_test_scaled, decider=MinimizeDecider(0.5),
                 use_subtraction=False):
        if averagine is None:
            averagine = peptide
        else:
            averagine = Averagine(averagine)
        if isinstance(decider, float):
            decider = MinimizeDecider(decider)
        self.peaklist = peaklist.clone()
        self.averagine = averagine
        self.scorer = scorer
        self.decider = decider
        self._deconvoluted_peaks = []
        self.use_subtraction = use_subtraction

    def charge_state_determination(self, peak, error_tolerance=2e-5, charge_range=(1, 8)):
        results = []
        for charge in charge_range_(*charge_range):
            tid = self.averagine.isotopic_cluster(peak.mz, charge)
            eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance)
            if len(eid) < 2 and len(tid) > 1:
                continue
            score = self.scorer(eid, tid, error_tolerance=error_tolerance)
            if np.isnan(score):
                continue

            results.append((score, charge, tid, eid))
        try:
            result = min(results, key=operator.itemgetter(0))
            return result
        except ValueError:
            return None

    def deconvolute_peak(self, peak, error_tolerance=2e-5, charge_range=(1, 8)):
        charge_det = self.charge_state_determination(peak, charge_range=charge_range)
        if charge_det is None:
            return
        score, charge, tid, eid = charge_det
        if self.decider(score):
            total_abundance = sum(p.intensity for p in eid)
            monoisotopic_mass = neutral_mass(eid[0].mz, charge)
            peak = DeconvolutedPeak(monoisotopic_mass, total_abundance, charge, signal_to_noise=eid[0].signal_to_noise,
                                    index=len(self._deconvoluted_peaks),
                                    full_width_at_half_max=eid[0].full_width_at_half_max, score=score)
            self._deconvoluted_peaks.append(peak)
            self.scale_theoretical_distribution(tid, eid)
            self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), order_chooser='signal_to_noise'):
        for peak in sorted(self.peaklist, key=operator.attrgetter(order_chooser), reverse=True):
            self.deconvolute_peak(peak, error_tolerance=error_tolerance,  charge_range=charge_range)
        return DeconvolutedPeakSet(self._deconvoluted_peaks)._reindex()


class CompositionListDeconvoluter(DeconvoluterBase):
    def __init__(self, peaklist, composition_list, scorer=g_test_scaled, decider=MinimizeDecider(0.5),
                 use_subtraction=False):
        if isinstance(decider, float):
            decider = MinimizeDecider(decider)
        self.peaklist = peaklist.clone()
        self.composition_list = composition_list
        self.scorer = scorer
        self.decider = decider
        self.use_subtraction = use_subtraction
        self._deconvoluted_peaks = []

    def deconvolute_composition(self, composition, error_tolerance=2e-5, charge_range=(1, 8)):
        for charge in charge_range_(*charge_range):
            tid = isotopic_variants(composition, charge=charge)
            eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance)

            if len(eid) < 2 and len(tid) > 1:
                continue
            score = self.scorer(eid, tid)
            if np.isnan(score):
                continue
            if self.decider(score):
                total_abundance = sum(p.intensity for p in eid)
                monoisotopic_mass = neutral_mass(eid[0].mz, charge)
                peak = DeconvolutedPeak(monoisotopic_mass, total_abundance, charge,
                                        signal_to_noise=eid[0].signal_to_noise,
                                        index=len(self._deconvoluted_peaks),
                                        full_width_at_half_max=eid[0].full_width_at_half_max, score=score)
                self._deconvoluted_peaks.append(peak)
                self.scale_theoretical_distribution(tid, eid)
                self.subtraction(tid, error_tolerance)

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8)):
        for composition in self.composition_list:
            self.deconvolute_composition(composition, error_tolerance=error_tolerance,
                                         charge_range=charge_range)
        return DeconvolutedPeakSet(self._deconvoluted_peaks)._reindex()


class MultipleAveragineDiagnosticFitter(DeconvoluterBase):
    def __init__(self, peaklist, averagine=None, scorer=g_test_scaled, decider=MinimizeDecider(0.5),
                 use_subtraction=False):
        if averagine is None:
            averagine = multiple
        else:
            averagine = MultipleAveragine(averagine)
        if isinstance(decider, float):
            decider = MinimizeDecider(decider)
        self.peaklist = peaklist.clone()
        self.averagine = averagine
        self.scorer = scorer
        self.decider = decider
        self._deconvoluted_peaks = []
        self.use_subtraction = use_subtraction

    def charge_state_determination(self, peak, error_tolerance=2e-5, charge_range=(1, 8)):
        totals = []
        for charge in charge_range_(*charge_range):
            results = []
            for averagine, tid in self.averagine.isotopic_cluster(peak.mz, charge):
                eid = self.match_theoretical_isotopic_distribution(tid, error_tolerance)
                if len(eid) < 2 and len(tid) > 1:
                    results.append(float('inf'), charge, tid, eid)
                    continue
                score = self.scorer(eid, tid, error_tolerance=error_tolerance)
                if np.isnan(score):
                    results.append(float('inf'), charge, tid, eid)
                    continue

                results.append((score, charge, tid, eid))
            totals.append(results)
        return totals

    def deconvolute(self, error_tolerance=2e-5, charge_range=(1, 8), order_chooser='signal_to_noise'):
        for peak in sorted(self.peaklist, key=operator.attrgetter(order_chooser), reverse=True):
            scores = self.charge_state_determination(peak, error_tolerance=error_tolerance,  charge_range=charge_range)
            self._deconvoluted_peaks.append(scores)
        return self._deconvoluted_peaks


def deconvolute_peaks(peaklist, deconvoluter_type=AveragineDeconvoluter, decon_config=None, charge_range=(1, 8),
                      error_tolerance=2e-5):
    decon_config = decon_config or {}
    decon = deconvoluter_type(peaklist=peaklist, **decon_config)
    deconvoluted_peaks = decon.deconvolute
    return deconvoluted_peaks
