import operator
from .averagine import Averagine, peptide, neutral_mass
from .peak_set import DeconvolutedPeak, DeconvolutedPeakSet
from .scoring import g_test, g_test_scaled, decon2ls_chisqr_test


class DeconvoluterBase(object):

    def match_theoretical_isotopic_distribution(self, theoretical_distribution, error_tolerance_ppm=2e-5):
        experimental_distribution = [self.peaklist.has_peak(
            p.mz, error_tolerance_ppm) for p in theoretical_distribution]
        experimental_distribution = [p if p is not None else FittedPeak(mz=None, intensity=1.0, 0, 0, 0, 0)
                                     for p in experimental_distribution]
        return experimental_distribution


class AveragineDeconvoluter(DeconvoluterBase):
    def __init__(self, peaklist, averagine=None, scorer=g_test_scaled):
        if averagine is None:
            averagine = peptide
        else:
            averagine = Averagine(averagine)
        self.peaklist = peaklist.clone()
        self.averagine = averagine
        self.scorer = scorer
        self._deconvoluted_peaks = []

    def charge_state_determination(self, peak, charge_range=(1, 8)):
        results = []
        for charge in range(*charge_range):
            tid = self.averagine.isotopic_cluster(peak.mz, charge)
            eid = self.match_theoretical_isotopic_distribution(tid)
            if len(eid) < 2:
                continue
            score = self.scorer(eid, tid)
            results.append((score, charge, tid, eid))
        return min(results, key=operator.itemgetter(0))

    def subtraction(self, isotopic_cluster):
        for peak in isotopic_cluster:
            match = self.peaklist.has_peak(peak.mz)
            if match is not None:
                match.intensity -= peak.intensity

    def deconvolute_peak(self, peak, charge_range=(1, 8)):
        score, charge, tid, eid = self.charge_state_determination(peak, charge_range=charge_range)
        total_abundance = sum(p.intensity for p in eid)
        monoisotopic_mass = neutral_mass(eid[0].mz, charge)
        peak = DeconvolutedPeak(monoisotopic_mass, total_abundance, charge, signal_to_noise=eid[0].signal_to_noise,
                                index=len(self._deconvoluted_peaks),
                                full_width_at_half_max=eid[0].full_width_at_half_max)
        self._deconvoluted_peaks.append(peak)
        self.subtraction(tid)

    def deconvolute(self, charge_range=(1, 8)):
        for peak in self.peaklist:
            self.deconvolute_peak(peak)
        return DeconvolutedPeakSet(self._deconvoluted_peaks)
