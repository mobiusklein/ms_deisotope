import unittest

import numpy as np

import brainpy

from ms_deisotope.data_source import common, mzml
from ms_deisotope.averagine import peptide, TheoreticalIsotopicPattern
from ms_peak_picker import reprofile
from ms_deisotope.deconvolution import (
    deconvolute_peaks, AveragineDeconvoluter,
    AveraginePeakDependenceGraphDeconvoluter,
    CompositionListDeconvoluter,
    CompositionListPeakDependenceGraphDeconvoluter)
from ms_deisotope.scoring import PenalizedMSDeconVFitter
from brainpy import neutral_mass
from ms_deisotope.test.test_scan import (
    make_profile, points, fwhm, gaussian_shape,
    FittedPeak)


class TestAveragineDeconvolution(unittest.TestCase):
    def make_scan(self):
        mz, intensity = make_profile(points, fwhm)
        scan = common.Scan(
            {
                "id": "test-scan",
                "index": 0,
                "m/z array": mz,
                "ms level": 1,
                "scan time": 0.0,
                "intensity array": intensity,
                "profile spectrum": "",
                "positive scan": "",
            },
            mzml.MzMLDataInterface())
        return scan

    def test_deconvolution(self):
        scan = self.make_scan()
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        algorithm_type = AveragineDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, {
                "averagine": peptide,
                "scorer": PenalizedMSDeconVFitter(5., 1.),
                "use_subtraction": False
            }, deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        for point in points:
            peak = dpeaks.has_peak(neutral_mass(point[0], point[1]))
            self.assertIsNotNone(peak)

    def test_graph_deconvolution(self):
        scan = self.make_scan()
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        algorithm_type = AveraginePeakDependenceGraphDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, {
                "averagine": peptide,
                "scorer": PenalizedMSDeconVFitter(5., 1.)
            }, deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        deconvoluter = deconresult.deconvoluter
        for point in points:
            peak = dpeaks.has_peak(neutral_mass(point[0], point[1]))
            self.assertIsNotNone(peak)
            fp = scan.has_peak(peak.mz)
            self.assertAlmostEqual(
                deconvoluter.peak_dependency_network.find_solution_for(fp).mz,
                peak.mz, 3)


class TestCompositionListDeconvolution(unittest.TestCase):
    compositions = [
        brainpy.parse_formula('C84H138N6O62'),
        brainpy.parse_formula('C90H148N6O66'),
        brainpy.parse_formula('C90H148N6O69S1'),
        brainpy.parse_formula('C62H94N4O69S6'),
    ]

    charges = [
        ((2, 100), (3, 250),),
        ((2, 75), (3, 200),),
        ((2, 90), (3, 270),),
        ((3, 20), (4, 90), (5, 45))
    ]

    def make_scan(self):
        peaks = []
        for comp, charges in zip(self.compositions, self.charges):
            for charge, abundance in charges:
                tid = brainpy.isotopic_variants(comp, charge=-charge)
                tid = TheoreticalIsotopicPattern(tid, tid[0].mz)
                tid._scale_raw(abundance * 100)
                peaks.extend(tid)
        peaks.sort(key=lambda x: x.mz)

        mz = np.array([0])
        intensity = np.array([0])

        fpeaks = []
        for p in peaks:
            fpeaks.append(FittedPeak(
                mz=p.mz, intensity=p.intensity, signal_to_noise=p.intensity, peak_count=-1, index=-1,
                full_width_at_half_max=fwhm, area=p.intensity))
        mz, intensity = reprofile(fpeaks)
        scan = common.Scan(
            {
                "id": "test-scan",
                "index": 0,
                "m/z array": mz,
                "intensity array": intensity,
                "ms level": 1,
                "scan time": 0.0,
                "profile spectrum": "",
                "negative scan": "",
            },
            mzml.MzMLDataInterface())
        return scan

    def test_deconvolution(self):
        scan = self.make_scan()
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        algorithm_type = CompositionListDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, {
                "composition_list": self.compositions,
                "scorer": PenalizedMSDeconVFitter(5., 2.),
                "use_subtraction": True
            }, charge_range=(-1, -8), deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        n_cases = sum(map(len, self.charges))
        assert len(dpeaks) == n_cases

    def test_graph_deconvolution(self):
        scan = self.make_scan()
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        algorithm_type = CompositionListPeakDependenceGraphDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, {
                "composition_list": self.compositions,
                "scorer": PenalizedMSDeconVFitter(5., 2.),
                "use_subtraction": True
            }, charge_range=(-1, -8), deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        n_cases = sum(map(len, self.charges))
        assert len(dpeaks) == n_cases


if __name__ == '__main__':
    unittest.main()
