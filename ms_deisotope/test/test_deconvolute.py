import unittest

import numpy as np

import brainpy

from ms_deisotope.data_source import common, mzml, MSFileLoader
from ms_deisotope.averagine import peptide, glycopeptide, TheoreticalIsotopicPattern
from ms_peak_picker import reprofile
from ms_deisotope.deconvolution import (
    deconvolute_peaks, AveragineDeconvoluter,
    AveraginePeakDependenceGraphDeconvoluter,
    CompositionListDeconvoluter,
    CompositionListPeakDependenceGraphDeconvoluter)
from ms_deisotope.scoring import PenalizedMSDeconVFitter
from brainpy import neutral_mass

from ms_deisotope.test.test_scan import (
    make_profile, points, fwhm,
    FittedPeak)
from ms_deisotope.test.common import datafile


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
        scan.pick_peaks()
        for point in points:
            self.assertIsNotNone(scan.peak_set.has_peak(point[0]))
        return scan

    def test_deconvolution(self):
        scan = self.make_scan()
        algorithm_type = AveragineDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, {
                "averagine": peptide,
                "scorer": PenalizedMSDeconVFitter(5., 1.),
                "use_subtraction": False
            }, deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        assert len(dpeaks) == 6
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
        assert len(dpeaks) == 2
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
                tid.scale_raw(abundance * 100)
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


class TestSolutionRetrieval(unittest.TestCase):
    def make_scan(self):
        complex_compressed_mzml = datafile("20150710_3um_AGP_001_29_30.mzML.gz")
        reader = MSFileLoader(complex_compressed_mzml)
        bunch = next(reader)
        return bunch

    def target_envelopes(self):
        envelopes = [
            [1161.00927734375,
             1161.25830078125,
             1161.508544921875,
             1161.7586669921875,
             1162.0111083984375,
             1162.2607421875,
             1162.5128173828125],
            [1197.52197265625,
             1197.7735595703125,
             1198.024169921875,
             1198.2742919921875,
             1198.524658203125,
             1198.774169921875,
             1199.02490234375]
            [929.0088500976562,
             929.2091674804688,
             929.4107055664062,
             929.6107788085938,
             929.8109741210938,
             930.01171875,
             930.213134765625]
            [1088.234619140625,
             1088.4853515625,
             1088.7354736328125,
             1088.987060546875,
             1089.237060546875,
             1089.4873046875,
             1089.7379150390625],
            [1002.0338745117188,
             1002.2335815429688,
             1002.4354858398438,
             1002.63671875,
             1002.837646484375,
             1003.0364990234375,
             1003.2337646484375],
        ]
        return envelopes

    def test_retrieve_deconvolution_solution(self):
        bunch = self.make_scan()
        scan = bunch.precursor
        scan.pick_peaks()
        ms1_deconvolution_args = {
            "averagine": glycopeptide,
            "scorer": PenalizedMSDeconVFitter(20., 2.),
        }
        priorities = []
        for product in bunch.products:
            priorities.append(scan.has_peak(product.precursor_information.mz))
        algorithm_type = AveraginePeakDependenceGraphDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, ms1_deconvolution_args,
            priority_list=priorities,
            deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        deconvoluter = deconresult.deconvoluter
        priority_results = deconresult.priorities
        reference_deconvoluter = algorithm_type(scan.peak_set.clone(), **ms1_deconvolution_args)
        for i, result in enumerate(priority_results):
            query = priorities[i].mz
            if result is None:
                raw_peaks = scan.peak_set.between(query - 2, query + 3)
                anchor_peak = scan.peak_set.has_peak(query)
                deconvoluted_peaks = dpeaks.between(query - 2, query + 3, use_mz=True)
                possible_solutions = reference_deconvoluter._fit_all_charge_states(anchor_peak)
                assert anchor_peak is not None and raw_peaks and possible_solutions and not deconvoluted_peaks
                assert deconvoluter.peak_dependency_network.find_solution_for(anchor_peak) is not None
                assert dpeaks.has_peak(query, use_mz=True)
            else:
                assert 0 <= abs(result.mz - query) < 1
                anchor_peak = scan.peak_set.has_peak(query)
                assert deconvoluter.peak_dependency_network.find_solution_for(anchor_peak) is not None


if __name__ == '__main__':
    unittest.main()
