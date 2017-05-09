import unittest

import numpy as np

from ms_deisotope.data_source import common, mzml
from ms_deisotope.averagine import peptide
from ms_deisotope.deconvolution import (
    deconvolute_peaks, AveragineDeconvoluter,
    AveraginePeakDependenceGraphDeconvoluter)
from ms_deisotope.scoring import PenalizedMSDeconVFitter
from brainpy import neutral_mass
from ms_deisotope.test.test_scan import make_profile, points, fwhm


class TestDeconvolution(unittest.TestCase):
    def make_scan(self):
        mz, intensity = make_profile(points, fwhm)
        scan = common.Scan(
            {
                "m/z array": mz,
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
                "scorer": PenalizedMSDeconVFitter(5., 1.)
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


if __name__ == '__main__':
    unittest.main()
