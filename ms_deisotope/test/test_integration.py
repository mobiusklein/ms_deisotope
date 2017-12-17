import unittest
import gzip
import pickle

from ms_deisotope.test.common import datafile

from ms_deisotope.data_source import common, mzml
from ms_deisotope.averagine import peptide
from ms_deisotope.scoring import MSDeconVFitter


class TestIntegration(unittest.TestCase):
    def get_scan(self):
        scan_data = pickle.load(gzip.open(datafile("test_scan.pkl.gz")))
        scan = common.Scan(scan_data, mzml.MzMLDataInterface())
        return scan

    def get_reference(self):
        processed_scan = pickle.load(gzip.open(datafile("test_scan_results.pkl.gz")))
        return processed_scan

    def test_scan(self):
        scan = self.get_scan()
        scan.pick_peaks()
        scan.deconvolute(averagine=peptide, scorer=MSDeconVFitter(10.0),
                         ignore_below=0.05, truncate_after=0.8, charge_range=(1, 4))
        reference = self.get_reference()
        for peak in scan.deconvoluted_peak_set:
            peaks = reference.deconvoluted_peak_set.all_peaks_for(peak.neutral_mass, 1e-6)
            for ref in peaks:
                if ref.charge == peak.charge:
                    self.assertTrue(
                        abs(peak.neutral_mass - ref.neutral_mass) < 1e-3,
                        peak.neutral_mass - ref.neutral_mass)
                    self.assertTrue(
                        abs(peak.intensity - ref.intensity) < 1e-3,
                        peak.intensity - ref.intensity)
                    self.assertTrue(
                        abs(peak.score - ref.score) < 1e-3,
                        peak.score - ref.score)
                    for at, bt in zip(peak.envelope, ref.envelope):
                        self.assertTrue(
                            abs(at.mz - bt.mz) < 1e-3,
                            at.mz - bt.mz)
                        self.assertTrue(
                            abs(at.intensity - bt.intensity) < 1e-3,
                            at.intensity - bt.intensity)


if __name__ == '__main__':
    unittest.main()
