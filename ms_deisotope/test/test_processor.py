import unittest

from ms_deisotope import processor
from ms_deisotope.averagine import glycopeptide, peptide
from ms_deisotope.scoring import PenalizedMSDeconVFitter, MSDeconVFitter

from ms_deisotope.test.common import datafile


class TestScanProcessor(unittest.TestCase):
    mzml_path = datafile("three_test_scans.mzML")
    missing_charge_mzml = datafile("has_missing_charge_state_info.mzML")
    complex_compressed_mzml = datafile("20150710_3um_AGP_001_29_30.mzML.gz")

    def test_processor(self):
        proc = processor.ScanProcessor(self.mzml_path, ms1_deconvolution_args={
            "averagine": glycopeptide,
            "scorer": PenalizedMSDeconVFitter(5., 2.)
        })
        for scan_bunch in iter(proc):
            self.assertIsNotNone(scan_bunch)
            self.assertIsNotNone(scan_bunch.precursor)
            self.assertIsNotNone(scan_bunch.products)

    def test_averaging_processor(self):
        proc = processor.ScanProcessor(self.mzml_path, ms1_deconvolution_args={
            "averagine": glycopeptide,
            "scorer": PenalizedMSDeconVFitter(5., 2.)
        }, ms1_averaging=1)
        for scan_bunch in iter(proc):
            self.assertIsNotNone(scan_bunch)
            self.assertIsNotNone(scan_bunch.precursor)
            self.assertIsNotNone(scan_bunch.products)

    def test_missing_charge_processing(self):
        proc = processor.ScanProcessor(self.missing_charge_mzml, ms1_deconvolution_args={
            "averagine": glycopeptide,
            "scorer": PenalizedMSDeconVFitter(5., 2.)
        })
        for scan_bunch in iter(proc):
            self.assertIsNotNone(scan_bunch)
            self.assertIsNotNone(scan_bunch.precursor)
            self.assertIsNotNone(scan_bunch.products)

    def test_complex_processor(self):
        proc = processor.ScanProcessor(self.complex_compressed_mzml, ms1_deconvolution_args={
            "averagine": glycopeptide,
            "scorer": PenalizedMSDeconVFitter(20., 2.),
        }, msn_deconvolution_args={
            "averagine": peptide,
            "scorer": MSDeconVFitter(10.),
        })
        bunch = next(proc)
        assert len(bunch.products) == 5
        for product in bunch.products:
            assert not product.precursor_information.defaulted
        product = bunch.products[0]
        diff_mz = product.precursor_information.mz - product.precursor_information.extracted_mz
        assert 0 < abs(diff_mz) < 1


if __name__ == '__main__':
    unittest.main()
