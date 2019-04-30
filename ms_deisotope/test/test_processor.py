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
            for product in scan_bunch.products:
                if product.precursor_information.defaulted:
                    candidates = scan_bunch.precursor.peak_set.between(
                        product.precursor_information.mz - 1, product.precursor_information.mz + 1)
                    assert len(candidates) == 0

    def test_complex_processor(self):
        proc = processor.ScanProcessor(self.complex_compressed_mzml, ms1_deconvolution_args={
            "averagine": glycopeptide,
            "scorer": PenalizedMSDeconVFitter(20., 2.),
            "truncate_after": 0.95
        }, msn_deconvolution_args={
            "averagine": peptide,
            "scorer": MSDeconVFitter(10.),
            "truncate_after": 0.8
        })
        bunch = next(proc)
        assert len(bunch.products) == 5
        for product in bunch.products:
            assert not product.precursor_information.defaulted
        recalculated_precursors = {
            'scanId=1740086': 4640.00074242012,
            'scanId=1740149': 4786.05878475792,
            'scanId=1740226': 4640.007868154431,
            'scanId=1740344': 4348.90894554512,
            'scanId=1740492': 5005.1329902247435
        }
        for product in bunch.products:
            mass = product.precursor_information.extracted_neutral_mass
            self.assertAlmostEqual(mass, recalculated_precursors[product.id], 2)

        proc.start_from_scan("scanId=1760847")
        bunch = next(proc)
        recalculated_precursors = {
            'scanId=1761168': 4640.01972225792,
            'scanId=1761235': 4640.019285920238,
            'scanId=1761325': 4786.07251976387,
            'scanId=1761523': 4696.016295197582,
            'scanId=1761804': 986.58798612896
        }
        for product in bunch.products:
            mass = product.precursor_information.extracted_neutral_mass
            self.assertAlmostEqual(mass, recalculated_precursors[product.id], 2)


if __name__ == '__main__':
    unittest.main()
