import unittest


from ms_deisotope.data_source import infer_type
from .common import datafile


class TestInferType(unittest.TestCase):
    mzml_path = datafile("three_test_scans.mzML")
    mzmlgz_path = datafile("three_test_scans.mzML.gz")

    def test_infer_loader(self):
        reader = infer_type.MSFileLoader(self.mzml_path)
        self.assertIsNotNone(next(reader))

    def test_guess_from_path(self):
        loader_t = infer_type.guess_type_from_path(self.mzml_path)
        self.assertEqual(loader_t, infer_type.MzMLLoader)

    def test_guess_from_file_sniffing(self):
        loader_t = infer_type.guess_type_from_file_sniffing(self.mzml_path)
        self.assertEqual(loader_t, infer_type.MzMLLoader)

    def test_infer_loader_compressed(self):
        reader = infer_type.MSFileLoader(self.mzmlgz_path)
        assert reader.source_file_name is not None
        self.assertIsNotNone(next(reader))
