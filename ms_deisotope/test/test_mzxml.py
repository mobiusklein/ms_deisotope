import unittest

from ms_deisotope.data_source import MzXMLLoader
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type


class TestMzXMLLoaderScanBehavior(unittest.TestCase):
    path = datafile("microscans.mzXML")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path)

    @property
    def first_scan(self):
        return self.reader.next().precursor

    def test_id(self):
        loader = self.reader
        scan = next(loader).precursor
        self.assertEqual(scan.id, "210")
        scan = loader.get_scan_by_id("210")
        self.assertEqual(scan.id, "210")

    def test_start_from_scan(self):
        loader = self.reader
        time = 0.4856916666666667
        bunch = next(loader.start_from_scan(rt=time))
        self.assertAlmostEqual(bunch.precursor.scan_time, time, 3)
        ix = bunch.precursor.index
        assert next(loader.start_from_scan(index=ix)).precursor == bunch.precursor

    def test_polarity(self):
        self.assertEqual(self.first_scan.polarity, 1)

    def test_index(self):
        self.assertEqual(self.first_scan.polarity, 1)

    def test_arrays(self):
        self.assertEqual(len(self.first_scan.arrays), 2)

    def test_precursor_info(self):
        self.assertEqual(self.first_scan.precursor_information, None)

    def test_file_description(self):
        file_info = self.reader.file_description()
        source_file = file_info.source_files[0]
        assert source_file.name == "AGP_tryptic_300ng_3microscans_glycoproteomics_nCE_27-35.raw"
        assert "location" not in source_file.parameters

    def test_data_processing(self):
        proc_info = self.reader.data_processing()
        assert len(proc_info) == 2


if __name__ == '__main__':
    unittest.main()
