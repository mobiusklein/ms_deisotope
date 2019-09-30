import unittest

from ms_deisotope.data_source import MGFLoader, Scan
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type


class TestMGFLoaderScanBehavior(unittest.TestCase):
    path = datafile("small.mgf")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path)

    def test_source_file_name(self):
        reader = self.reader
        assert reader.source_file_name.endswith("small.mgf")

    def test_index(self):
        reader = self.reader
        assert len(reader.index) == 34
        scan = reader.get_scan_by_id(
            'small.10.10')
        assert scan.id ==\
            'small.10.10'
        scan = reader[10]
        assert scan.index == 10

    def test_get_time(self):
        reader = self.reader
        scan = reader.get_scan_by_time(0.3)
        assert scan.id == 'small.31.31'
        scan = next(reader.start_from_scan(rt=0.3, grouped=False))
        assert scan.id == 'small.31.31'

    def test_annotations(self):
        scan = self.reader[10]
        assert scan.annotations == {}


    def test_scan_interface(self):
        reader = self.reader
        scan = next(reader)
        assert isinstance(scan, Scan)
        assert not scan.is_profile
        assert scan.precursor_information.precursor_scan_id is None


if __name__ == '__main__':
    unittest.main()
