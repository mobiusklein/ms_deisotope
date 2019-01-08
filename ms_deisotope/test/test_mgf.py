import unittest

from ms_deisotope.data_source import MGFLoader, Scan
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type


class TestMGFLoaderScanBehavior(unittest.TestCase):
    path = datafile("small.mgf")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path)

    def test_index(self):
        reader = self.reader
        assert len(reader.index) == 34
        scan = reader.get_scan_by_id(
            'small.10.10')
        assert scan.id ==\
            'small.10.10'

    def test_scan_interface(self):
        reader = self.reader
        scan = next(reader)
        assert isinstance(scan, Scan)
        assert not scan.is_profile
        assert scan.precursor_information.precursor_scan_id is None


if __name__ == '__main__':
    unittest.main()
