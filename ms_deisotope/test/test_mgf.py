import unittest

from ms_deisotope.data_source import MGFLoader, Scan
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type


class TestMGFLoaderScanBehavior(unittest.TestCase):
    path = datafile("test_utf_16.mgf")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path, encoding='utf-16')

    def test_index(self):
        reader = self.reader
        assert len(reader.index) == 287
        scan = reader.get_scan_by_id(
            '20150710_3um_AGP_001.1700581.1700581.4 File:"20150710_3um_AGP_001.d", NativeID:"scanId=1700581"')
        assert scan.id ==\
            '20150710_3um_AGP_001.1700581.1700581.4 File:"20150710_3um_AGP_001.d", NativeID:"scanId=1700581"'

    def test_scan_interface(self):
        reader = self.reader
        scan = next(reader)
        assert isinstance(scan, Scan)
        assert not scan.is_profile
        assert scan.precursor_information.precursor_scan_id is None


if __name__ == '__main__':
    unittest.main()
