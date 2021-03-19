import unittest

from ms_deisotope.data_source import memory
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type

scan_ids = [
    "controllerType=0 controllerNumber=1 scan=10014",
    "controllerType=0 controllerNumber=1 scan=10015",
    "controllerType=0 controllerNumber=1 scan=10016"
]


class TestMemoryScanSource(unittest.TestCase):
    path = datafile("three_test_scans.mzML")

    @property
    def source_reader(self):
        return infer_type.MSFileLoader(self.path)

    @property
    def prepare_source(self):
        source = self.source_reader
        loader = memory.ScanCollection.build(source)
        return loader

    def test_iteration(self):
        g = iter(scan_ids)
        bunch = next(self.prepare_source)
        assert bunch.precursor.id == next(g)
        for product in bunch.products:
            assert product.id == next(g)

    def test_source_file_name_none(self):
        source = self.prepare_source
        assert source.source_file_name is None


if __name__ == '__main__':
    unittest.main()
