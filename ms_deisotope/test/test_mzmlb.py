import unittest
import os

from ms_deisotope.data_source.mzmlb import MzMLbLoader, determine_if_available
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type


@unittest.skipIf(not determine_if_available(), "mzMLb libraries not available")
class TestMzMLbLoaderScanBehavior(unittest.TestCase):
    path = datafile("20150710_3um_AGP_001_29_30.mzMLb")
    ref_path = datafile("20150710_3um_AGP_001_29_30.mzMLb")
    reader = None
    reference_reader = None

    @classmethod
    def setUpClass(cls):
        cls.reader = MzMLbLoader(cls.path)
        cls.reference_reader = infer_type.MSFileLoader(cls.ref_path)
        super(TestMzMLbLoaderScanBehavior, cls).setUpClass()

    @classmethod
    def tearDownClass(cls):
        cls.reader.close()
        cls.reference_reader.close()
        super(TestMzMLbLoaderScanBehavior, cls).tearDownClass()

    def test_infer(self):
        reader = infer_type.MSFileLoader(self.path)
        assert isinstance(reader, MzMLbLoader)

    def test_get_by_id_equiv(self):
        reader = self.reader
        reference_reader = self.reference_reader
        scan = reader.get_scan_by_id("scanId=1740226")
        ref = reference_reader.get_scan_by_id("scanId=1740226")
        assert scan == ref

    def test_start_from_equiv(self):
        reader = self.reader
        reference_reader = self.reference_reader

        n = len(reader)
        mid = reader[n // 2].scan_time
        reader.start_from_scan(rt=mid)
        reference_reader.start_from_scan(rt=mid)
        i = 0
        for a, b in zip(reader, reference_reader):
            assert a == b
            i += 1
            if i > 5:
                break

