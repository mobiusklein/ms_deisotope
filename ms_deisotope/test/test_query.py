import unittest

import ms_deisotope
from ms_deisotope.data_source import query
from ms_deisotope.test.common import datafile



class TestFancyIterator(unittest.TestCase):
    complex_compressed_mzml = datafile("20150710_3um_AGP_001_29_30.mzML.gz")

    def _get_reader(self):
        return ms_deisotope.MSFileLoader(self.complex_compressed_mzml)

    def test_time_interval_iterator(self):
        reader = self._get_reader()
        tii = query.TimeIntervalIterator(reader, 29.5, 31)
        assert tii.has_ms1_scans()
        assert tii.has_msn_scans()
        n = 0
        for precursor, products in tii:
            if n == 0:
                assert abs(precursor.scan_time - 29.5) < 1e-2
            else:
                assert precursor.scan_time >= 29.5
            n += 1
        assert n == 24

    def test_index_interval_iterator(self):
        reader = self._get_reader()
        iii = query.IndexIntervalIterator(reader, end=31)
        assert iii.has_ms1_scans()
        assert iii.has_msn_scans()
        assert iii.start == 0

        n = 0
        for precursor, products in iii:
            n += 1
            assert precursor.index <= 31
        assert n == 6

    def test_ms_level_filter(self):
        reader = self._get_reader()
        flt = query.MSLevelFilter(reader, 2)
        n = 0
        for batch in flt:
            assert batch.precursor is None
            n += 1
        assert n == 52

    def test_ms1_merger(self):
        reader = self._get_reader()
        trf = query.MS1MergingTransformer(reader)
        n = 0
        n2 = 0
        for batch in trf:
            n += 1
            n2 += len(batch.products)
        assert n == 10
        assert n2 == 260
