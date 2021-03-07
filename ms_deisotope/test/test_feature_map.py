import unittest

import ms_deisotope
from ms_deisotope.feature_map import feature_map
from ms_deisotope.test.common import datafile

complex_compressed_mzml = datafile("20150710_3um_AGP_001_29_30.mzML.gz")


class LCMSFeatureMapTest(unittest.TestCase):

    features = None

    @classmethod
    def setUpClass(cls):
        reader = ms_deisotope.MSFileLoader(complex_compressed_mzml)
        features = feature_map.LCMSFeatureForest.from_reader(reader)
        cls.features = features

    @classmethod
    def tearDownClass(cls):
        cls.features = None

    def test_forest(self):
        features = self.features
        assert len(features) == 4167
        f = features.search(1161.50875)
        assert f is not None

    def test_search(self):
        features = self.features
        f = features.search(1161.50875)
        assert f is not None
        assert abs(f.mz - 1161.50875) / f.mz < 1e-6

    def test_between(self):
        features = self.features
        between = (features.between(920, 940))
        assert len(between) == 180
