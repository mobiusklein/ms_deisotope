import unittest

import pytest

from ms_deisotope.output import ProcessedMzMLDeserializer
from ms_deisotope.data_source import get_opener

from ms_deisotope.clustering import scan_clustering

from ms_deisotope.test.common import datafile


@pytest.mark.slow
class TestScanClustering(unittest.TestCase):
    path = datafile("AGP_tryptic_300ng_2microscans_glycoproteomics_nCE_27-30.preprocessed.mzML.gz")

    @property
    def reader(self):
        reader = ProcessedMzMLDeserializer(get_opener(self.path))
        return reader

    def load_msms_scans(self, reader):
        products = list(map(reader.get_scan_by_id, reader.extended_index.msn_ids.keys()))
        return products

    def cluster_scans(self, scans):
        clusters = scan_clustering.cluster_scans(scans)
        return clusters

    def test_cluster_scans(self):
        reader = self.reader
        scans = self.load_msms_scans(reader)
        clusters = self.cluster_scans(scans)
        assert len(clusters) == 1124


if __name__ == "__main__":
    unittest.main()
