import unittest

from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type
from ms_deisotope.data_source.metadata import scan_traits


class TestScanTraits(unittest.TestCase):
    path = datafile("three_test_scans.mzML")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path)

    def test_traits(self):
        bunch = next(self.reader)
        scan = bunch.precursor
        acquisition = scan.acquisition_information
        assert len(acquisition) == 1
        scan_event = acquisition[0]
        assert not scan_event.has_ion_mobility()
        assert len(scan_event) == 1
        scan_window = scan_event[0]
        assert scan_window.lower == 350
        assert scan_window.upper == 1500
        assert not scan_window.is_empty()
        assert scan_traits.ScanWindow(0, 0).is_empty()
        assert scan_traits.ScanWindow(350, 1500) == scan_window
        assert scan_window == scan_event.total_scan_window()
        assert scan_event == scan_event
        assert acquisition == acquisition
        assert scan.tic.raw() - 1.8161617e+10 + 104 == 0.0
        scan.pick_peaks()
        assert abs(scan.tic.centroided() - 4531158140.658203) < 1e-3
        assert abs(scan.tic() - 4531158140.658203) < 1e-3

        scan = bunch.products[0]
        isolation = scan.isolation_window
        assert not isolation.is_empty()
        assert scan_traits.IsolationWindow(None, 200, None).is_empty()
        assert scan_traits.IsolationWindow(0, 200, 0).is_empty()
        assert isolation.spans(isolation.target)
        assert isolation == isolation


