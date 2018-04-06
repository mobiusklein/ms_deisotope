import unittest
import platform

import numpy as np

from ms_deisotope.data_source.thermo_raw import (
    ThermoRawLoader, filter_line_parser)
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type


@unittest.skipIf(not platform.platform().lower().startswith("windows"), "Requires Windows COM")
class TestThermoRawLoaderScanBehavior(unittest.TestCase):
    path = datafile("small.RAW")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path)

    def test_iteration(self):
        reader = self.reader
        reader.start_from_scan('controllerType=0 controllerNumber=1 scan=10')
        bunch = next(reader)
        assert bunch.precursor.id == 'controllerType=0 controllerNumber=1 scan=9'
        bunch = next(reader)
        assert bunch.precursor.id == 'controllerType=0 controllerNumber=1 scan=15'

    def test_file_level_metadata(self):
        reader = self.reader
        desc = reader.file_description()
        assert desc.has_content("MS1 spectrum")
        assert desc.has_content("MSn spectrum")

        inst_config = reader.instrument_configuration()
        assert inst_config[0].analyzers[0] == 'orbitrap'

    def test_scan_level_data(self):
        reader = self.reader
        reader.start_from_scan('controllerType=0 controllerNumber=1 scan=10')
        bunch = next(reader)
        assert np.isclose(bunch.precursor.scan_time, 0.077788333333)
        assert len(bunch.precursor.pick_peaks(signal_to_noise_threshold=1.5).peak_set) == 3110
        scan_window = bunch.precursor.acquisition_information.scan_list[0][0]
        assert scan_window.lower == 200.0 and scan_window.upper == 2000.0
        product = bunch.products[0]
        assert product.ms_level == 2
        assert product.index == 9
        assert product.activation.energy == 35.0
        assert np.isclose(product.precursor_information.mz, 810.7528)
        annotations = {'[Thermo Trailer Extra]Micro Scan Count': 3.0,
                       '[Thermo Trailer Extra]Scan Event': 3.0,
                       '[Thermo Trailer Extra]Scan Segment': 1.0,
                       'filter_line': 'ITMS + c ESI d Full ms2 810.75@cid35.00 [210.00-1635.00]'}
        assert product.annotations == annotations
        assert np.isclose(product.isolation_window.target, 810.752807)
        assert product.isolation_window.lower == 1.0
        assert not product.is_profile


if __name__ == '__main__':
    unittest.main()
