import unittest
import platform

import numpy as np

from ms_deisotope.data_source.thermo_raw import (
    determine_if_available)
from ms_deisotope.data_source.thermo_raw_net import (
    determine_if_available as determine_if_available_net)
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type

not_windows = not platform.platform().lower().startswith("windows")
missing_reader_dll = not determine_if_available()
missing_net_dll = not determine_if_available_net()


class ThermoRawLoaderScanBehaviorBase(object):
    path = datafile("small.RAW")

    reference_mzml = datafile("small.mzML")
    reference_mgf = datafile("small.mgf")

    @property
    def reader(self):  # pragma: no cover
        # Required for implementation
        raise NotImplementedError()

    def test_iteration(self):
        reader = self.reader
        reader.start_from_scan('controllerType=0 controllerNumber=1 scan=10')
        bunch = next(reader)
        assert bunch.precursor.id == 'controllerType=0 controllerNumber=1 scan=9'
        bunch = next(reader)
        assert bunch.precursor.id == 'controllerType=0 controllerNumber=1 scan=15'
        reader.start_from_scan(rt=0.077788333333)
        bunch = next(reader)
        assert np.isclose(bunch.precursor.scan_time, 0.077788333333)

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
                       'filter string': 'ITMS + c ESI d Full ms2 810.75@cid35.00 [210.00-1635.00]'}
        assert product.annotations == annotations
        assert np.isclose(product.isolation_window.target, 810.752807)
        assert product.isolation_window.lower == 1.0
        assert not product.is_profile

    def test_size(self):
        reader = self.reader
        n = len(reader)
        assert n == 48
        x = reader[-1]
        y = reader.get_scan_by_time(float('inf'))
        assert x == y

    def test_compat(self):
        raw_reader = self.reader
        mzml_reader = infer_type.MSFileLoader(self.reference_mzml)
        mgf_reader = infer_type.MSFileLoader(self.reference_mgf)

        mgf_scan = next(mgf_reader)
        mzml_scan = mzml_reader[2]
        raw_scan = raw_reader[2]

        self.assertEqual(mzml_scan, raw_scan)

        mgf_scan.pick_peaks()
        raw_scan.pick_peaks()
        self.assertEqual(raw_scan.peak_set, mgf_scan.peak_set)

        self.assertEqual(raw_scan.precursor_information.precursor,
                         mzml_scan.precursor_information.precursor)

    def test_samples(self):
        reader = self.reader
        samples = reader.samples()
        assert len(samples) == 1
        sample = samples[0]
        assert sample.name == '1'
        assert sample.id == '1'
        assert sample.parameters['sample vial'] == '1a1'

    def test_source_file_name(self):
        reader = self.reader
        assert reader.source_file_name.lower().endswith("small.raw")


@unittest.skipIf(not_windows or missing_reader_dll, "Requires Windows COM and MSFileReader.dll")
class TestCOMThermoRawLoaderScanBehavior(ThermoRawLoaderScanBehaviorBase, unittest.TestCase):

    @property
    def reader(self):
        from ms_deisotope.data_source.thermo_raw import ThermoRawLoader
        return ThermoRawLoader(self.path)


@unittest.skipIf(missing_net_dll, "Requires .NET Runtime and RawFileReader.dll")
class TestNETThermoRawLoaderScanBehavior(ThermoRawLoaderScanBehaviorBase, unittest.TestCase):

    @property
    def reader(self):
        from ms_deisotope.data_source.thermo_raw_net import ThermoRawLoader
        return ThermoRawLoader(self.path)


if __name__ == '__main__':
    unittest.main()
