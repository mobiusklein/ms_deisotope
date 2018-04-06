import unittest

from ms_deisotope.data_source import MzMLLoader
from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type, xml_reader

scan_ids = [
    "controllerType=0 controllerNumber=1 scan=10014",
    "controllerType=0 controllerNumber=1 scan=10015",
    "controllerType=0 controllerNumber=1 scan=10016"
]


class TestMzMLLoaderScanBehavior(unittest.TestCase):
    path = datafile("three_test_scans.mzML")
    only_ms2_path = datafile("only_ms2_mzml.mzML")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path)

    def test_index_building(self):
        try:
            MzMLLoader.prebuild_byte_offset_file(self.path)
            parser = MzMLLoader._parser_cls(self.path)
            assert parser._check_has_byte_offset_file()
            assert isinstance(parser._offset_index,
                              xml_reader.PrebuiltOffsetIndex)
        except OSError:
            pass

    def test_iteration(self):
        reader = self.reader
        i = 0
        bunch = next(reader)
        if bunch.precursor:
            i += 1
        i += len(bunch.products)
        self.assertEqual(i, 3)

        reader.reset()
        reader.make_iterator(grouped=False)

        scan = next(reader)
        scan._load()
        self.assertEqual(scan.index, 0)

        reader.close()

    def test_index(self):
        reader = self.reader
        bunch = next(reader)
        self.assertEqual(bunch.precursor.index, 0)
        for i, scan in enumerate(bunch.products, 1):
            self.assertEqual(scan.index, i)
        reader.close()

    def test_ms_level(self):
        reader = self.reader
        bunch = next(reader)
        self.assertEqual(bunch.precursor.ms_level, 1)
        for i, scan in enumerate(bunch.products, 1):
            self.assertEqual(scan.ms_level, 2)
        reader.close()

    def test_polarity(self):
        reader = self.reader
        bunch = next(reader)
        self.assertEqual(bunch.precursor.polarity, 1)
        reader.close()

    def test_activation(self):
        reader = self.reader
        bunch = next(reader)
        self.assertEqual(bunch.precursor.activation, None)
        for product in bunch.products:
            self.assertNotEqual(product.activation, None)
            self.assertEqual(product.activation.method, "beam-type collision-induced dissociation")
        reader.close()

    def test_precursor(self):
        reader = self.reader
        bunch = next(reader)
        self.assertEqual(bunch.precursor.precursor_information, None)
        for product in bunch.products:
            self.assertNotEqual(product.precursor_information, None)
            self.assertEqual(product.precursor_information.precursor, bunch.precursor)
        reader.close()

    def test_pick_peaks(self):
        reader = self.reader
        bunch = next(reader)
        scan = bunch.precursor.pick_peaks()
        self.assertEqual(len(scan.peak_set), 2108)
        reader.close()

    def test_pack(self):
        reader = self.reader
        bunch = next(reader)
        bunch.precursor.pick_peaks()
        self.assertEqual(bunch.precursor.pack().title, bunch.precursor.title)
        reader.close()

    def test_get_scan_by_id(self):
        reader = self.reader
        precursor = reader.get_scan_by_id(scan_ids[0])
        self.assertEqual(precursor.id, scan_ids[0])
        self.assertEqual(precursor.index, 0)

        product = reader.get_scan_by_id(scan_ids[2])
        self.assertEqual(product.id, scan_ids[2])
        self.assertEqual(product.index, 2)

        self.assertEqual(product.precursor_information.precursor_scan_id, scan_ids[0])
        self.assertIs(precursor, reader.get_scan_by_id(scan_ids[0]))
        reader.close()

    def test_get_scan_by_index(self):
        reader = self.reader

        precursor = reader.get_scan_by_index(0)
        self.assertEqual(precursor.index, 0)
        self.assertEqual(precursor.id, scan_ids[0])
        reader.close()

    def test_get_scan_by_time(self):
        reader = self.reader
        precursor = reader.get_scan_by_time(22.12829)
        self.assertEqual(precursor.id, scan_ids[0])

        product = reader.get_scan_by_time(22.132753)
        self.assertEqual(product.index, 1)
        reader.close()

    def test_instrument_configuration(self):
        reader = self.reader
        bunch = next(reader)
        precursor = bunch.precursor
        config = precursor.instrument_configuration
        self.assertEqual(config.id, "IC1")
        assert "orbitrap" in config.analyzers

    def test_file_description(self):
        file_info = self.reader.file_description()
        assert "MS1 spectrum" in file_info.contents
        assert "MSn spectrum" in file_info.contents
        source_file = file_info.source_files[0]
        assert source_file.name == "three_test_scans.mzML"
        assert "location" not in source_file.parameters

    def test_acquisition_information(self):
        reader = self.reader
        bunch = next(reader)
        precursor = bunch.precursor
        acquisition = precursor.acquisition_information
        self.assertTrue(
            abs(acquisition[0].start_time - precursor.scan_time) < 1e-3)
        self.assertEqual(len(acquisition), 1)
        window = acquisition[0].total_scan_window()
        if window:
            self.assertIn(precursor.arrays.mz[len(precursor.arrays.mz) // 2], window)

    def test_annotations(self):
        reader = self.reader
        bunch = next(reader)
        precursor = bunch.precursor
        assert len(precursor.annotations) > 0

    def test_iteration_mode_detection(self):
        reader = infer_type.MSFileLoader(self.only_ms2_path)
        assert reader.iteration_mode == 'single'

    def test_source_file_parsing(self):
        reader = self.reader
        finfo = reader.file_description()
        sf = finfo.source_files[0]
        assert sf.name == 'three_test_scans.mzML'
        assert isinstance(sf.path, str)

        reader = infer_type.MSFileLoader(self.only_ms2_path)
        finfo = reader.file_description()
        sf = finfo.source_files[0]
        assert sf.name == 'analysis.baf'
        assert isinstance(sf.path, str)

    def test_data_processing_parsing(self):
        reader = infer_type.MSFileLoader(self.only_ms2_path)
        assert len(reader.data_processing()[0]) == 3


if __name__ == '__main__':
    unittest.main()
