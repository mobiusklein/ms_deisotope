import unittest
import os
import tempfile
import gzip

import ms_deisotope
from ms_deisotope.data_source import MzMLLoader
from ms_deisotope.data_source import _compression
from ms_deisotope.test.common import datafile

from ms_deisotope.output.mzml import MzMLSerializer, ProcessedMzMLDeserializer


class TestMzMLSerializer(unittest.TestCase):
    source_data_path = datafile("three_test_scans.mzML")

    def test_writer(self):
        source_reader = MzMLLoader(self.source_data_path)
        fd, name = tempfile.mkstemp()
        with open(name, 'wb') as fh:
            writer = MzMLSerializer(fh, n_spectra=len(source_reader.index), deconvoluted=True)
            description = source_reader.file_description()
            writer.add_file_information(description)
            writer.add_file_contents("profile spectrum")
            writer.add_file_contents("centroid spectrum")
            writer.remove_file_contents("profile spectrum")

            instrument_configs = source_reader.instrument_configuration()
            for config in instrument_configs:
                writer.add_instrument_configuration(config)

            software_list = source_reader.software_list()
            for software in software_list:
                writer.add_software(software)

            data_processing_list = source_reader.data_processing()
            for dp in data_processing_list:
                writer.add_data_processing(dp)

            processing = writer.build_processing_method()
            writer.add_data_processing(processing)
            bunch = next(source_reader)
            bunch.precursor.pick_peaks()
            bunch.precursor.deconvolute()
            for product in bunch.products:
                product.pick_peaks()
                product.deconvolute()
            writer.save(bunch)
            writer.complete()
            fh.flush()
            writer.format()
        source_reader.reset()
        processed_reader = ProcessedMzMLDeserializer(_compression.get_opener(writer.handle.name))

        for a, b in zip(source_reader.instrument_configuration(), processed_reader.instrument_configuration()):
            assert a.analyzers == b.analyzers
        for a, b in zip(source_reader, processed_reader):
            assert a.precursor.id == b.precursor.id
            assert (a.precursor.acquisition_information == b.precursor.acquisition_information)
            for an, bn in zip(a.products, b.products):
                assert an.id == bn.id
                assert abs(an.precursor_information.neutral_mass - bn.precursor_information.neutral_mass) < 1e-6
        processed_reader.reset()
        description = processed_reader.file_description()
        assert "profile spectrum" not in description.contents
        assert "centroid spectrum" in description.contents
        sf = description.source_files[0]
        assert 'location' not in sf.parameters
        assert sf.parameters['SHA-1'] == 'a2a091b82f27676da87a6c7d17cc90d2d90b8fbf'
        assert processed_reader.sample_run.name == "sample_1"
        index = processed_reader.extended_index
        pinfo = index.find_msms_by_precursor_mass(ms_deisotope.neutral_mass(562.7397, 2))
        assert len(pinfo) > 0

        software_list = processed_reader.software_list()
        for sw in software_list:
            if sw.name == 'ms_deisotope' and sw.id == 'ms_deisotope_1':
                break
        else:
            raise AssertionError("Could not find \"ms_deisotope\" entry in software list %r" % (software_list, ))

        processed_reader.close()
        try:
            os.remove(name)
            os.remove(processed_reader._index_file_name)
        except OSError as _err:
            pass

if __name__ == '__main__':
    unittest.main()
