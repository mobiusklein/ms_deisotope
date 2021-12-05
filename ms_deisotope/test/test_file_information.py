import unittest

from ms_deisotope.test.common import datafile
from ms_deisotope.data_source import infer_type
from ms_deisotope.data_source.metadata import file_information


class TestFileMetadata(unittest.TestCase):
    path = datafile("three_test_scans.mzML")

    @property
    def reader(self):
        return infer_type.MSFileLoader(self.path)

    def test_file_information(self):
        reader = self.reader
        finfo = reader.file_description()
        assert "MS1 spectrum" in finfo
        assert reader.id_format == "no nativeID format"

    def test_source_file(self):
        id_fmt, fmt = file_information.SourceFile.guess_format(self.path)
        assert fmt == "mzML format"
        assert id_fmt == "no nativeID format"

        sf = file_information.SourceFile.from_path(self.path)
        assert not sf.has_checksum()
        sf.add_checksum('sha1')
        assert sf.has_checksum('sha1') and sf.has_checksum()
        assert sf.parameters["SHA-1"] == sf.checksum("sha1")
        assert sf.validate_checksum()

        other = sf.copy()
        assert sf == other


if __name__ == '__main__':
    unittest.main()
