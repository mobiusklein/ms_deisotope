from psims.mzml.writer import COMPRESSION_ZLIB

try:
    from psims.mzmlb.writer import MzMLbWriter as _MzMLbWriter
    is_available = True
except ImportError:
    _MzMLbWriter = None
    is_available = False

from ms_deisotope.data_source.mzmlb import MzMLbLoader

from .common import LCMSMSQueryInterfaceMixin, ScanDeserializerBase, SampleRun
from .mzml import MzMLSerializer as _MzMLSerializer, PeakSetDeserializingMixin


class MzMLbSerializer(_MzMLSerializer):

    file_extensions = {
        "mzmlb",
    }


    def _make_writer(self, handle):
        compression = self.compression
        compression_opts = None
        if self.compression == COMPRESSION_ZLIB:
            compression = 'gzip'
            compression_opts = 4
        self.compression = 'none'
        return _MzMLbWriter(self.handle, h5_compression=compression, h5_compression_options=compression_opts)



class ProcessedMzMLbLoader(PeakSetDeserializingMixin, MzMLbLoader, ScanDeserializerBase, LCMSMSQueryInterfaceMixin):

    file_extensions = {
        "mzmlb",
    }

    def __init__(self, source_file, use_index=True, use_extended_index=True):
        super(ProcessedMzMLbLoader, self).__init__(
            source_file, use_index=use_index, decode_binary=True)
        self.extended_index = None
        self._scan_id_to_rt = dict()
        self._sample_run = None
        self._use_extended_index = use_extended_index
        if self._use_index:
            if self._use_extended_index:
                self.require_extended_index()

    def _dispose(self):
        self._scan_id_to_rt.clear()
        self.extended_index.clear()
        super(ProcessedMzMLbLoader, self)._dispose()

    def __reduce__(self):
        return self.__class__, (self.source_file, self._use_index, self._use_extended_index)

    def _make_sample_run(self):
        samples = self.samples()
        sample = samples[0]
        return SampleRun(name=sample.name, uuid=sample['SampleRun-UUID'], **dict(sample.items()))

    @property
    def sample_run(self):
        if self._sample_run is None:
            self._sample_run = self._make_sample_run()
        return self._sample_run


ProcessedMzMLbDeserializer = ProcessedMzMLbLoader
