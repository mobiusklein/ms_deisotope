"""
Writing mzMLb
-------------

Using the :mod:`psims` library, :mod:`ms_deisotope.output.mzmlb` can write an mzMLb
file with all associated metadata. This uses the same machinery as :mod:`ms_deisotope.output.mzml`
to accomplish this.

This module also contains a specialized version of :class:`~.MzMLbLoader`,
:class:`~.ProcessedMzMLbLoader`, which can directly reconstruct each
deconvoluted peak list and provides fast access to an extended index of
metadata that :class:`~.MzMLbSerializer` writes to an external file.

All behavior available for the :class:`~.ProcessedMzMLLoader` are available for
:class:`~.ProcessedMzMLbLoader`.


.. code:: python

    import ms_deisotope
    from ms_deisotope.test.common import datafile
    from ms_deisotope.output.mzmlb import MzMLbSerializer

    reader = ms_deisotope.MSFileLoader(datafile("small.mzML"))
    with MzMLbSerializer('small.deconvoluted.mzMLb', n_spectra=len(reader)) as writer:
        writer.copy_metadata_from(reader)
        for bunch in reader:
            bunch.precursor.pick_peaks()
            bunch.precursor.deconvolute()
            for product in bunch.products:
                product.pick_peaks()
                product.deconvolute()
            writer.save(bunch)

"""

import numpy as np

from psims.mzml.writer import COMPRESSION_ZLIB

try:
    from psims.mzmlb.writer import MzMLbWriter as _MzMLbWriter
    from psims.mzmlb.writer import DEFAULT_COMPRESSOR
    is_available = True
except ImportError:
    DEFAULT_COMPRESSOR = COMPRESSION_ZLIB
    _MzMLbWriter = None
    is_available = False

from ms_deisotope.data_source.mzmlb import MzMLbLoader

from .common import LCMSMSQueryInterfaceMixin, ScanDeserializerBase, SampleRun
from .mzml import _IonMobility3DPeakPacker, MzMLSerializer as _MzMLSerializer, PeakSetDeserializingMixin, _IonMobility3DSerializerBase


class MzMLbSerializer(_MzMLSerializer):
    """
    Write :mod:`ms_deisotope` data structures to a file in mzMLb format.

    This type inherits most of its functionality from :class:`~ms_deisotope.output.mzml.MzMLSerializer`.
    """

    file_extensions = {
        "mzmlb",
    }

    _format_conversion_term = "Conversion to mzMLb"

    default_compression = DEFAULT_COMPRESSOR

    def _make_writer(self, handle):
        compression = self.compression
        compression_opts = None
        if self.compression == COMPRESSION_ZLIB:
            compression = 'gzip'
            compression_opts = 4
        self.compression = 'none'
        return _MzMLbWriter(
            handle,
            close=self._should_close,
            h5_compression=compression,
            h5_compression_options=compression_opts
        )

    @classmethod
    def get_scan_packer_type(cls):
        return super().get_scan_packer_type()



class IonMobilityAware3DMzMLbSerializer(_IonMobility3DSerializerBase, MzMLbSerializer):
    default_data_encoding = MzMLbSerializer.default_data_encoding.copy()
    default_data_encoding.update({
        "feature id array": np.int32,
    })

    @classmethod
    def get_scan_packer_type(cls):
        return _IonMobility3DPeakPacker


class ProcessedMzMLbLoader(PeakSetDeserializingMixin, MzMLbLoader, ScanDeserializerBase, LCMSMSQueryInterfaceMixin):
    """
    Extends :class:`.MzMLbLoader` to support deserializing preprocessed data and to provide indexing information.

    Attributes
    ----------
    extended_index: :class:`~.ExtendedIndex`
        Holds the additional indexing information
        that may have been generated with the data
        file being accessed.
    sample_run: :class:`SampleRun`
    """

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
