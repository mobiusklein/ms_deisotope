Writing mzML
------------

Using the :mod:`psims` library, :mod:`ms_deisotope.output.mzml` can write an mzML
file with all associated metadata, including deconvoluted peak arrays, chromatograms,
and data transformations. The :class:`~.MzMLSerializer` class handles all facets of
this process.

This module also contains a specialized version of :class:`~.MzMLLoader`,
:class:`~.ProcessedMzMLDeserializer`, which can directly reconstruct each
deconvoluted peak list and provides fast access to an extended index of
metadata that :class:`~.MzMLSerializer` writes to an external file.


.. automodule:: ms_deisotope.output.mzml
    
    .. autoclass:: ms_deisotope.output.mzml.MzMLSerializer
        :members: save_scan_bunch, add_sample, add_instrument_configuration, add_source_file,
                  add_file_information, add_software, add_file_contents, add_processing_parameter,
                  save_chromatogram, complete, format

    .. autoclass:: ms_deisotope.output.mzml.ProcessedMzMLDeserializer
        :members: get_index_information_by_scan_id, has_index_file,
                  msms_for, precursor_information
