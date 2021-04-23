.. automodule:: ms_deisotope.output.mzmlb

    .. autoclass:: ms_deisotope.output.mzml.MzMLbSerializer
        :members: save, save_scan_bunch, save_scan, add_sample, add_instrument_configuration, add_source_file,
                  add_file_information, add_software, add_file_contents, add_processing_parameter,
                  add_data_processing, add_software, save_chromatogram, close


    .. autoclass:: ms_deisotope.output.mzml.ProcessedMzMLbLoader
        :members: get_index_information_by_scan_id, has_index_file,
                  msms_for, precursor_information
