Deconvolution Pipeline
----------------------

The deconvolution process from start to finish can be pipelined from start to finish
using the :class:`~.ScanProcessor` class.

.. automodule:: ms_deisotope.processor
    

    .. autoclass:: ScanProcessor
        :members: pick_precursor_scan_peaks, pick_product_scan_peaks, process_scan_group,
                  process, next, start_from_scan, _get_envelopes, _average_ms1, deconvolute_precursor_scan,
                  deconvolute_product_scan


    .. autoclass:: PriorityTarget
        :members: charge_range_hint
