.. automodule:: ms_deisotope.processor

    .. autoclass:: ScanProcessor
        :members: pick_precursor_scan_peaks, pick_product_scan_peaks, process_scan_group,
                  process, next, start_from_scan, deconvolute_precursor_scan, deconvolute_product_scan

    .. autofunction:: process


Supporting Types
================

When an ion has been selected with a known m/z and charge, it may be wrapped in an
instance of :class:`PriorityTarget`.

    .. autoclass:: PriorityTarget
        :members: charge_range_hint
