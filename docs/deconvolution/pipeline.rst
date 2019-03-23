Deconvolution Pipeline
----------------------

The deconvolution process from start to finish can be pipelined from start to finish
using the :class:`~.ScanProcessor` class. This includes precursor recalculation and
coisolation detection.

.. code:: python

    from ms_deisotope import ScanProcessor, glycopeptide, peptide
    from ms_deisotope.scoring import PenalizedMSDeconVFitter, MSDeconVFitter
    from ms_deisotope.test.common import datafile

    proc = processor.ScanProcessor(datafile("20150710_3um_AGP_001_29_30.mzML.gz"), ms1_deconvolution_args={
        "averagine": glycopeptide,
        "scorer": PenalizedMSDeconVFitter(20., 2.),
        "truncate_after": 0.95
    }, msn_deconvolution_args={
        "averagine": peptide,
        "scorer": MSDeconVFitter(10.),
        "truncate_after": 0.8
    })

    bunch = next(proc)
    print(bunch)
    print(bunch.precursor.deconvoluted_peak_set)


.. automodule:: ms_deisotope.processor

    .. autoclass:: ScanProcessor
        :members: pick_precursor_scan_peaks, pick_product_scan_peaks, process_scan_group,
                  process, next, start_from_scan, deconvolute_precursor_scan, deconvolute_product_scan

Supporting Types
================

When an ion has been selected with a known m/z and charge, it may be wrapped in an
instance of :class:`PriorityTarget`.

    .. autoclass:: PriorityTarget
        :members: charge_range_hint
