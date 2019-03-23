Deconvolution Pipeline
----------------------

The deconvolution process from start to finish can be pipelined from start to finish
using the :class:`~.ScanProcessor` class.


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



.. automodule:: ms_deisotope.processor


    .. autoclass:: ScanProcessor
        :members: pick_precursor_scan_peaks, pick_product_scan_peaks, process_scan_group,
                  process, next, start_from_scan, _get_envelopes, _average_ms1, deconvolute_precursor_scan,
                  deconvolute_product_scan


    .. autoclass:: PriorityTarget
        :members: charge_range_hint
