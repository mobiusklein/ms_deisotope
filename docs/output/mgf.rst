Writing MGF
-----------

:mod:`ms_deisotope.output.mgf` can write an MGF for centroided or deconvoluted spectra.


.. code:: python

    from ms_deisotope.output.mgf import MGFSerializer

    def to_mgf(reader, outstream, msn_filters=None):
        if not msn_filters:
            msn_filters = []
        reader.make_iterator(grouped=False)
        writer = MGFSerializer(outstream, deconvoluted=False)
        with writer:
            for scan in reader:
                if scan.ms_level == 1:
                    continue
                if msn_filters:
                    scan = scan.transform(msn_filters)
                if scan.peak_set is None:
                    scan.pick_peaks()
                writer.save_scan(scan)


.. automodule:: ms_deisotope.output.mgf

    .. autoclass:: ms_deisotope.output.mgf.MGFSerializer
        :members: save_scan_bunch, add_parameter, add_global_parameter

    .. autoclass:: ms_deisotope.output.mgf.ProcessedMGFDeserializer


