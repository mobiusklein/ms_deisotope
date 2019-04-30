Visualizing Mass Spectra
------------------------

The :mod:`ms_deisotope.plot` module contains several plotting functions for visualizing
mass spectra using :mod:`matplotlib`.

.. automodule:: ms_deisotope.plot


    .. plot::
        :include-source:

        import ms_deisotope
        from ms_deisotope import plot
        from ms_deisotope.test.common import datafile

        reader = ms_deisotope.MSFileLoader(datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
        bunch = next(reader)

        # create a profile spectrum
        for peak in bunch.precursor.pick_peaks():
            peak.full_width_at_half_max = 0.02
        scan = bunch.precursor.reprofile(dx=0.001)
        ax = plot.draw_raw(scan.arrays, color='black', lw=0.5)
        ax.set_xlim(1160, 1165)
        ax.figure.set_figwidth(12)
        ax.set_title("Raw Profile Plot", fontsize=16)

        scan.pick_peaks()
        ax = plot.draw_peaklist(scan.peak_set, color='black')
        ax.set_xlim(1160, 1165)
        ax.figure.set_figwidth(12)
        ax.set_title("Centroid Peak Plot", fontsize=16)


Basic Spectrum Drawing
======================

    .. autofunction:: draw_raw

    .. autofunction:: draw_peaklist


Annotating Peaks and Envelopes
==============================

    .. autofunction:: annotate_isotopic_peaks

    .. autofunction:: label_peaks

    .. autofunction:: annotate_scan_single

    .. autofunction:: annotate_scan
