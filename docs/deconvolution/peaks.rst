Deconvoluted Peak Sets
----------------------

.. module:: ms_deisotope.peak_set

.. autoclass:: DeconvolutedPeak

.. class:: DeconvolutedPeakSet

    A :class:`Sequence` of :class:`DeconvolutedPeak` objects, ordered by neutral mass. This class is
    intended to support fast searches by mass.

    The recommended query method is :meth:`all_peaks_for` when you care about all distinct charge states,
    while :meth:`has_peak` can be used when you only care about the closest peak matching the query mass.

    .. note::

        :meth:`has_peak` and :meth:`between` support searching for m/z values instead of intensity values, but
        this is purely for convenience. The PPM error scale on m/z is different than on neutral mass, and should
        not be assumed to give identical results as the neutral mass query methods.

    If you wish to use the resulting peaks in another tool that assumes all peaks are singly charged, see
    :func:`decharge`.

    .. automethod:: has_peak

    .. automethod:: all_peaks_for

    .. automethod:: between




Transformations
===============

.. autofunction:: decharge

.. autofunction:: merge

.. autofunction:: envelopes_to_peak_set

.. autofunction:: window_peak_set


Peak Components
===============
Deconvoluted peaks may contain extra information beyond point-estimates from the raw spectrum. These data structures
hold some of that information.

.. autoclass:: Envelope

.. autoclass:: _Index


Other Deconvolution Product Types
=================================

.. autoclass:: DeconvolutedPeakSolution

