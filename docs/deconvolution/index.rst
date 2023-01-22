Deconvolution
-------------

:mod:`ms_deisotope`'s primary function is to deisotope and charge state deconvolute
mass spectra. It contains several approaches to this problem class. This can be broken
into two cooperative procedures, *envelope search* and *envelope scoring*. Envelope
search can be targeted, using a list of compositions, or exhaustive using an average
monomer isotopic model or :title-reference:`averagine` [Senko]_.

.. toctree::
    :caption: Topics

    Envelope Search <deconvolution>
    Isotopic Pattern Generation <averagine>
    Envelope Scoring <envelope_scoring>
    Deconvolution Pipeline <pipeline>
    peaks