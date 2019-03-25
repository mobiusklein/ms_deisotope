Envelope Search
---------------

The :mod:`ms_deisotope.deconvolution` module contains algorithms for deisotoping
centroided mass spectra, computing the monoisotopic neutral mass and charge of
each fitted envelope. There are two types of searches that can be used, envelope
search can be targeted, using a list of compositions, or exhaustive using an average
monomer isotopic model or :title-reference:`averagine` [Senko]_.

These strategies are implemented as large object hierarchies, but the high level
function :func:`~.deconvolute_peaks` encapsulates the process for most use cases.

.. code:: python

    import ms_deisotope

    deconvoluted_peaks, targeted = ms_deisotope.deconvolute_peaks(peaks, averagine=ms_deisotope.peptide,
                                                                  scorer=ms_deisotope.MSDeconVFitter(10.))

.. contents:: Table of Contents
    :local:

.. automodule:: ms_deisotope.deconvolution

High Level API
==============

    .. autofunction:: deconvolute_peaks


Algorithms for Complex Spectra
==============================

The following algorithms are appropriate for deconvoluting complex, partially
overlapping isotopic envelopes. The default algorithm used by all functions is
:class:`.AveraginePeakDependenceGraphDeconvoluter`.

    .. autoclass:: AveraginePeakDependenceGraphDeconvoluter
        :members: deconvolute, charge_state_determination, targeted_deconvolution,
                  populate_graph

    .. autoclass:: CompositionListPeakDependenceGraphDeconvoluter
        :members: deconvolute, deconvolute_composition, populate_graph

    If your data do not conform to a single averagine, :class:`MultipleAveraginePeakDependenceGraphDeconvoluter`
    can take a :class:`list` of :class:`~.Averagine` objects, selecting the best averagine
    for each experimental isotopic pattern.

Algorithms for Simple Spectra
=============================

These algorithms do not take into account the optimal peak assignment across fits
and should be used only for spot-checks or for simple spectra where best fit resolution
across multiple fits is not required.

    .. autoclass:: AveragineDeconvoluter
        :members: deconvolute, charge_state_determination, targeted_deconvolution

    .. autoclass:: CompositionListDeconvoluter
        :members: deconvolute, deconvolute_composition

Base Classes
============

    .. autoclass:: DeconvoluterBase
        :members:

    .. autoclass:: AveragineDeconvoluterBase
        :members:

    .. autoclass:: CompositionListDeconvoluterBase
        :members:

    .. autoclass:: ExhaustivePeakSearchDeconvoluterBase
        :members:

    .. autoclass:: PeakDependenceGraphDeconvoluterBase
        :members:


Accounting For Lower Quality Data
=================================

.. automodule:: ms_deisotope.deconvolution.peak_retention_strategy

    .. autoclass:: PeakRetentionStrategyBase
        :members:

        .. automethod:: __call__


    .. autoclass:: TopNRetentionStrategy
        :members:

    :data:`ms_deisotope.deconvolution.peak_retention_strategy.simple_peak_retention` is an instance of
    :class:`TopNRetentionStrategy` with `n_peaks=50, base_peak_coefficient=0.05, max_mass=850.0`, the default
    values.


Utilities and Implementation Details
====================================

.. automodule:: ms_deisotope.utils

    .. autoclass:: DeconvolutionProcessResult
        :members:

    .. autoclass:: TargetedDeconvolutionResultBase
        :members:


.. automodule:: ms_deisotope.peak_dependency_network

    .. autoclass:: NetworkedTargetedDeconvolutionResult
        :members:

    .. autoclass:: PeakDependenceGraph


.. [Senko]
    Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
    for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
    Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8


