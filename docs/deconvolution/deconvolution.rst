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

Parameter Recomendations and Commentary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Because deconvolution can be a fussy problem to solve when dealing with disparate instrument types and
signal qualities, here are some notes when tuning parameters.


Missing High Mass Ions? Understanding how and when to use `truncate_after`
**************************************************************************

.. note:: Short Answer: Use 0.95 for MS1 and 0.8 for MSn spectra, or experiment with ``incremental_truncation``.

Theoretical isotopic patterns can simulate any number of isotopic peaks, though outside of native MS of large
molecules, most of these isotopic peaks are trivially small and undetectable. :mod:`ms_deisotope` uses a truncation
strategy to stop including isotopic peaks after ``truncate_after`` * 100% of the signal has been observed in an
isotopic pattern. By default ``truncate_after`` is 0.95, which I've found works well for MS1 spectra on all instruments.
This is trivial for most biomolecules below 800 Da, but as molecules grow larger, the theoretical pattern starts to
exceed the instrument's ability to detect all peaks within 95% of the isotopic pattern, especially on highly processed
spectra like FTMS from Orbitraps. For MSn spectra, ions often appear on both sides of this line, making it especially
important to account for. For these spectra, I have found ``truncate_after`` at 0.8 works well.

There is another option, which is to use the experimental ``incremental_truncation`` strategy which takes each isotopic
pattern after ``truncate_after`` and considers variants incrementally dropping trailing peaks until ``incremental_truncation``
* 100% of the signal from the starting pattern remains. Starting with ``truncate_after`` set to 0.999 and ``incremental_truncation``
set to 0.8 will allow you to match both ends of the spectrum, but at the cost of doing a bit more work per peak.


`quick_charge` To Speed Things Up
*********************************
Use ``quick_charge`` = :const:`True` when missing peaks around the monoisotopic peak are not an issue to drastically cut
down on the number of theoretical pattern fits that need to be performed per peak.


Searching Left and Right
************************
When using one of the graph-based algorithms for complex spectra (the default), these parameters are superfluous. When using
one of the much simpler algorithms suitable for deconvolution of limited complexity data or targeted deconvolution only, these
parameters benefit from being set to a value between 1 and 3, depending upon the distance form the base peak of an isotopic pattern
to the monoisotopic peak.


When To Use `priority_list`
***************************
The ``priority_list`` parameter is really only useful when you need to tie a specific deconvolution result to a specific experimental
peak, as in the case of extracting the deconvoluting a precursor ion while also deconvolving the entire spectrum. Unless you have these
types of scenarios, this argument is not generally necessary. It does get used heavily by :class:`ms_deisotope.processor.ScanProcessor`.


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

.. automodule:: ms_deisotope.deconvolution.utils

    .. autofunction:: prepare_peaklist

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


