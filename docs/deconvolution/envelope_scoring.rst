Envelope Scoring
----------------

The :mod:`ms_deisotope.scoring` module contains classes for evaluating the goodness-of-fit
of isotopic pattern matches. It is used by deconvoluters defined in :mod:`ms_deisotope.deconvolution`
to decide which pattern fit is best. Each instance of :class:`IsotopicFitterBase` takes a score
threshold, a :class:`float`, and additional configuration arguments. The score threshold is used
to filter out noise matches.

When the isotopic patterns of the target data are of high quality and intensity can be used
to effectively filter out noise peaks, :class:`~.PenalizedMSDeconVFitter` often works best. If
the isotopic patterns are not high quality, but intensity can be used to discriminate noise,
:class:`~.MSDeconVFitter` is more forgiving. If intensity is not important, one of either
:class:`~.ScaledGTestFitter` or :class:`~.LeastSquaresFitter` may work.

.. automodule:: ms_deisotope.scoring

Scoring Functions
=================

    .. autoclass:: IsotopicFitterBase
        :members: evaluate, reject, is_maximizing, __call__

Minimizing Fitters
~~~~~~~~~~~~~~~~~~

Minimizing envelope scoring methods aim to minimize some "fit error" criterion, and discard solutions which don't fit well.
They tend to be quite sensitive and can be independent of the magnitude of the signal. This also means they do not handle
detector noise or interference well. These methods can work well when targeting a list of compositions instead of exhaustively
deconvolving an entire spectrum.

    .. autoclass:: GTestFitter

    .. autoclass:: ScaledGTestFitter

    .. autoclass:: LeastSquaresFitter

Maximizing Fitters
~~~~~~~~~~~~~~~~~~

Maximizing envelope scoring methods aim to maximize some "goodness-of-fit" criterion, and discards solutions that don't score
highly enough. While not universally true, many of the maximizing scoring functions here are a function of the magnitude of
the signal, which means that the threshold selected is signal magnitude dependent. This has the advantage of making the score
threshold also a detector noise filter, but the threshold would now be instrument type-dependent.

A maximizing fitter with a well-chosen threshold is good at exhaustively deconvoluting a spectrum because it can more easily
eliminate bad fits, but as isotopic pattern becomes more homogenous around the base peak of the isotopic pattern, they can more
easily make mistakes identifying the monoisotopic peak. This is where :class:`MSDeconVFitter` is less accurate than :class:`PenalizedMSDeconVFitter`,
albeit neither is perfect.

    .. autoclass:: MSDeconVFitter

    .. autoclass:: PenalizedMSDeconVFitter

Other Fitters
~~~~~~~~~~~~~

    .. autoclass:: FunctionScorer

Support Structures
==================

    .. autoclass:: IsotopicFitRecord

    .. autoclass:: FitSelectorBase

    .. autoclass:: MaximizeFitSelector

    .. autoclass:: MinimizeFitSelector