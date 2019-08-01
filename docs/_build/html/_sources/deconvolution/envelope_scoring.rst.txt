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

    .. autoclass:: GTestFitter

    .. autoclass:: ScaledGTestFitter

    .. autoclass:: LeastSquaresFitter

Maximizing Fitters
~~~~~~~~~~~~~~~~~~

    .. autoclass:: MSDeconVFitter

    .. autoclass:: PenalizedMSDeconVFitter

Support Structures
==================

    .. autoclass:: IsotopicFitRecord

    .. autoclass:: FitSelectorBase

    .. autoclass:: MaximizeFitSelector

    .. autoclass:: MinimizeFitSelector