Isotopic Pattern Generation
---------------------------

The :mod:`ms_deisotope.averagine` module contains several pre-created isotopic models whose
name is derived from the original peptide :title-reference:`averagine` [Senko]_, as well
as isotopic models for other common molecular classes. Isotopic patterns are generated from
a chemical composition using :mod:`brainpy`, which is an implementation of the
:title-reference:`Baffling Recursive Algorithm for Isotopic distributioN calculations` [Dittwald2014]_.

.. code:: python

    import ms_deisotope

    # Senko peptide averagine
    isotopic_pattern = ms_deisotope.peptide.isotopic_cluster(966.12, 2)
    for peak in isotopic_pattern:
        print(peak.mz, peak.intensity)


.. automodule:: ms_deisotope.averagine

    .. autoclass:: Averagine
        :members: scale, isotopic_cluster

    .. autoclass:: TheoreticalIsotopicPattern
        :members: scale, shift, truncate_after, ignore_below

    .. autoclass:: AveragineCache
        :members: isotopic_cluster


.. [Senko]
    Senko, M. W., Beu, S. C., & McLafferty, F. W. (1995). Determination of monoisotopic masses and ion populations
    for large biomolecules from resolved isotopic distributions. Journal of the American Society for Mass
    Spectrometry, 6(4), 229–233. http://doi.org/10.1016/1044-0305(95)00017-8

.. [Dittwald2014]
    Dittwald, P., & Valkenborg, D. (2014). BRAIN 2.0: time and memory complexity improvements in the algorithm for calculating the isotope distribution. Journal of the American Society for Mass Spectrometry, 25(4), 588–94. https://doi.org/10.1007/s13361-013-0796-5