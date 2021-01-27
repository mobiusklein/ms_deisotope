.. image:: https://raw.githubusercontent.com/mobiusklein/ms_deisotope/master/docs/_static/logo.png
    :target: https://mobiusklein.github.io/ms_deisotope

`Documentation <https://mobiusklein.github.io/ms_deisotope>`_ | |PYPIBADGE| | |GHAB|


A Library for Deisotoping and Charge State Deconvolution For Mass Spectrometry
------------------------------------------------------------------------------

This library combines `brainpy` and `ms_peak_picker` to build a toolkit for
MS and MS/MS data. The goal of these libraries is to provide pieces of the puzzle
for evaluating MS data modularly. The goal of this library is to combine the modules
to streamline processing raw data.

Deconvolution
=============

The general-purpose averagine-based deconvolution procedure can be called by using the high level
API function `deconvolute_peaks`, which takes a sequence of peaks, an averagine model, and a isotopic
goodness-of-fit scorer:

.. code:: python

    import ms_deisotope

    deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(peaks, averagine=ms_deisotope.peptide,
                                                           scorer=ms_deisotope.MSDeconVFitter(10.))

The result is a deisotoped and charge state deconvoluted peak list where each peak's neutral mass is known
and the fitted charge state is recorded along with the isotopic peaks that gave rise to the fit.

Refer to the `Documentation <https://mobiusklein.github.io/ms_deisotope>`_ for a deeper description of isotopic pattern fitting.

Averagine
=========

An "Averagine" model is used to describe the composition of an "average amino acid",
which can then be used to approximate the composition and isotopic abundance of a
combination of specific amino acids. Given that often the only solution available is
to guess at the composition of a particular *m/z* because there are too many possible
elemental compositions, this is the only tractable solution.

This library supports arbitrary Averagine formulae, but the Senko Averagine is provided
by default: `{"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}`

.. code:: python

    from ms_deisotope import Averagine
    from ms_deisotope import plot

    peptide_averagine = Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})

    plot.draw_peaklist(peptide_averagine.isotopic_cluster(1266.321, charge=1))


`ms_deisotope` includes several pre-defined averagines (or "averagoses" as may be more appropriate):
    1. Senko's peptide - `ms_deisotope.peptide`
    2. Native *N*- and *O*-glycan - `ms_deisotope.glycan`
    3. Permethylated glycan - `ms_deisotope.permethylated_glycan`
    4. Glycopeptide - `ms_deisotope.glycopeptide`
    5. Sulfated Glycosaminoglycan - `ms_deisotope.heparan_sulfate`
    6. Unsulfated Glycosaminoglycan - `ms_deisotope.heparin`


Please see the `Documentation <https://mobiusklein.github.io/ms_deisotope>`_ for more information on mass spectrum data file reading/writing, peak sets, and lower-level signal processing tools.


Installing
----------

``ms_deisotope`` uses PEP 517 and 518 build system definition and isolation to ensure all of its
compile-time dependencies are installed prior to building. Normal installation should work with `pip`,
and pre-built wheels are available for Windows.

.. code:: bash

    $ pip install ms_deisotope

C Extensions
============

``ms_deisotope`` and several of its dependencies use C extensions to make iterative operations *much*
faster. If you plan to use this library on a large amount of data, I highly recommend you ensure they
are installed:

.. code:: python

    >>> import ms_deisotope
    >>> ms_deisotope.DeconvolutedPeak
    <type 'ms_deisotope._c.peak_set.DeconvolutedPeak'>

Building C extensions from source requires a version of Cython >= 0.27.0

Compiling C extensions requires that ``numpy``, ``brain-isotopic-distribution``, and ``ms_peak_picker``
be compiled and installed prior to building ``ms_deisotope``:

.. code:: bash

    pip install numpy
    pip install -v brain-isotopic-distribution ms_peak_picker
    pip install -v ms_deisotope

If these libraries are not installed, ``ms_deisotope`` will fall back to using pure Python implementations,
which are much slower.


.. |PYPIBADGE| image:: https://badge.fury.io/py/ms-deisotope.svg
    :target: https://badge.fury.io/py/ms-deisotope
.. |GHAB| image:: https://github.com/mobiusklein/ms_deisotope/workflows/tests/badge.svg