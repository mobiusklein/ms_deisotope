.. ms_deisotope documentation master file, created by
   sphinx-quickstart on Mon Feb 05 17:40:57 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ms_deisotope's documentation!
========================================

:mod:`ms_deisotope` contains a variety of implementations for deisotoping and charge
state deconvolution of mass spectra. It uses :mod:`ms_peak_picker` to perform the low
level signal processing needed to translate profile spectra into centroided *m/z* peak lists
and :mod:`brainpy` to generate theoretical isotopic patterns that it fits to these centroids
to produce deconvoluted *neutral mass* peak lists.

:mod:`ms_deisotope` also includes a data access layer for reading mass spectra from several
common formats, including instrument information and and file metadata.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Data Access <data_source/index>
   Deconvolution <deconvolution/index>



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
