Data Access
-----------

:mod:`ms_deisotope` can read from mzML, mzXML and MGF files directly, using the :mod:`pyteomics` library.
On Windows, it can also use :mod:`comtypes` to access Thermo Fisher's `MSFileReader
<https://thermo.flexnetoperations.com/control/thmo/index>`_ to read RAW files and Agilent's
MassSpecDataReader to read .d directories. Whenever possible, the library provides a
common interface to all supported formats.

To read from a supported file format, do one of the following:

.. code:: python

    import ms_deisotope

    # automatically determine the appropriate reader type
    # by inspecting the file
    reader = ms_deisotope.MSFileReader("path/to/file.mzML")

    # explicitly select the reader type to use if the file
    # is compressed or ambiguous
    reader = ms_deisotope.MzMLLoader("path/to/file.mzML")



.. toctree::
    :caption: Topics

    common_scan
    common_reader
    mzml
    mzxml
    formats
