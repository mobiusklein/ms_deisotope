File Metadata
-------------

Mass spectrometry data files contain global metadata describing the type of data they contain, from
the types of spectra or MS level to the component of the machine used to acquire the data. These
structures can describe these properties.

Most of the time these are accessed through methods on a :class:`ScanDataSource` or
attribute on a :class:`Scan` object.

.. toctree::
    :caption: Topics

    file_metadata
    instrument_configuration
    activation
    cv