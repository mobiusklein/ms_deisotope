Common MS File Model
--------------------

    :mod:`ms_deisotope.data_source` uses a set of common interfaces for reading
    mass spectrometry data files so that code written for one format should work
    for all formats which implement the same interfaces.

    :mod:`ms_deisotope.data_source.metadata` defines a set of data structures and
    a collection of controlled vocabulary terms describing mass spectrometers and
    mass spectrometry data files.


.. automodule:: ms_deisotope.data_source.common

    Abstract Base Classes
    =======================

    :mod:`ms_deisotope` supports reading from many different file formats. While
    the file format is abstracted away as much as possible, the modes of access are
    built into the type hierarchy.

    All of the currently implemented formats implement both :class:`ScanIterator` and
    :class:`RandomAccessScanSource`.

    .. autoclass:: ScanDataSource
        :members:


    .. autoclass:: ScanIterator
        :members:


    .. autoclass:: RandomAccessScanSource
        :members:


File-Level Metadata
===================
    
File Description
~~~~~~~~~~~~~~~~

        A mass spectrometry data file contains heterogenous types of information
        derived from one or more provenance sources. Some formats, like `mzML
        <https://doi.org/10.1074/mcp.R110.000133>`_, track this
        information. This information can be queried to make the object select
        behavior appropriate to its contents and to provide a clear chai of
        evidence back to the original raw data.


        .. automodule:: ms_deisotope.data_source.metadata.file_information

            .. autoclass:: FileInformation
                :members:

            .. autoclass:: SourceFile
                :members:


Instrument Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

        A mass spectrometer may be configured to operate in multiple modes during
        a single data acquisition run. A program may need to know which configurations
        were available, and formats like :title-reference:`mzML` are able record this
        information. To retrieve the configuration used for a particular scan, see
        :attr:`.Scan.instrument_configuration`

        .. automodule:: ms_deisotope.data_source.metadata.instrument_components

            .. autoclass:: Component
                :members:
                
                .. automethod:: is_a

            .. autoclass:: ComponentGroup
                :members:

            .. autoclass:: InstrumentInformation
                :members:
