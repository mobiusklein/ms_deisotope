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
