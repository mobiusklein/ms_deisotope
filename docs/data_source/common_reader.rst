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


    Iteratation Strategies
    ======================

    :class:`ScanIterator` instances may iterate over scans in :term:`single` or :term:`grouped` strategies.
    :term:`single` mode produces a single instance of :class:`~.Scan` on each iteration, while :term:`grouped`
    produces a :class:`~.ScanBunch` containing an MS1 :class:`~.Scan` (may be :const:`None`) and 0 or more
    related MSn :class:`~.Scan` instances which are derived from the MS1 :class:`~.Scan`. The default
    mode for a given :class:`ScanIterator` depends upon both the file format and available metadata.

    You can force the iteration strategy to be :term:`grouped` when calling :meth:`ScanIterator.make_iterator`
    by passing ``grouped=True``, and :term:`single` by passing ``grouped=False``. The same applies to
    :meth:`RandomAccessScanSource.start_from_scan`. When :term:`grouped` mode is requested but cannot be
    fulfilled, :class:`~.ScanBunch` objects are still produced, but the :attr:`~.ScanBunch.precursor` may be
    :const:`None` or :attr:`~.ScanBunch.products` may be empty.

    The iteration mode of a :class:`ScanIterator` is always available through it's :attr:`iteration_mode`
    attribute, which should have the value ``"single"`` or ``"grouped"`` accordingy.