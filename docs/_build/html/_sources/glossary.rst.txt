Glossary
--------

Library Concept Glossary
========================

.. glossary::

    single
        Single iteration mode iterators will produce a single :class:`~ms_deisotope.data_source.Scan`-like object
        per iteration, regardless of :attr:`ms_deisotope.data_source.Scan.ms_level`. This is the default mode for
        formats like :class:`~ms_deisotope.data_source.mgf.MGFLoader`.

    grouped
        Grouped iteration mode iterators will produce a :class:`~ms_deisotope.data_source.ScanBunch` object per
        iteration. A :class:`~ms_deisotope.data_source.ScanBunch` usually contains a precursor MS1 spectrum (but not
        always) plus a list of zero or more MSn spectra which are derived from it or which have no locate-able
        precursor.


File Format Glossary
====================

.. glossary::

    mzML : mzml
    mzml : mzml
        A standard rich XML-format for raw mass spectrometry data storage, working as both
        a centroid :term:`peak list format` or continuous profiles. See `http://www.psidev.info/ <http://www.psidev.info/index.php?q=node/257>`_
        for additional information. This is the preferred open format for data exchange.

    MGF : mgf
    mgf : mgf
        The "Mascot Generic File" format, a plain text multiple :term:`peak list format`.

    peak list format
        A file or data serialization format the encodes a list of mass spectrum peaks with
        or without additional metadata. A format may either encode multiple spectra or a single
        spectrum, depending upon design.



