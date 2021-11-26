# -*- coding: utf-8 -*-
'''mzMLb is a standard rich HDF5-based format for raw mass spectrometry data storage.
This module provides :class:`MzMLbLoader`, a :class:`~.RandomAccessScanSource`
implementation. It is based upon the mzML XML file format, re-using a subset of the
features. The original design for mzMLb is described in [Bhamber]_.

The parser is based on :mod:`pyteomics.mzmlb`. It requires :mod:`h5py` to be installed
for reading, and :mod:`hdf5plugin` to use the faster, non-zlib-based compressors.

References
==========
.. [Bhamber] Bhamber, R. S., Jankevics, A., Deutsch, E. W., Jones, A. R., & Dowsey, A. W. (2021).
    MzMLb: A Future-Proof Raw Mass Spectrometry Data Format Based on Standards-Compliant
    mzML and Optimized for Speed and Storage Requirements. Journal of Proteome Research,
    20(1), 172–183. https://doi.org/10.1021/acs.jproteome.0c00192
'''

import logging
logging.getLogger("hdf5plugin").addHandler(logging.NullHandler())

try:
    from pyteomics import mzmlb
    _BaseParser = mzmlb.MzMLb
except ImportError:
    mzmlb = None
    _BaseParser = object

from .mzml import MzMLLoader as _MzMLLoader
from ._compression import DefinitelyFastRandomAccess


class _MzMLbParser(_BaseParser):

    def _handle_param(self, element, **kwargs):
        try:
            element.attrib["value"]
        except KeyError:
            element.attrib["value"] = ""
        return super(_MzMLbParser, self)._handle_param(element, **kwargs)


class MzMLbLoader(_MzMLLoader):
    """Reads scans from PSI-HUPO mzMLb HDF5 files. Provides both iterative and
    random access.

    Attributes
    ----------
    source_file: str
        Path to file to read from.
    source: pyteomics.mzmlb.MzMLb
        Underlying scan data source
    """
    _parser_cls = _MzMLbParser

    @property
    def has_fast_random_access(self):
        return DefinitelyFastRandomAccess

    @classmethod
    def prebuild_byte_offset_file(cls, path):
        """A stub method. MzMLb does not require an external index.

        Parameters
        ----------
        path : :class:`str` or file-like
            The path to the file to index, or a file-like object with a name attribute.
        """
        return None


def is_mzmlb_file(path):
    '''Detect whether or not the file referenced by ``path``
    is a mzMLb file.

    Parameters
    ----------
    path: :class:`str`
        The path to test

    Returns
    -------
    :class:`bool`:
        Whether or not the file is a mzMLb file.
    '''
    try:
        import h5py
        if _BaseParser == object:
            raise ImportError('pyteomics.mzmlb')
    except ImportError:
        return False
    try:
        source = h5py.File(path, 'r')
        source['mzML']
        return True
    except KeyError:
        return False


def infer_reader(path):
    '''If the file referenced by ``path`` is a mzMLb
    file, return the callable (:class:`MzMLbLoader`) to
    open it, otherwise raise an exception.

    Parameters
    ----------
    path: :class:`str`
        The path to test

    Returns
    -------
    :class:`type`:
        The type to use to open the file

    Raises
    ------
    :class:`ValueError`:
        If the file is not a mzMLb file
    '''
    if is_mzmlb_file(path):
        return MzMLbLoader
    raise ValueError("Not mzMLb File")


def determine_if_available():
    '''Checks whether or not the mzMLb HDF5-based
    file reading feature is available.

    Returns
    -------
    :class:`bool`:
        Whether or not the feature is enabled.
    '''
    try:
        import h5py
        if _BaseParser == object:
            raise ImportError('pyteomics.mzmlb')
        return True
    except (OSError, ImportError):
        return False
