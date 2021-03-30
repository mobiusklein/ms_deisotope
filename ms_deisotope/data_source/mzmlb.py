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
