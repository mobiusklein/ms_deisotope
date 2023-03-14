import os
import io

from typing import Dict, Union, List, Type, Callable

from .mzml import MzMLLoader
from .mzxml import MzXMLLoader
from .mgf import MGFLoader
from .scan.loader import RandomAccessScanSource, ScanIterator
from . import _compression


PathLikeOrFile = Union[io.IOBase, os.PathLike, str]
PathLikeOrFileToTypeFn = Callable[[
    PathLikeOrFile], Union[RandomAccessScanSource, ScanIterator]]


class FormatGuesser(object):
    """Infer which file format a given file path or file-like object.

    Attributes
    ----------
    guessers : List[Callable[[PathLikeOrFile], Type[:class:`~.RandomAccessScanSource`]]]
    reader_types : List[:class:`~.RandomAccessScanSource`]
    """

    guessers: List[PathLikeOrFileToTypeFn]
    reader_types: List[Union[RandomAccessScanSource, ScanIterator]]

    def __init__(self, guessers: List[PathLikeOrFileToTypeFn], reader_types: List[Type]):
        self.guessers = list(guessers or [])
        self.reader_types = list(reader_types or [])

    def register_type_guesser(self, reader_guesser: PathLikeOrFileToTypeFn) -> PathLikeOrFileToTypeFn:
        """
        Add a type determining function to the internal registry.

        This method can be used as a decorator.

        Parameters
        ----------
        reader_guesser : Callable[[Union[io.IOBase, os.PathLike, str], Union[RandomAccessScanSource, ScanIterator]]
            The function predicate that takes a path or file-like object and returns
            an appropriate type for reading it, or raises an error if not compatible.

        Returns
        -------
        Callable[[Union[io.IOBase, os.PathLike, str], Union[RandomAccessScanSource, ScanIterator]]
        """
        self.guessers.append(reader_guesser)
        return reader_guesser

    def guess_type(self, file_path: PathLikeOrFile) -> Type[Union[RandomAccessScanSource, ScanIterator]]:
        """Test the path or file-like objects using each registered predicate, returning the first type that works."""
        for guesser in self.guessers:
            try:
                reader_type = guesser(file_path)
                return reader_type
            except (ValueError, IOError, ImportError, TypeError, AttributeError):
                continue
        raise ValueError("Cannot determine ScanLoader type from %r" % (file_path, ))

    def add_reader_type(self, reader: Union[RandomAccessScanSource, ScanIterator]):
        self.reader_types.append(reader)

    def guess_from_file_extension(self, file_path: PathLikeOrFile, ext_map: Dict[str, Type[Union[RandomAccessScanSource, ScanIterator]]]):
        if hasattr(file_path, 'name'):
            file_path = file_path.name
        if file_path.endswith(".gz"):
            file_path = file_path[:-3]
        ext = os.path.splitext(file_path)[1].lower()
        if ext in ext_map:
            return ext_map[ext]
        else:
            raise ValueError("Cannot determine opener type from file path")

    def add_file_extension_guesser(self, ext_map):
        def guesser(file_path):
            return self.guess_from_file_extension(file_path, ext_map)

        self.register_type_guesser(guesser)
        return guesser

    def open_file(self, file_path: PathLikeOrFile, *args, **kwargs) -> Union[RandomAccessScanSource, ScanIterator]:
        """
        Create an object that reads scans from any supported data file format.

        Provides both iterative and random access.
        """
        reader_type = self.guess_type(file_path)
        try:
            is_file = os.path.isfile(file_path)
        except TypeError:
            is_file = False
        if is_file:
            is_gz_compressed = _compression.test_gzipped(file_path)
            if is_gz_compressed:
                fobj = _compression.get_opener(file_path)
            else:
                fobj = file_path
        else:
            fobj = file_path
        return reader_type(fobj, *args, **kwargs)

    def __call__(self, file_path: PathLikeOrFile, *args, **kwargs) -> Union[RandomAccessScanSource, ScanIterator]:
        """Proxy for :meth:`open_file`.

        Parameters
        ----------
        file_path : str or :class:`io.IOBase`

        Returns
        -------
        :class:`~.RandomAccessScanSource` or :class:`~.ScanIterator`
        """
        return self.open_file(file_path, *args, **kwargs)


MSFileLoader = FormatGuesser([], [MzMLLoader, MzXMLLoader, MGFLoader])

register_type_guesser = MSFileLoader.register_type_guesser
guessers = MSFileLoader.guessers
reader_types = MSFileLoader.reader_types

guess_type_from_path = MSFileLoader.add_file_extension_guesser({
    ".mzml": MzMLLoader,
    ".mzxml": MzXMLLoader,
    ".mgf": MGFLoader
})


@register_type_guesser
def guess_type_from_file_sniffing(file_path):
    is_random_access_file = is_random_access(file_path)
    if is_random_access_file:
        handle = file_path
        header = handle.read(1000)
        handle.seek(0)
    else:
        with open(file_path, 'rb') as handle:
            header = handle.read(1000)

    if _compression.starts_with_gz_magic(header):
        if is_random_access_file:
            handle = _compression.GzipFile(fileobj=file_path, mode='rb')
            header = handle.read(1000)
            handle.seek(0)
        else:
            with _compression.GzipFile(file_path, mode='rb') as handle:
                header = handle.read(1000)

    if b"mzML" in header:
        return MzMLLoader
    elif b"mzXML" in header:
        return MzXMLLoader
    elif b"BEGIN IONS" in header:
        return MGFLoader
    else:
        raise ValueError("Cannot determine ScanLoader type from header")


guess_type = MSFileLoader.guess_type


def is_random_access(fp):
    return hasattr(fp, 'read') and hasattr(fp, 'seek')


try:
    from .mzmlb import MzMLbLoader, infer_reader as _infer_mzmlb, determine_if_available

    if determine_if_available():
        reader_types.append(MzMLbLoader)
        register_type_guesser(_infer_mzmlb)
except ImportError:
    pass


try:
    from .thermo_raw_net import (
        ThermoRawLoader as ThermoRawNetLoader, infer_reader as _check_is_thermo_raw_net,
        register_dll as register_thermo_net_dll)
    reader_types.append(ThermoRawNetLoader)
    register_type_guesser(_check_is_thermo_raw_net)

except ImportError:  # pragma: no cover
    def register_thermo_net_dll(*args, **kwargs):
        pass

try:
    from .thermo_raw import (
        ThermoRawLoader, infer_reader as _check_is_thermo_raw,
        register_dll as register_thermo_dll)

    reader_types.append(ThermoRawLoader)
    register_type_guesser(_check_is_thermo_raw)

except ImportError:  # pragma: no cover
    def register_thermo_dll(*args, **kwargs):
        pass


try:
    from .agilent_d import (
        AgilentDLoader,
        infer_reader as _check_is_agilent_d,
        register_dll_dir as register_agilent_dll_dir)

    reader_types.append(AgilentDLoader)
    register_type_guesser(_check_is_agilent_d)
except ImportError:  # pragma: no cover
    def register_agilent_dll_dir(*args, **kwargs):
        pass


try:
    from .masslynx import (
        MassLynxRawLoader, infer_reader as _check_mass_lynx_raw,
        register_waters_masslynx_dll)

    reader_types.append(MassLynxRawLoader)
    register_type_guesser(_check_mass_lynx_raw)
except ImportError:
    def register_waters_masslynx_dll(*args, **kwargs):
        pass
