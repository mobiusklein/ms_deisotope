import os
from .mzml import MzMLLoader
from .mzxml import MzXMLLoader
from .mgf import MGFLoader
from . import _compression

guessers = []
reader_types = [MzMLLoader, MzXMLLoader, MGFLoader]


def register_type_guesser(reader_guesser):
    guessers.append(reader_guesser)
    return reader_guesser


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


@register_type_guesser
def guess_type_from_path(file_path):
    if hasattr(file_path, 'name'):
        file_path = file_path.name
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.mzml':
        return MzMLLoader
    elif ext == '.mzxml':
        return MzXMLLoader
    elif ext == '.mgf':
        return MGFLoader
    else:
        raise ValueError("Cannot determine ScanLoader type from file path")


@register_type_guesser
def guess_type_from_file_sniffing(file_path):
    with open(file_path, 'rb') as handle:
        header = handle.read(1000)

    if _compression.starts_with_gz_magic(header):
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


def guess_type(file_path):
    for guesser in guessers:
        try:
            reader_type = guesser(file_path)
            return reader_type
        except (ValueError, IOError, ImportError, TypeError, AttributeError):
            continue
    raise ValueError("Cannot determine ScanLoader type")


def MSFileLoader(file_path, *args, **kwargs):
    """Factory function to create an object that reads scans from
    any supported data file format. Provides both iterative and
    random access.
    """
    reader_type = guess_type(file_path)
    is_gz_compressed = _compression.test_gzipped(file_path)
    if is_gz_compressed:
        fobj = _compression.get_opener(file_path)
    else:
        fobj = file_path
    return reader_type(fobj, *args, **kwargs)
