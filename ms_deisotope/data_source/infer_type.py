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


def guess_type(file_path):
    for guesser in guessers:
        try:
            reader_type = guesser(file_path)
            return reader_type
        except (ValueError, IOError, ImportError, TypeError, AttributeError):
            continue
    raise ValueError("Cannot determine ScanLoader type from %r" % (file_path, ))


def MSFileLoader(file_path, *args, **kwargs):
    """Factory function to create an object that reads scans from
    any supported data file format. Provides both iterative and
    random access.
    """
    reader_type = guess_type(file_path)
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
