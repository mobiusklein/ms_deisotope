import os
from .mzml import MzMLLoader
from .mzxml import MzXMLLoader


guessers = []


def register_type_guesser(reader_guesser):
    guessers.append(reader_guesser)
    return reader_guesser


try:
    from .thermo_raw import (
        ThermoRawLoader, infer_reader as check_is_thermo_raw,
        register_dll as register_thermo_dll)

    register_type_guesser(check_is_thermo_raw)

except ImportError:
    pass


@register_type_guesser
def guess_type_from_path(file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.mzml':
        return MzMLLoader
    elif ext == '.mzxml':
        return MzXMLLoader
    else:
        raise ValueError("Cannot determine ScanLoader type from file path")


@register_type_guesser
def guess_type_from_file_sniffing(file_path):
    with open(file_path, 'rb') as handle:
        header = handle.read(1000)
        if b"mzML" in header:
            return MzMLLoader
        elif b"mzXML" in header:
            return MzXMLLoader
        else:
            raise ValueError("Cannot determine ScanLoader type from header")


def guess_type(file_path):
    for guesser in guessers:
        try:
            reader_type = guesser(file_path)
            return reader_type
        except (ValueError, IOError):
            continue
    raise ValueError("Cannot determine ScanLoader type")

    # try:
    #     return guess_type_from_path(file_path)
    # except ValueError:
    #     try:
    #         return guess_type_from_file_sniffing(file_path)
    #     except ValueError:
    #         raise ValueError("Cannot determine ScanLoader type")


def MSFileLoader(file_path, *args, **kwargs):
    """Factory function to create an object that reads scans from
    any supported data file format. Provides both iterative and
    random access.
    """
    reader_type = guess_type(file_path)
    return reader_type(file_path, *args, **kwargs)
