import os
from .mzml import MzMLLoader
from .mzxml import MzXMLLoader


def guess_type_from_path(file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.mzml':
        return MzMLLoader
    elif ext == '.mzxml':
        return MzXMLLoader
    else:
        raise ValueError("Cannot determine ScanLoader type from file path")


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
    try:
        return guess_type_from_path(file_path)
    except ValueError:
        try:
            return guess_type_from_file_sniffing(file_path)
        except ValueError:
            raise ValueError("Cannot determine ScanLoader type")
