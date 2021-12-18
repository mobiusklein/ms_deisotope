from ms_deisotope.data_source.infer_type import FormatGuesser, is_random_access
from ms_deisotope.data_source import _compression
from ms_deisotope.data_source.mzmlb import is_mzmlb_file as _is_mzmlb_file

from ms_deisotope.output.mzml import ProcessedMzMLLoader
from ms_deisotope.output.mzmlb import ProcessedMzMLbLoader
from ms_deisotope.output.mgf import ProcessedMGFLoader

ProcessedMSFileLoader = FormatGuesser([], [ProcessedMzMLLoader, ProcessedMzMLbLoader, ProcessedMGFLoader])

ProcessedMSFileLoader.add_file_extension_guesser({
    ".mzml": ProcessedMzMLLoader,
    ".mzmlb": ProcessedMzMLbLoader,
    ".mgf": ProcessedMGFLoader
})


@ProcessedMSFileLoader.register_type_guesser
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
        return ProcessedMzMLLoader
    elif b"BEGIN IONS" in header:
        return ProcessedMGFLoader
    else:
        raise ValueError("Cannot determine ScanLoader type from header")


@ProcessedMSFileLoader.register_type_guesser
def is_mzmlb_file(file_path):
    if _is_mzmlb_file(file_path):
        return ProcessedMzMLbLoader
    raise ValueError("Not an mzMLb file")
