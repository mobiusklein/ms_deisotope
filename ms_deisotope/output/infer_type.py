import io
import os

from ms_deisotope.data_source.infer_type import FormatGuesser, is_random_access
from ms_deisotope.data_source import _compression
from ms_deisotope.data_source.mzmlb import is_mzmlb_file as _is_mzmlb_file

from ms_deisotope.output.mzml import ProcessedMzMLLoader, MzMLSerializer
from ms_deisotope.output.mzmlb import ProcessedMzMLbLoader, MzMLbSerializer
from ms_deisotope.output.mgf import ProcessedMGFLoader, MGFSerializer

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


serializers = {}

for tp in [MzMLSerializer, MzMLbSerializer, MGFSerializer]:
    for ext in tp.file_extensions:
        serializers[ext] = tp


def get_writer(filename, **kwargs):
    '''Open a writer for a provided filename, inferring the format from
    the file extension.

    .. warning::
        If using a file-like object, do not use a compressed writer or else
        the stream will be doubly-compressed.

    Parameters
    ----------
    filename : :class:`str`, :class:`os.PathLike`, or file-like object
        The path to the file to open, or a file-like object with a "name"
        attribute.
    **kwargs
        Keyword arguments forwarded to the writer

    Returns
    -------
    MzMLSerializer or MGFSerializer or MzMLbSerializer
    '''
    if hasattr(filename, 'name'):
        handle = filename
        name = handle.name
    else:
        handle = None
        name = os.path.basename(filename)

    is_gzipped = False
    if name.lower().endswith(".gz"):
        is_gzipped = True
        name = name[:-3]

    name, ext = os.path.splitext(name.lower())
    serializer_cls = serializers[ext[1:]]

    if serializer_cls == MzMLbSerializer:
        return serializer_cls(filename, **kwargs)

    if is_gzipped:
        if handle is None:
            handle = _compression.GzipFile(filename, 'wb')
        else:
            handle = _compression.GzipFile(fileobj=filename, mode='wb')

    else:
        if handle is None:
            handle = open(filename, 'wb')

    return serializer_cls(handle, **kwargs)
