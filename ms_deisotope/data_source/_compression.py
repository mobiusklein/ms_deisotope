import io
import gzip

try:
    import idzip
    GzipFile = idzip.IdzipFile
except (ImportError, AttributeError):
    GzipFile = gzip.GzipFile


DEFAULT_BUFFER_SIZE = int(2e6)


def test_gzipped(f):
    """Checks the first two bytes of the
    passed file for gzip magic numbers

    Parameters
    ----------
    f : file-like or path-like

    Returns
    -------
    bool
    """
    if isinstance(f, basestring):
        f = io.open(f, 'rb')
    current = f.tell()
    f.seek(0)
    magic = f.read(2)
    f.seek(current)
    return magic == b'\037\213'


def starts_with_gz_magic(bytestring):
    return bytestring.startswith(b'\037\213')


def get_opener(f, buffer_size=None):
    if buffer_size is None:
        buffer_size = DEFAULT_BUFFER_SIZE
    if not hasattr(f, 'read'):
        f = io.open(f, 'rb')
    buffered_reader = io.BufferedReader(f, buffer_size)
    if test_gzipped(f):
        handle = GzipFile(fileobj=buffered_reader, mode='rb')
    else:
        handle = buffered_reader
    return handle
