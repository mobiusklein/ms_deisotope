'''Helpers for dealing with generic compressed files.
'''

import io
import gzip

from six import string_types as basestring, PY2

from ms_deisotope.utils import Constant

GzipFile = _GzipFile = gzip.GzipFile

slow_random_access_file_types = [_GzipFile]
try:
    import bz2
    slow_random_access_file_types.append(bz2.BZ2File)
except ImportError:
    pass


try:
    import idzip

    class IdzipFile(idzip.IdzipFile):
        '''Fixes io.BufferedWriter compatibility until merged upstream
        '''
        def write(self, b):
            self._check_can_write()
            value = self._impl.write(b)
            return value

    idzip.compressor.IdzipWriter.enforce_extension = False

    GzipFile = IdzipFile
    WRITE_BUFFER_SIZE = idzip.MAX_MEMBER_SIZE
    has_idzip = True
except (ImportError, AttributeError):
    GzipFile = gzip.GzipFile
    WRITE_BUFFER_SIZE = 2 ** 16
    has_idzip = False

if PY2:
    file_like_object_bases = (file, io.IOBase) # pylint: disable=undefined-variable
else:
    file_like_object_bases = (io.IOBase, )


def test_if_file_has_fast_random_access(file_obj):
    """Determine whether or not the passed file-like
    object supports fast random access.

    Parameters
    ----------
    file_obj : :class:`io.IOBase`
        The file-like object to test

    Returns
    -------
    bool
    """
    # If we have a normal gzip.GzipFile, then we definitely don't have fast random access
    if isinstance(file_obj, tuple(slow_random_access_file_types)):
        return DefinitelyNotFastRandomAccess
    # If we have an idzip.IdzipFile, then we need to query its _impl attribute (QUESTION: Should this be made
    # a part of the IdzipFile API and be pushed upstream?)
    elif has_idzip and isinstance(file_obj, GzipFile):
        # If the _impl attribute is an ordinary gzip.GzipFile, then it is just a normal gzip file without
        # an index.
        if isinstance(file_obj._impl, _GzipFile):
            return DefinitelyNotFastRandomAccess
        # Otherwise it's an idzip file and we can use fast random access
        return DefinitelyFastRandomAccess
    else:
        # We're looking at a file-like object of some sort. It could be a compressed file not caught
        # by the earlier checks. The only good test would be to examine the file's raw contents, but
        # this is not an option here. Assume that we're looking at an uncompressed stream.
        if isinstance(file_obj, file_like_object_bases):
            return DefinitelyFastRandomAccess
        # It's not a file of any sort. It could be a regular file, it could be something else entirely.
        return MaybeFastRandomAccess


DEFAULT_BUFFER_SIZE = int(2e6)
GZIP_MAGIC = b'\037\213'

DefinitelyNotFastRandomAccess = Constant("DefinitelyNotFastRandomAccess", False)
MaybeFastRandomAccess = Constant("MaybeFastRandomAccess", True)
DefinitelyFastRandomAccess = Constant("DefinitelyFastRandomAccess", True)


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
    return magic == GZIP_MAGIC


def starts_with_gz_magic(bytestring):
    '''Tests whether or not a byte string starts with
    the GZIP magic bytes.

    Parameters
    ----------
    bytestring : bytes
        The bytes to test.

    Returns
    -------
    bool
    '''
    return bytestring.startswith(GZIP_MAGIC)


ZSTD_MAGIC = b'(\xb5/\xfd'

def starts_with_zstd_magic(bytestring):
    '''Tests whether or not a byte string starts with
    the ZSTD magic bytes.

    Parameters
    ----------
    bytestring : bytes
        The bytes to test.

    Returns
    -------
    bool
    '''
    return bytestring.startswith(ZSTD_MAGIC)


def get_opener(f, buffer_size=None):
    '''Select the file reading type for the given path or stream.

    Detects whether the file is gzip encoded.
    '''
    if buffer_size is None:
        buffer_size = DEFAULT_BUFFER_SIZE
    if not hasattr(f, 'read'):
        f = io.open(f, 'rb')
    # On Py2, dill doesn't behave correctly with io-derived objects, so we have to
    # patch it below. Don't try to wrap an io.TextIOWrapper on Py3.
    if not isinstance(f, io.BufferedReader) and not isinstance(f, io.TextIOWrapper):
        buffered_reader = io.BufferedReader(f, buffer_size)
    else:
        buffered_reader = f
    if test_gzipped(buffered_reader):
        handle = GzipFile(fileobj=buffered_reader, mode='rb')
    else:
        handle = buffered_reader
    return handle



if PY2:
    # Begin Monkeypatching Dill
    import dill

    @dill.register(io.BufferedReader)
    def save_file(pickler, obj):
        '''Monkeypatch Py2 Dill
        '''
        dill._dill._save_file(pickler, obj, io.open)
