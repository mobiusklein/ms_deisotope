import io

DEFAULT_BUFFER_SIZE = 2 ** 18


class PreBufferedStreamReader(io.IOBase):
    '''A file-like object that can wrap a non-rewindable (e.g. un-seekable)
    file stream like :obj:`sys.stdin` and support seeking over the first **N**
    bytes, but blocking seek operations once the stream has been read beyond
    that point.

    This type exists to adapt :obj:`sys.stdin` or socket files to work with the
    behavior of readers like:class:`~ms_deisotope.data_source.mzml.MzMLLoader`
    that rely on early random access to extract metadata even when full random
    access indexing isn't enabled.

    Attributes
    ----------
    stream : file-like object
        The object to actually read from. Must support :meth:`io.IOBase.read`,
        :meth:`io.IOBase
    '''
    def __init__(self, stream, buffer_size=DEFAULT_BUFFER_SIZE):
        self.stream = stream
        dat = self.stream.read(buffer_size)
        self.buffer = io.BytesIO(dat)
        self._cursor = 0
        self._total_cursor = len(dat)
        self._buffer_size = buffer_size

    def read(self, n=-1):
        if n == -1:
            buffer = bytearray()
            buffer.extend(self.buffer.read())
            buffer.extend(self.stream.read())
            self._cursor = self._buffer_size
            self._total_cursor = len(buffer) - self._cursor
            return bytes(buffer)
        buffer = bytearray(n)
        view = memoryview(buffer)
        k = self.buffer.readinto(view)
        self._cursor += k
        j = 0
        if k < n:
            j = self.stream.readinto(view[k:])
            self._total_cursor += j
            k += j
        return bytes(view[:k])

    def fileno(self):
        return self.stream.fileno()

    def seek(self, offset, whence=io.SEEK_SET):
        if (whence == io.SEEK_SET and offset > self._buffer_size):
            raise io.UnsupportedOperation(
                "Cannot seek beyond the pre-loaded buffer")
        elif whence == io.SEEK_CUR or whence == io.SEEK_END:
            raise io.UnsupportedOperation(
                "Cannot seek not from the start of the file")
        if self._total_cursor > self._buffer_size:
            raise io.UnsupportedOperation(
                "Already read beyond the pre-loaded buffer")
        else:
            self.buffer.seek(offset, whence)
            self._cursor = offset

    def truncate(self, size=None):
        raise io.UnsupportedOperation("Read-only")

    def close(self):
        self.buffer.close()
        return self.stream.close()

    def closed(self):
        return self.stream.closed()

    def seekable(self):
        return True

    def tell(self):
        if self._cursor < self._buffer_size:
            return self._cursor
        return self._total_cursor

    def writable(self):
        return False

    def __enter__(self):
        self.stream.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()