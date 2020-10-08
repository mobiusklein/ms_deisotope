import threading
import logging

from .scan.scan import Scan
from .scan.loader import RandomAccessScanSource
from .infer_type import MSFileLoader
from ._compression import get_opener


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class PerThreadFileHandle(object):
    def __init__(self, path, mode='rb', encoding=None, opener=None):
        if opener is None:
            opener = get_opener
        self.local = threading.local()
        self._path = path
        self._mode = mode
        self._encoding = encoding
        self._opener = opener

    @property
    def file(self):
        try:
            return self.local.file
        except AttributeError:
            return self._open()

    def __len__(self):
        return len(self.file)

    def __iter__(self):
        return iter(self.file)

    def __getitem__(self, i):
        return self.file[i]

    def __reduce__(self):
        return self.__class__, (self._path, self._mode, self._encoding, self._opener)

    def _open(self):
        self.local.file = self._opener(self._path)
        logger.debug("Opening %r on Thread %r", self._path,
                     threading.current_thread())
        return self.local.file

    def __getattr__(self, attrib):
        if attrib in ('local', '_path', '_mode', '_encoding', '_opener'):
            raise AttributeError(attrib)
        return object.__getattribute__(self.file, attrib)



class ThreadsafeScanSource(object):

    def __init__(self, source_file, opener=None):
        if opener is None:
            opener = MSFileLoader
        self._source_file_name = source_file
        self._opener = opener
        self._threadlocal_store = threading.local()

    def _reader(self):
        try:
            return self._threadlocal_store.source
        except AttributeError:
            self._threadlocal_store.source = self._patch(self._opener(self._source_file_name))
            return self._threadlocal_store.source

    def _make_scan(self, data):
        return Scan(data, self).bind(self)

    def _patch(self, reader):
        '''Replace the _make_scan method of `reader`
        '''
        reader._make_scan = self._make_scan
        return reader

    def __getattr__(self, name):
        reader = self._reader()
        value = getattr(reader, name)
        return value

    def __reduce__(self):
        return self.__class__, (self._source_file_name, self._opener)
