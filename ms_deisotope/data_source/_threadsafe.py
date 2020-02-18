import threading

from .scan.scan import Scan
from .scan.loader import RandomAccessScanSource
from .infer_type import MSFileLoader


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
