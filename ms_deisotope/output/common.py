from ..data_source.common import Scan, ScanBunch, ProcessedScan, PrecursorInformation


class ScanSerializerBase(object):
    def __init__(self, *args, **kwargs):
        pass

    def save_scan_bunch(self, bunch, **kwargs):
        raise NotImplementedError()

    def save(self, bunch, **kwargs):
        self.save_scan_bunch(bunch, **kwargs)


class ScanDeserializerBase(object):
    def __init__(self, *args, **kwargs):
        pass

    def __next__(self):
        raise NotImplementedError()

    def next(self):
        return self.__next__()

    def get_scan_by_id(self, x):
        raise NotImplementedError()

    def get_scan_by_time(self, x):
        raise NotImplementedError()

    def get_scan_by_index(self, x):
        raise NotImplementedError()
