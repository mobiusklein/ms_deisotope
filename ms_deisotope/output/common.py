# pragma: no cover
from ms_deisotope.utils import Base
from ms_deisotope.data_source.common import ScanBunch, ScanIterator


class ScanSerializerBase(object):
    def __init__(self, *args, **kwargs):
        super(ScanSerializerBase, self).__init__(*args, **kwargs)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def save_scan_bunch(self, bunch, **kwargs):
        self.save_scan(bunch.precursor, **kwargs)
        for prod in bunch.products:
            self.save_scan(prod, **kwargs)

    def save_scan(self, scan, **kwargs):
        raise NotImplementedError()

    def save(self, bunch, **kwargs):
        if isinstance(bunch, ScanBunch):
            self.save_scan_bunch(bunch, **kwargs)
        else:
            self.save_scan(bunch, **kwargs)

    def complete(self):
        pass


class ScanDeserializerBase(object):
    def __init__(self, *args, **kwargs):
        super(ScanDeserializerBase, self).__init__(*args, **kwargs)

    def __next__(self):
        return self.next()

    def next(self):
        raise NotImplementedError()

    def get_scan_by_id(self, scan_id):
        raise NotImplementedError()

    def get_scan_by_time(self, rt, require_ms1=False):
        raise NotImplementedError()

    def get_scan_by_index(self, index):
        raise NotImplementedError()


ScanIterator.register(ScanDeserializerBase)


class SampleRun(Base):

    def __init__(self, name, uuid, completed=True, sample_type=None, **kwargs):
        self.name = name
        self.uuid = uuid
        self.sample_type = sample_type
        self.completed = completed
        self.parameters = kwargs
