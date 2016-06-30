from pyteomics import mzml
from .common import PrecursorInformation, ScanIteratorBase, ChargeNotProvided
from weakref import WeakValueDictionary


def _yield_from_index(self, start=None):
    offset_provider = self._offset_index.offsets
    keys = offset_provider.keys()
    if start is not None:
        if isinstance(start, basestring):
            start = keys.index(start)
        elif isinstance(start, int):
            start = start
        else:
            raise TypeError("Cannot start from object %r" % start)
    else:
        start = 0
    for key in keys[start:]:
        yield self.get_by_id(key)


class PaddedBuffer(object):
    def __init__(self, content, start='<pad>', end="</pad>"):
        self.content = content
        self.start = start
        self.end = end
        self.position = 0
        self.diff = 0

    def read(self, n):
        if self.position < len(self.start):
            out = "".join([self.start, self.content.read(n - len(self.start))])
            self.position += n
            return out
        else:
            out = self.content.read(n)
            if len(out) < n:
                diff = n - len(out)
                if self.diff == 0:
                    self.diff = diff
                    return "".join([out, self.end[:diff]])
                else:
                    return self.end[self.diff:diff]


class MzMLLoader(ScanIteratorBase):
    def __init__(self, mzml_file, use_index=True):
        self.mzml_file = mzml_file
        self._source = mzml.MzML(mzml_file, read_schema=True, iterative=True, use_index=use_index)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._use_index = use_index

    def __reduce__(self):
        return PyteomicsMzMLLoader, (self.mzml_file, self._use_index)

    def reset(self):
        self._make_iterator(None)
        self._scan_cache = WeakValueDictionary()

    def _make_iterator(self, iterator=None):
        self._producer = self._scan_group_iterator(iterator)

    def _scan_group_iterator(self, iterator=None):
        if iterator is None:
            iterator = iter(self._source)
        precursor_scan = None
        product_scans = []

        current_level = 1

        _make_scan = self._make_scan

        for scan in iterator:
            packed = _make_scan(scan)
            self._scan_cache[packed.id] = packed
            if scan['ms level'] == 2:
                if current_level < 2:
                    current_level = 2
                product_scans.append(packed)
            elif scan['ms level'] == 1:
                if current_level > 1:
                    precursor_scan.product_scans = list(product_scans)
                    yield precursor_scan, product_scans
                else:
                    if precursor_scan is not None:
                        precursor_scan.product_scans = list(product_scans)
                        yield precursor_scan, product_scans
                precursor_scan = packed
                product_scans = []
            else:
                raise Exception("This library is not able to handle MS levels higher than 2")

    def next(self):
        return self._producer.next()

    def get_scan_by_id(self, scan_id):
        try:
            return self._scan_cache[scan_id]
        except KeyError:
            packed = self._make_scan(self._source.get_by_id(scan_id))
            self._scan_cache[packed.id] = packed
            return packed

    def start_from_scan(self, scan_id):
        iterator = _yield_from_index(self._source, scan_id)
        self._make_iterator(iterator)

    # Begin ScanDataSourceBase API

    def scan_arrays(self, scan):
        return scan['m/z array'], scan["intensity array"]

    def precursor_information(self, scan):
        pinfo_dict = scan["precursorList"]['precursor'][0]["selectedIonList"]['selectedIon'][0]
        precursor_scan_id = scan["precursorList"]['precursor'][0]['spectrumRef']
        pinfo = PrecursorInformation(
            mz=pinfo_dict['selected ion m/z'],
            intensity=pinfo_dict.get('peak intensity', 0.0),
            charge=pinfo_dict.get('charge state', ChargeNotProvided),
            precursor_scan_id=precursor_scan_id,
            _source=self)
        return pinfo

    def scan_title(self, scan):
        return scan["spectrum title"]

    def scan_id(self, scan):
        return scan["id"]

    def scan_index(self, scan):
        return scan['index']

    def ms_level(self, scan):
        return scan['ms level']

    def scan_time(self, scan):
        return scan['scanList']['scan'][0]['scan start time']


PyteomicsMzMLLoader = MzMLLoader
