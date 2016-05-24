from pyteomics import mzml
from .common import Scan, PrecursorInformation
from weakref import WeakValueDictionary


class PyteomicsMzMLLoader(object):
    def __init__(self, mzml_file, use_index=False):
        self.mzml_file = mzml_file
        self._source = mzml.MzML(mzml_file, read_schema=True, iterative=True, use_index=use_index)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._use_index = use_index

    def __reduce__(self):
        return PyteomicsMzMLLoader, (self.mzml_file, self._use_index)

    def scan_arrays(self, scan):
        return scan['m/z array'], scan["intensity array"]

    def precursor_information(self, scan):
        pinfo_dict = scan["precursorList"]['precursor'][0]["selectedIonList"]['selectedIon'][0]
        precursor_scan_id = scan["precursorList"]['precursor'][0]['spectrumRef']
        pinfo = PrecursorInformation(
            mz=pinfo_dict['selected ion m/z'],
            intensity=pinfo_dict['peak intensity'],
            charge=pinfo_dict['charge state'],
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

    def _make_scan(self, data):
        return Scan(data, self)

    def _scan_group_iterator(self):
        precursor_scan = None
        product_scans = []

        current_level = 1

        _make_scan = self._make_scan

        for scan in self._source:
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

    __next__ = next

    def __iter__(self):
        return self
