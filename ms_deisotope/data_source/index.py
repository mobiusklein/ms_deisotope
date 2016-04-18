import sqlitedict

from common import Scan, PrecursorInformation
from mzml import PyteomicsMzMLLoader


class IndexedScanLoader(object):
    def __init__(self, path, journal_mode='OFF'):
        self.store = sqlitedict.open(path, journal_mode=journal_mode)
        self.scan_index = self.store.get("_index_", [])
        self._iter_index = 0

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

    def ms_level(self, scan):
        return scan['ms level']

    def _make_scan(self, data):
        return Scan(data, self)

    def next(self):
        try:
            scan = self.store[self.scan_index[self._iter_index]]
            self._iter_index += 1
            return scan
        except IndexError:
            raise StopIteration()

    __next__ = next

    def __iter__(self):
        self._iter_index = 0
        return self

    def get_scan_by_id(self, scan_id):
        return self.store[scan_id]

    def _scan_group_iterator(self):
        precursor_scan = None
        product_scans = []

        current_level = 1

        _make_scan = self._make_scan

        for scan in self:
            if scan['ms level'] == 2:
                if current_level < 2:
                    current_level = 2
                product_scans.append(_make_scan(scan))
            elif scan['ms level'] == 1:
                if current_level > 1:
                    precursor_scan.product_scans = list(product_scans)
                    yield precursor_scan, product_scans
                else:
                    if precursor_scan is not None:
                        precursor_scan.product_scans = list(product_scans)
                        yield precursor_scan, product_scans
                precursor_scan = _make_scan(scan)
                product_scans = []
            else:
                raise Exception("This library is not able to handle MS levels higher than 2")