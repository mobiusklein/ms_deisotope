from pyteomics import mzml
from .common import (
    PrecursorInformation, ScanIterator, ScanDataSource, ChargeNotProvided,
    ScanBunch)
from weakref import WeakValueDictionary
from lxml.etree import XMLSyntaxError


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


class MzMLDataInterface(ScanDataSource):
    """Provides implementations of all of the methods needed to implement the
    :class:`ScanDataSource` for mzML files. Not intended for direct instantiation.
    """
    def _scan_arrays(self, scan):
        try:
            return scan['m/z array'], scan["intensity array"]
        except KeyError:
            return mzml.np.array([]), mzml.np.array([])

    def _precursor_information(self, scan):
        pinfo_dict = scan["precursorList"]['precursor'][0]["selectedIonList"]['selectedIon'][0]
        precursor_scan_id = scan["precursorList"]['precursor'][0]['spectrumRef']
        pinfo = PrecursorInformation(
            mz=pinfo_dict['selected ion m/z'],
            intensity=pinfo_dict.get('peak intensity', 0.0),
            charge=pinfo_dict.get('charge state', ChargeNotProvided),
            precursor_scan_id=precursor_scan_id,
            source=self)
        return pinfo

    def _scan_title(self, scan):
        return scan["spectrum title"]

    def _scan_id(self, scan):
        return scan["id"]

    def _scan_index(self, scan):
        return scan['index']

    def _ms_level(self, scan):
        return scan['ms level']

    def _scan_time(self, scan):
        return scan['scanList']['scan'][0]['scan start time']

    def _is_profile(self, scan):
        return "profile spectrum" in scan

    def _polarity(self, scan):
        if "positive scan" in scan:
            return 1
        elif "negative scan" in scan:
            return -1


class MzMLLoader(MzMLDataInterface, ScanIterator):

    __data_interface__ = MzMLDataInterface

    def __init__(self, mzml_file, use_index=True):
        self.mzml_file = mzml_file
        self._source = mzml.MzML(mzml_file, read_schema=True, iterative=True, use_index=use_index)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._use_index = use_index

    def __reduce__(self):
        return MzMLLoader, (self.mzml_file, self._use_index)

    @property
    def index(self):
        return self._source._offset_index

    @property
    def source(self):
        return self._source

    def reset(self):
        self._make_iterator(None)
        self._scan_cache = WeakValueDictionary()

    def _make_iterator(self, iterator=None):
        self._producer = self._scan_group_iterator(iterator)

    def _validate(self, scan_dict):
        return "m/z array" in scan_dict

    def _scan_group_iterator(self, iterator=None):
        if iterator is None:
            iterator = iter(self._source)
        precursor_scan = None
        product_scans = []

        current_level = 1

        _make_scan = self._make_scan

        for scan in iterator:
            if not self._validate(scan):
                continue
            packed = _make_scan(scan)
            self._scan_cache[packed.id] = packed
            if scan['ms level'] == 2:
                if current_level < 2:
                    current_level = 2
                product_scans.append(packed)
            elif scan['ms level'] == 1:
                if current_level > 1:
                    precursor_scan.product_scans = list(product_scans)
                    yield ScanBunch(precursor_scan, product_scans)
                else:
                    if precursor_scan is not None:
                        precursor_scan.product_scans = list(product_scans)
                        yield ScanBunch(precursor_scan, product_scans)
                precursor_scan = packed
                product_scans = []
            else:
                raise Exception("This object is not able to handle MS levels higher than 2")

    def next(self):
        try:
            return self._producer.next()
        except XMLSyntaxError:
            raise StopIteration(
                "This iterator may need to be reset by calling `reset` to continue using it after"
                " using a random-access function like `get_by_id`")

    def __next__(self):
        return self.next()

    def get_scan_by_id(self, scan_id):
        """Retrieve the scan object for the specified scan id. If the
        scan object is still bound and in memory somewhere, a reference
        to that same object will be returned. Otherwise, a new object will
        be created.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        Scan
        """
        try:
            return self._scan_cache[scan_id]
        except KeyError:
            packed = self._make_scan(self._source.get_by_id(scan_id))
            self._scan_cache[packed.id] = packed
            return packed

    def get_scan_by_time(self, time):
        scan_ids = tuple(self.index)
        lo = 0
        hi = len(scan_ids)
        while hi != lo:
            mid = (hi + lo) / 2
            sid = scan_ids[mid]
            scan = self.get_scan_by_id(sid)
            if not self._validate(scan._data):
                sid = scan_ids[mid - 1]
                scan = self.get_scan_by_id(sid)

            scan_time = scan.scan_time
            if scan_time == time:
                return scan
            elif (hi - lo) == 1:
                return scan
            elif scan_time > time:
                hi = mid
            else:
                lo = mid
        if hi == 0 and not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")

    def get_scan_by_index(self, index):
        if not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")
        return self.get_scan_by_id(tuple(self.index)[index])

    def _locate_ms1_scan(self, scan):
        while scan.ms_level != 1:
            if scan.index <= 0:
                raise IndexError("Cannot search backwards with a scan index <= 0 (%r)" % scan.index)
            scan = self.get_scan_by_index(scan.index - 1)
        return scan

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True):
        if scan_id is None:
            if rt is not None:
                scan = self.get_scan_by_time(rt)
            elif index is not None:
                scan = self.get_scan_by_index(index)
            else:
                raise ValueError("Must provide a scan locator, one of (scan_id, rt, index)")

            scan_id = scan.id
        else:
            scan = self.get_scan_by_id(scan_id)

        # We must start at an MS1 scan, so backtrack until we reach one
        if require_ms1:
            scan = self._locate_ms1_scan(scan)
            scan_id = scan.id

        iterator = _yield_from_index(self._source, scan_id)
        self._make_iterator(iterator)
        return self


PyteomicsMzMLLoader = MzMLLoader
