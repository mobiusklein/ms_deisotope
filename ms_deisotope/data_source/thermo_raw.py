# pragma: no cover
import re

from weakref import WeakValueDictionary
from collections import OrderedDict

import logging

import numpy as np

from ms_deisotope.data_source.common import (
    ScanDataSource, ScanIterator, RandomAccessScanSource,
    Scan, PrecursorInformation, ScanBunch, ChargeNotProvided,
    ActivationInformation)

from ms_deisotope.utils import Base


try:
    from ms_deisotope.data_source._vendor.MSFileReader import (
        ThermoRawfile as _ThermoRawFileAPI, register_dll,
        log as _api_logger)

    comtypes_logger = logging.getLogger("comtypes")
    comtypes_logger.setLevel("INFO")
    _api_logger.setLevel("INFO")

    def is_thermo_raw_file(path):
        try:
            _ThermoRawFileAPI(path)
            return True
        except (WindowsError, IOError, ImportError):
            return False

    def infer_reader(path):
        if is_thermo_raw_file(path):
            return ThermoRawLoader
        raise ValueError("Not Thermo Raw File")

    def determine_if_available():
        try:
            _ThermoRawFileAPI.create_com_object()
            return True
        except ImportError:
            return False


except ImportError as e:
    message = str(e)

    def is_thermo_raw_file(path):
        return False

    def infer_reader(path):
        raise ValueError(message)

    def register_dll(paths):
        print("no-op: %s" % (message,))
        return False

    def determine_if_available():
        print("no-op: %s" % (message,))
        return False

try:
    range = xrange
except NameError:
    pass


analyzer_pat = re.compile(r"(?P<mass_analyzer_type>ITMS|TQMS|SQMS|TOFMS|FTMS|SECTOR)")
polarity_pat = re.compile(r"(?P<polarity>[\+\-])")
point_type_pat = re.compile(r"(?P<point_type>[CP])")
ionization_pat = re.compile(r"(?P<ionization_type>EI|CI|FAB|APCI|ESI|APCI|NSI|TSP|FD|MALDI|GD)")
scan_type_pat = re.compile(r"(?P<scan_type>FULL|SIM|SRM|CRM|Z|Q1MS|Q3MS)")
activation_pat = re.compile(
    r"""ms(?P<ms_level>\d*)\s
        (?:(?P<isolation_mz>\d+\.\d*)@
        (?P<activation_type>[a-z]+)
        (?P<activation_energy>\d*\.?\d*))?""", re.VERBOSE)


class FilterLine(str):
    def __init__(self, value):
        self.data = filter_line_parser(self)

    def get(self, key):
        return self.data.get(key)


def filter_line_parser(line):
    words = line.upper().split(" ")
    values = dict()
    i = 0
    activation_info = activation_pat.search(line)
    if activation_info is not None:
        activation_info = activation_info.groupdict()
        if activation_info['ms_level'] != "":
            values["ms_level"] = int(activation_info['ms_level'])
            values["isolation_mz"] = float(activation_info['isolation_mz'])
            values["activation_type"] = activation_info['activation_type']
            values["activation_energy"] = float(activation_info['activation_energy'])
    try:
        word = words[i]
        i += 1
        analyzer_info = analyzer_pat.search(word)
        if analyzer_info is not None:
            values['analyzer'] = analyzer_info.group(0)
            word = words[i]
            i += 1
        polarity_info = polarity_pat.search(word)
        if polarity_info is not None:
            polarity_sigil = polarity_info.group(0)
            if polarity_sigil == "+":
                polarity = 1
            elif polarity_sigil == "-":
                polarity = -1
            else:
                polarity = 0
            values["polarity"] = polarity
            word = words[i]
            i += 1

        return values
    except IndexError:
        return values


_id_template = "controllerType=0 controllerNumber=1 scan="


class ThermoRawScanPtr(Base):
    def __init__(self, scan_number):
        self.scan_number = scan_number
        self.filter_line = None


def _make_id(scan_number):
    return "%s%d" % (_id_template, (scan_number))


def _parse_id(scan_id):
    return int(scan_id.replace(_id_template, ""))


class ThermoRawDataInterface(ScanDataSource):
    def _scan_index(self, scan):
        return scan.scan_number - 1

    def _scan_id(self, scan):
        return _make_id(scan.scan_number)

    def _scan_time(self, scan):
        return self._source.RTFromScanNum(
            scan.scan_number)

    def _ms_level(self, scan):
        return self._source.GetMSOrderForScanNum(
            scan.scan_number)

    def _is_profile(self, scan):
        return self._source.IsProfileScanForScanNum(
            scan.scan_number)

    def _polarity(self, scan):
        filter_line = self._filter_line(scan)
        return filter_line.data['polarity']

    def _filter_line(self, scan):
        if scan.filter_line is None:
            scan.filter_line = FilterLine(self._source.GetFilterForScanNum(scan.scan_number))
        return scan.filter_line

    def _scan_title(self, scan):
        return "scan=%d" % scan.scan_number

    def _scan_arrays(self, scan):
        arrays, flags = self._source.GetMassListFromScanNum(
            scan.scan_number)
        mz, intensity = arrays
        return np.array(mz), np.array(intensity)

    def _precursor_information(self, scan):
        if self._ms_level(scan) == 1:
            return None
        scan_number = scan.scan_number
        pinfo_struct = self._source.GetPrecursorInfoFromScanNum(scan_number)
        precursor_scan_number = None
        labels, flags, _ = self._source.GetAllMSOrderData(scan_number)
        if pinfo_struct:
            mz = pinfo_struct.monoIsoMass
            charge = pinfo_struct.chargeState
            intensity = float(labels.intensity[0])
            precursor_scan_number = pinfo_struct.scanNumber
        else:
            mz = labels.mass[0]
            intensity = float(labels.intensity[0])
            charge = labels.charge[0]
        if not charge:
            charge = ChargeNotProvided
        trailer = self._source.GetTrailerExtraForScanNum(scan_number)
        _mz = trailer.get('Monoisotopic M/Z', 0.0)
        if _mz > 0:
            mz = _mz
        _charge = trailer.get('Charge State', 0)
        if _charge != 0:
            charge = _charge
        # Guess which previous scan was the precursor by iterating
        # backwards until a scan is found with a lower MS level
        if precursor_scan_number is None:
            last_index = self._scan_index(scan) - 1
            current_level = self._ms_level(scan)
            i = 0
            while last_index > 0 and i < 100:
                prev_scan = self.get_scan_by_index(last_index)
                if prev_scan.ms_level >= current_level:
                    last_index -= 1
                else:
                    precursor_scan_number = prev_scan._data.scan_number
                    break
                i += 1
        if intensity is None:
            intensity = 0.0
        if mz is None:
            mz = 0.0
        if charge is None or charge == 0:
            charge = ChargeNotProvided
        pinfo = PrecursorInformation(
            mz, intensity, charge, _make_id(precursor_scan_number),
            0, 0, 0,
            product_scan_id=_make_id(scan.scan_number))
        return pinfo

    def _activation(self, scan):
        filter_line = self._filter_line(scan)
        activation_type = filter_line.get("activation_type")
        if activation_type is not None:
            energy = filter_line.get("activation_energy")
            return ActivationInformation(activation_type, energy)
        return None


class ThermoRawLoader(ThermoRawDataInterface, RandomAccessScanSource, ScanIterator):
    def __init__(self, source_file, **kwargs):
        self.source_file = source_file
        self._source = _ThermoRawFileAPI(self.source_file)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._index = self._pack_index()
        self._first_scan_time = self.get_scan_by_index(0).scan_time
        self._last_scan_time = self.get_scan_by_id(self._source.LastSpectrumNumber).scan_time

    def __reduce__(self):
        return self.__class__, (self.source_file,)

    @property
    def index(self):
        return self._index

    def __repr__(self):
        return "ThermoRawLoader(%r)" % (self.source_file)

    def _pack_index(self):
        index = OrderedDict()
        for sn in range(1, self._source.NumSpectra):
            index[_make_id(sn)] = sn
        return index

    def close(self):
        self._source.Close()

    def reset(self):
        self.make_iterator(None)
        self._scan_cache = WeakValueDictionary()

    def get_scan_by_id(self, scan_id):
        """Retrieve the scan object for the specified scan id.

        If the scan object is still bound and in memory somewhere,
        a reference to that same object will be returned. Otherwise,
        a new object will be created.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        Scan
        """
        scan_number = int(str(scan_id).replace(_id_template, ''))
        try:
            return self._scan_cache[scan_number]
        except KeyError:
            package = ThermoRawScanPtr(scan_number)
            scan = Scan(package, self)
            self._scan_cache[scan_number] = scan
            return scan

    def get_scan_by_index(self, index):
        """Retrieve the scan object for the specified scan index.

        This internally calls :meth:`get_scan_by_id` which will
        use its cache.

        Parameters
        ----------
        index: int
            The index to get the scan for

        Returns
        -------
        Scan
        """
        scan_number = int(index) + 1
        try:
            return self._scan_cache[scan_number]
        except KeyError:
            package = ThermoRawScanPtr(scan_number)
            scan = Scan(package, self)
            self._scan_cache[scan_number] = scan
            return scan

    def get_scan_by_time(self, time):
        """Retrieve the scan object for the specified scan time.

        This internally calls :meth:`get_scan_by_id` which will
        use its cache.

        Parameters
        ----------
        time : float
            The time to get the nearest scan from

        Returns
        -------
        Scan
        """
        if time < self._first_scan_time:
            time = self._first_scan_time
        elif time > self._last_scan_time:
            time = self._last_scan_time
        scan_number = self._source.ScanNumFromRT(time)
        try:
            return self._scan_cache[scan_number]
        except KeyError:
            package = ThermoRawScanPtr(scan_number)
            scan = Scan(package, self)
            self._scan_cache[scan_number] = scan
            return scan

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True):
        if scan_id is not None:
            scan_number = int(str(scan_id).replace(_id_template, '')) - 1
        elif index is not None:
            scan_number = int(index)
        elif rt is not None:
            start_index = self._source.ScanNumFromRT(rt)
            if require_ms1:
                while start_index != 0:
                    scan = self.get_scan_by_index(start_index)
                    if scan.ms_level > 1:
                        start_index -= 1
                    else:
                        break
                scan_number = start_index
        self._producer = self._scan_group_iterator(
            self._make_scan_index_producer(
                start_index=scan_number))
        return self

    def make_iterator(self, iterator=None, grouped=True):
        if grouped:
            self._producer = self._scan_group_iterator(iterator)
        else:
            self._producer = self._single_scan_iterator(iterator)

    def _make_scan_index_producer(self, start_index=None, start_time=None):
        if start_index is not None:
            return range(start_index + 1, self._source.NumSpectra + 1)
        elif start_time is not None:
            start_index = self._source.ScanNumFromRT(start_time)
            while start_index != 1:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            return range(start_index, self._source.NumSpectra + 1)
        else:
            return range(1, self._source.NumSpectra + 1)

    def _single_scan_iterator(self, iterator=None):
        if iterator is None:
            iterator = self._make_scan_index_producer()
        for ix in iterator:
            packed = self.get_scan_by_id(ix)
            self._scan_cache[packed._data.scan_number] = packed
            yield packed

    def _scan_group_iterator(self, iterator=None):
        if iterator is None:
            iterator = self._make_scan_index_producer()

        precursor_scan = None
        product_scans = []

        current_level = 1

        for ix in iterator:
            packed = self.get_scan_by_id(ix)
            self._scan_cache[packed._data.scan_number] = packed
            if packed.ms_level == 2:
                if current_level < 2:
                    current_level = 2
                product_scans.append(packed)
            elif packed.ms_level == 1:
                if current_level > 1 and precursor_scan is not None:
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
        if precursor_scan is not None:
            yield ScanBunch(precursor_scan, product_scans)

    def next(self):
        return next(self._producer)
