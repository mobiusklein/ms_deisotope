# pragma: no cover
import re
import os

from weakref import WeakValueDictionary
from collections import OrderedDict, defaultdict

import logging

import numpy as np

from ms_deisotope.data_source.common import (
    ScanDataSource, ScanIterator, RandomAccessScanSource,
    Scan, PrecursorInformation, ScanBunch, ChargeNotProvided,
    ActivationInformation, IsolationWindow,
    component, ComponentGroup, InstrumentInformation,
    FileInformation, SourceFile)

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


analyzer_map = {
    'FTMS': component("orbitrap"),
    "ITMS": component("ion trap"),
    "SQMS": component("quadrupole"),
    "TQMS": component("quadrupole"),
    "TOFMS": component("time-of-flight"),
    "SECTOR": component("magnetic sector")
}


ionization_map = {
    "EI": component("electron ionization"),
    "CI": component("chemical ionization"),
    "FAB": component("fast atom bombardment ionization"),
    "ESI": component("electrospray ionization"),
    "NSI": component("nanoelectrospray"),
    "APCI": component("atmospheric pressure chemical ionization"),
    "TSP": component("thermospray ionization"),
    "FD": component("field desorption"),
    "MALDI": component("matrix assisted laser desorption ionization"),
    "GD": component("glow discharge ionization"),
}


inlet_map = {
    "FAB": component("continuous flow fast atom bombardment"),
    "ESI": component("electrospray inlet"),
    "NSI": component("nanospray inlet"),
    "TSP": component("thermospray inlet"),
}


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
        if word in "PC":
            if word == 'P':
                values['peak_mode'] = 'profile'
            else:
                values['peak_mode'] = 'centroid'
            word = words[i]
            i += 1
        ionization_info = ionization_pat.search(word)
        if ionization_info is not None:
            values['ionization'] = ionization_info.group(0)
            word = words[i]
            i += 1

        return values
    except IndexError:
        return values


_id_template = "controllerType=0 controllerNumber=1 scan="


class _RawFileMetadataLoader(object):
    def _build_scan_type_index(self):
        self.make_iterator(grouped=False)
        index = defaultdict(int)
        analyzer_counter = 1
        analyzer_confs = dict()
        for scan in self:
            index[scan.ms_level] += 1
            fline = self._filter_line(scan._data)
            analyzer = analyzer_map[fline.data['analyzer']]
            try:
                analyzer_confs[analyzer]
            except KeyError:
                analyzer_confs[analyzer] = analyzer_counter
                analyzer_counter += 1
        self.reset()
        self._scan_type_index = index
        self._analyzer_to_configuration_index = analyzer_confs

    def _get_instrument_info(self):
        scan = self.get_scan_by_index(0)
        filter_line = self._filter_line(scan._data)
        ionization_label = filter_line.data.get("ionization")
        try:
            ionization = ionization_map[ionization_label]
        except KeyError:
            ionization = ionization_map['ESI']
        try:
            inlet = inlet_map[ionization_label]
        except KeyError:
            inlet = None

        source_group = ComponentGroup("source", [], 1)
        source_group.add(ionization)
        if inlet is not None:
            source_group.add(inlet)
        configs = []
        for analyzer, counter in sorted(self._analyzer_to_configuration_index.items(), key=lambda x: x[1]):
            analyzer_group = ComponentGroup('analyzer', [analyzer], 2)
            configs.append(InstrumentInformation(counter, [source_group, analyzer_group]))
        self._instrument_config = {
            c.id: c for c in configs
        }
        return configs

    def instrument_configuration(self):
        return sorted(self._instrument_config.values(), key=lambda x: x.id)

    def file_description(self):
        fi = FileInformation({}, [])
        fi.add_file(self.source_file)
        if 1 in self._scan_type_index:
            fi.add_content("MS1 spectrum")
        scan_types = sorted(self._scan_type_index, reverse=True)
        if scan_types:
            if scan_types[0] > 1:
                fi.add_content("MSn spectrum")
        return fi


class _InstrumentMethod(object):
    def __init__(self, method_text):
        self.text = method_text
        (self.isolation_width_by_segment_and_event,
         self.isolation_width_by_segment_and_ms_level) = method_parser(self.text)

    def isolation_width_for(self, segment, event=None, ms_level=None):
        if event is not None:
            try:
                width = self.isolation_width_by_segment_and_event[segment][event]
                return width
            except KeyError:
                return 0.0
        elif ms_level is not None:
            try:
                width = self.isolation_width_by_segment_and_ms_level[segment][ms_level]
                return width
            except KeyError:
                return 0.0
        else:
            raise ValueError("One of event or ms_level must not be None!")


def method_parser(method_text):
    scan_segment_re = re.compile(r"\s*Segment (\d+) Information\s*")
    scan_event_re = re.compile(r"\s*(\d+):.*")
    scan_event_isolation_width_re = re.compile(r"\s*Isolation Width:\s*(\S+)\s*")
    scan_event_iso_w_re = re.compile(r"\s*MS.*:.*\s+IsoW\s+(\S+)\s*")
    repeated_event_re = re.compile(r"\s*Scan Event (\d+) repeated for top (\d+)\s*")
    default_isolation_width_re = re.compile(r"\s*MS(\d+) Isolation Width:\s*(\S+)\s*")

    scan_segment = 1
    scan_event = 0
    scan_event_details = False
    data_dependent_settings = False

    isolation_width_by_segment_and_event = defaultdict(dict)
    isolation_width_by_segment_and_ms_level = defaultdict(dict)

    for line in method_text.splitlines():
        match = scan_segment_re.match(line)

        if match:
            scan_segment = int(match.group(1))
            continue

        if "Scan Event Details" in line:
            scan_event_details = True
            continue

        if scan_event_details:
            match = scan_event_re.match(line)
            if match:
                scan_event = int(match.group(1))
                continue

            match = scan_event_isolation_width_re.match(line)
            if match:
                isolation_width_by_segment_and_event[scan_segment][scan_event] = float(match.group(1))
                continue

            match = scan_event_iso_w_re.match(line)
            if match:
                isolation_width_by_segment_and_event[scan_segment][scan_event] = float(match.group(1))
                continue

            match = repeated_event_re.match(line)
            if match:
                repeated_event = int(match.group(1))
                repeat_count = int(match.group(2))
                repeated_width = isolation_width_by_segment_and_event[scan_segment][repeated_event]
                for i in range(repeated_width + 1, repeat_count + repeated_width):
                    isolation_width_by_segment_and_event[scan_segment][i] = repeated_width
                continue

            if not line.strip():
                scan_event_details = False

        if "Data Dependent Settings" in line:
            data_dependent_settings = True
            continue

        if data_dependent_settings:
            match = default_isolation_width_re.match(line)
            if match:
                ms_level = int(match.group(1))
                width = float(match.group(2))
                isolation_width_by_segment_and_ms_level[scan_segment][ms_level] = width
                continue

            if not line.strip():
                data_dependent_settings = False

    return isolation_width_by_segment_and_event, isolation_width_by_segment_and_ms_level


class ThermoRawScanPtr(Base):
    def __init__(self, scan_number):
        self.scan_number = scan_number
        self.filter_line = None


def _make_id(scan_number):
    try:
        return "%s%d" % (_id_template, (scan_number))
    except TypeError:
        return None


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
        return "%s %r" % (self._scan_id(scan), self._filter_line(scan))

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
            while last_index >= 0 and i < 100:
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
            source=self,
            product_scan_id=_make_id(scan.scan_number))
        return pinfo

    def _activation(self, scan):
        filter_line = self._filter_line(scan)
        activation_type = filter_line.get("activation_type")
        if activation_type is not None:
            energy = filter_line.get("activation_energy")
            return ActivationInformation(activation_type, energy)
        return None

    def _get_scan_segment(self, scan):
        trailer = self._source.GetTrailerExtraForScanNum(scan.scan_number)
        try:
            return int(trailer['Scan Segment'])
        except KeyError:
            return 1

    def _get_scan_event(self, scan):
        trailer = self._source.GetTrailerExtraForScanNum(scan.scan_number)
        try:
            return int(trailer['Scan Event'])
        except KeyError:
            return 1

    def _isolation_window(self, scan):
        ms_level = self._ms_level(scan)
        if ms_level == 1:
            return None
        isolation_width = 0
        trailer = self._source.GetTrailerExtraForScanNum(scan.scan_number)
        try:
            isolation_width = trailer['MS%d Isolation Width' % ms_level]
        except KeyError:
            segment = self._get_scan_segment(scan)
            event = self._get_scan_event(scan)
            isolation_width = self._method.isolation_width_for(segment, event=event)
            if not isolation_width:
                isolation_width = self._method.isolation_width_for(segment, ms_level=ms_level)
        if isolation_width == 0:
            return None
        isolation_width /= 2.
        isolation_mz = self._source.GetPrecursorMassForScanNum(scan.scan_number, ms_level)
        return IsolationWindow(isolation_width, isolation_mz, isolation_width)

    def _instrument_configuration(self, scan):
        fline = self._filter_line(scan)
        try:
            confid = self._analyzer_to_configuration_index[analyzer_map[fline.data.get("analyzer")]]
            return self._instrument_config[confid]
        except KeyError:
            return None


class ThermoRawLoader(ThermoRawDataInterface, RandomAccessScanSource, ScanIterator, _RawFileMetadataLoader):
    def __init__(self, source_file, **kwargs):
        self.source_file = source_file
        self._source = _ThermoRawFileAPI(self.source_file)
        self._producer = None
        self.make_iterator()
        self._scan_cache = WeakValueDictionary()
        self._index = self._pack_index()
        self._first_scan_time = self.get_scan_by_index(0).scan_time
        self._last_scan_time = self.get_scan_by_id(self._source.LastSpectrumNumber).scan_time
        self._method = _InstrumentMethod(self._source.GetInstMethod())
        self._build_scan_type_index()
        self._get_instrument_info()

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

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        if scan_id is not None:
            scan_number = int(str(scan_id).replace(_id_template, '')) - 1
        elif index is not None:
            scan_number = int(index)
        elif rt is not None:
            scan_number = self._source.ScanNumFromRT(rt)
        if require_ms1:
            start_index = scan_number
            while start_index != 0:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            scan_number = start_index
        iterator = self._make_scan_index_producer(start_index=scan_number)
        if grouped:
            self._producer = self._scan_group_iterator(iterator)
        else:
            self._producer = self._single_scan_iterator(iterator)
        return self

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
            if packed.ms_level > 1:
                # inceasing ms level
                if current_level < packed.ms_level:
                    current_level = packed.ms_level
                # decreasing ms level
                elif current_level > packed.ms_level:
                    current_level = packed.ms_level.ms_level
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
