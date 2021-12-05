import re
from collections import defaultdict

import numpy as np

from pyteomics.auxiliary import unitfloat

from ms_deisotope.utils import Base

from .metadata.instrument_components import (
    component, ComponentGroup, InstrumentInformation)
from .metadata.file_information import (
    FileInformation, MS_MS1_Spectrum, MS_MSn_Spectrum)

from .common import ScanFileMetadataBase

analyzer_pat = re.compile(r"(?P<mass_analyzer_type>ITMS|TQMS|SQMS|TOFMS|FTMS|SECTOR)")
polarity_pat = re.compile(r"(?P<polarity>[\+\-])")
point_type_pat = re.compile(r"(?P<point_type>[CP])")
ionization_pat = re.compile(r"(?P<ionization_type>EI|CI|FAB|APCI|ESI|APCI|NSI|TSP|FD|MALDI|GD)")
scan_type_pat = re.compile(r"(?P<scan_type>FULL|SIM|SRM|CRM|Z|Q1MS|Q3MS)")
ms_level_pat = re.compile(r" ms(?P<level>\d*) ")
compensation_voltage_pat = re.compile(r"CV=(-?\d+(?:\.\d+))?")

activation_pat = re.compile(
    r"""(?:(?P<isolation_mz>\d+\.\d*)@
        (?P<activation_type>[a-z]+)
        (?P<activation_energy>\d*\.?\d*))""", re.VERBOSE)
activation_mode_pat = re.compile(
    r"""(?P<activation_type>[a-z]+)
        (?P<activation_energy>\d*\.\d*)""", re.VERBOSE)
scan_window_pat = re.compile(
    r"""
    \[(?P<scan_start>[0-9\.]+)-(?P<scan_end>[0-9\.]+)\]
    """, re.VERBOSE)

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


class FilterString(str):
    def __init__(self, value):
        self.data = self._parse()

    def get(self, key):
        return self.data.get(key)

    def _parse(self):
        return filter_string_parser(self)


def filter_string_parser(line):
    """Parses instrument information from Thermo's filter string

    Parameters
    ----------
    line : str
        The filter string associated with a scan

    Returns
    -------
    dict
        Fields extracted from the filter string
    """
    words = line.upper().split(" ")
    values = dict()
    i = 0
    values['supplemental_activation'] = " sa " in line
    ms_level_info = ms_level_pat.search(line)
    if ms_level_info is not None:
        ms_level_data = ms_level_info.groupdict()
        level = ms_level_data.get("level")
        if level != "":
            parts = line[ms_level_info.end():].split(" ")
            tandem_sequence = []
            for part in parts:
                activation_info = activation_pat.search(part)
                if activation_info is not None:
                    activation_info = activation_info.groupdict()
                    activation_event = dict()
                    activation_event["isolation_mz"] = float(activation_info['isolation_mz'])
                    activation_event["activation_type"] = [activation_info['activation_type']]
                    activation_event["activation_energy"] = [float(activation_info['activation_energy'])]
                    if part.count("@") > 1:
                        act_events = activation_mode_pat.finditer(part)
                        # discard the first match which we already recorded
                        next(act_events)
                        for match in act_events:
                            act_type, act_energy = match.groups()
                            act_energy = float(act_energy)
                            activation_event["activation_type"].append(act_type)
                            activation_event['activation_energy'].append(act_energy)
                    tandem_sequence.append(activation_event)
            values['ms_level'] = int(level)
            values['tandem_sequence'] = tandem_sequence

    scan_window_info = scan_window_pat.search(line)
    if scan_window_info is not None:
        values['scan_window'] = (float(scan_window_info.group(1)), float(scan_window_info.group(2)))

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

        cv_info = compensation_voltage_pat.search(word)
        if cv_info is not None:
            values['compensation_voltage'] = unitfloat(cv_info.group(1), None)
            word = words[i]
            i += 1

        return values
    except IndexError:
        return values


_id_template = "controllerType=0 controllerNumber=1 scan="


def _make_id(scan_number):
    try:
        return "%s%d" % (_id_template, (scan_number))
    except TypeError:
        return None


def _parse_id(scan_id):
    return int(scan_id.replace(_id_template, ""))


class _RawFileMetadataLoader(ScanFileMetadataBase):
    def _get_instrument_serial_number(self):
        return ''

    def _get_instrument_model_name(self):
        return ''

    def _build_scan_type_index(self):
        self.make_iterator(grouped=False)
        index = defaultdict(int)
        analyzer_counter = 1
        analyzer_confs = dict()
        previous_ms_levels = dict()
        last_ms_levels = dict()
        for scan in self: # pylint: disable=not-an-iterable
            ms_level = scan.ms_level
            index[ms_level] += 1
            idx = scan.index
            previous_ms_levels[idx] = tuple(
                last_ms_levels.get(j) for j in range(1, ms_level + 1))

            last_ms_levels[ms_level] = idx

            fline = self._filter_string(scan._data)
            analyzer = analyzer_map[fline.data['analyzer']]
            try:
                analyzer_confs[analyzer]
            except KeyError:
                analyzer_confs[analyzer] = analyzer_counter
                analyzer_counter += 1
        self.reset()
        self._scan_type_index = index
        self._analyzer_to_configuration_index = analyzer_confs
        self._previous_ms_levels = previous_ms_levels

    def _get_previous_scan_index_for_ms_level(self, scan_index, ms_level):
        if ms_level < 1:
            raise ValueError("Cannot handle ms_level less than 1")
        levels = self._previous_ms_levels[scan_index]
        level_i = ms_level - 1
        prev_index = levels[level_i]
        return prev_index


    def _get_instrument_info(self):
        scan = self.get_scan_by_index(0)
        filter_string = self._filter_string(scan._data)
        ionization_label = filter_string.data.get("ionization")
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
            # this is the most common Thermo detector, but it is not universal. To get this right,
            # we'd need to reproduce the Proteowizard conversion table @
            # Thermo::(Reader_Thermo_Detail::)createInstrumentConfigurations
            detector_group = ComponentGroup("detector", [component('inductive detector')], 3)
            configs.append(InstrumentInformation(
                counter, [source_group, analyzer_group, detector_group],
                self._get_instrument_model_name(), self._get_instrument_serial_number())
            )
        self._instrument_config = {
            c.id: c for c in configs
        }
        return configs

    def instrument_configuration(self):
        return sorted(self._instrument_config.values(), key=lambda x: x.id)

    def file_description(self):
        fi = FileInformation({}, [])
        fi.add_file(self.source_file)
        sf = fi.source_files[0]
        sf.add_checksum("sha1")
        if 1 in self._scan_type_index:
            fi.add_content(MS_MS1_Spectrum)
        scan_types = sorted(self._scan_type_index, reverse=True)
        if scan_types:
            if scan_types[0] > 1:
                fi.add_content(MS_MSn_Spectrum)
        return fi

    def data_processing(self):
        return []


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
                for i in np.arange(repeated_width + 1, repeat_count + repeated_width):
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
        self.filter_string = None
        self.trailer_values = None

    def validate(self, source):
        try:
            source._scan_time(self)
            return True
        except IOError:
            return False
