import os
import glob
import re

from collections import defaultdict, OrderedDict

import numpy as np

import ms_deisotope
from ms_deisotope.utils import Base
from ms_deisotope.data_source import (
    Scan, ActivationInformation, PrecursorInformation, ChargeNotProvided, IsolationWindow)
from ms_deisotope.data_source.metadata.scan_traits import (
    ScanAcquisitionInformation, ScanWindow, ScanEventInformation,
    ion_mobility_drift_time)
from ms_deisotope.data_source.metadata.activation import HCD

from . import (MassLynxRawInfoReader,
               MassLynxRawScanReader,
               MassLynxRawChromatogramReader,
               MassLynxRawDefs, MassLynxParameters)


waters_id_pattern = re.compile("function=(\d+) process=(\d+) scan=(\d+)")


class IndexEntry(Base):
    __slots__ = ('function', 'process', 'block', 'scan', 'index', 'id')

    def __init__(self, function=None, process=None, block=None, scan=None, index=None, id=None):
        self.function = function
        self.process = process
        self.block = block
        self.scan = scan
        self.index = index
        self.id = id


class MassLynxReader(object):
    def __init__(self, raw_path):
        self.raw_path = raw_path
        self.info_reader = MassLynxRawInfoReader.MassLynxRawInfoReader(
            raw_path)
        self.chrom_reader = MassLynxRawChromatogramReader.MassLynxRawChromatogramReader(
            raw_path)
        self.scan_reader = MassLynxRawScanReader.MassLynxRawScanReader(
            raw_path)
        self._build_function_index()
        self._build_scan_index()

    # Vendor data access support methods
    def _read_header_properties(self):
        path = os.path.join(self.raw_path, "_HEADER.TXT")
        self.header_properties = {}
        with open(path, 'rt') as fh:
            for line in fh:
                if not line.startswith("$$ "):
                    continue
                name, value = line[3:].split(": ")
                self.header_properties[name] = value.strip()

    def _build_function_index(self):
        num_spectra = 0
        function_file_path_by_number = {}
        function_index_list = []
        for dat in glob.glob(os.path.join(self.raw_path, "_FUNC*.DAT")):
            fnum = int(os.path.basename(dat).replace(
                "_func", '').replace(".dat", ""))
            function_file_path_by_number[fnum - 1] = dat
            function_index_list.append(fnum - 1)
            num_spectra += self.info_reader.GetScansInFunction(fnum - 1)
        function_index_list.sort()
        self.num_spectra = num_spectra
        self.function_index_list = function_index_list
        self.function_file_path_by_number = function_file_path_by_number
        self.ion_mobility_by_function_index = OrderedDict()
        self.times_by_function_index = OrderedDict()
        self.tic_by_function_index = OrderedDict()
        self.sonar_enabled_by_function_index = OrderedDict()
        self.tic_by_function_index = OrderedDict()
        self.times_by_function_index = OrderedDict()
        for fnum, dat in function_file_path_by_number.items():
            cdt = dat.replace(".dat", ".cdt")
            self.ion_mobility_by_function_index[fnum] = False
            if os.path.exists(cdt) and self.info_reader.GetDriftScanCount(fnum) > 0:
                self.ion_mobility_by_function_index[fnum] = True
                self.sonar_enabled_by_function_index[fnum] = self.info_reader.GetScanItem(
                    fnum, 0, MassLynxRawDefs.MassLynxScanItem.SONAR_ENABLED.value)
            time, tic = self.chrom_reader.ReadTIC(fnum)
            self.tic_by_function_index[fnum] = np.array(tic)
            self.times_by_function_index[fnum] = np.array(time)
        self._read_header_properties()

    def _translate_function_type(self, fnum):
        fcode = self.info_reader.GetFunctionType(fnum)
        name = self.info_reader.GetFunctionTypeString(fcode)
        if name in ("MSMS", "TOFD", "MS2", ):
            return 2, "MSn spectrum"
        elif name in ("MS", "TOF MS", ):
            return 1, "MS1 spectrum"
        else:
            raise ValueError("Unknown function type")

    def _build_scan_index(self):
        function_and_scan_by_rt = []
        self.scan_time_to_function_block_map = defaultdict(list)
        num_scans_in_block = 0
        for fnum in self.function_index_list:
            _ms_level, spectrum_type = self._translate_function_type(fnum)
            if spectrum_type not in ("MSn spectrum", "MS1 spectrum"):
                continue
            scan_count = self.info_reader.GetScansInFunction(fnum)
            if self.ion_mobility_by_function_index[fnum]:
                num_scans_in_block = self.info_reader.GetDriftScanCount(fnum)
                for i in range(scan_count):
                    key = self.info_reader.GetRetentionTime(fnum, i)
                    self.scan_time_to_function_block_map[key *
                                                         60].append((fnum, i))
                    function_and_scan_by_rt.append((key, (fnum, i)))
            else:
                for i in range(scan_count):
                    function_and_scan_by_rt.append(
                        (self.info_reader.GetRetentionTime(fnum, i), (fnum, i)))
        function_and_scan_by_rt.sort(key=lambda x: x[0])
        self.index = []
        self.function_blocks = defaultdict(list)
        for rt, (fnum, i) in function_and_scan_by_rt:
            if self.ion_mobility_by_function_index[fnum]:
                # num_scans_in_block = self.info_reader.GetDriftScanCount(fnum)
                block_start = len(self.index)
                for j in range(num_scans_in_block):
                    self.index.append(IndexEntry())
                    ie = self.index[-1]
                    ie.function = fnum
                    ie.process = 0
                    ie.block = i
                    ie.scan = j
                    ie.index = len(self.index) - 1
                    ie.id = "function=%d process=%d scan=%d" % (ie.function + 1, ie.process,
                                                                num_scans_in_block * ie.block + ie.scan + 1)
                block_end = len(self.index)
                self.function_blocks[fnum].append((block_start, block_end))
            else:
                block_start = len(self.index)
                self.index.append(IndexEntry())
                ie = self.index[-1]
                ie.function = fnum
                ie.process = 0
                ie.block = -1
                ie.scan = i
                ie.index = len(self.index) - 1
                ie.id = "function=%d process=%d scan=%d" % (ie.function + 1, ie.process,
                                                            ie.scan + 1)
                block_end = len(self.index)
                self.function_blocks[fnum].append((block_start, block_end))

    # ScanSource methods

    def _make_scan(self, data):
        return Scan(data, self)

    def _scan_time(self, data):
        if data.block >= 0:
            return self.info_reader.GetRetentionTime(data.function, data.block)
        else:
            return self.info_reader.GetRetentionTime(data.function, data.scan)

    def _drift_time(self, data):
        if data.block >= 0:
            return self.info_reader.GetDriftTime(data.function, data.scan)
        else:
            return None

    def _scan_id(self, data):
        return data.id

    def _scan_index(self, data):
        return data.index

    def _ms_level(self, data):
        ms_level, _scan_type = self._translate_function_type(data.function)
        return ms_level

    def _scan_arrays(self, data):
        if data.block >= 0:
            mz, inten = self.scan_reader.ReadDriftScan(
                data.function, data.block, data.scan)
        else:
            mz, inten = self.scan_reader.ReadScan(data.function, data.scan)
        return np.array(mz), np.array(inten)

    def _precursor_information(self, data):
        if self._ms_level(data) == 1:
            return None
        set_mass_str = self.info_reader.GetScanItem(
            data.function,
            data.block if data.block >= 0 else data.scan,
            MassLynxRawDefs.MassLynxScanItem.SET_MASS.value
        )
        if set_mass_str:
            set_mass = float(set_mass_str)
        else:
            set_mass = 0.0
        if set_mass == 0:
            lower_bound, upper_bound = self.info_reader.GetAcquisitionMassRange(
                data.function)
            set_mass = (lower_bound + upper_bound) / 2.
        pinfo = PrecursorInformation(
            set_mass, 0, ChargeNotProvided, source=self, product_scan_id=data.id)
        return pinfo

    def _isolation_window(self, data):
        if self._ms_level(data) == 1:
            return None
        set_mass_str = self.info_reader.GetScanItem(
            data.function,
            data.block if data.block >= 0 else data.scan,
            MassLynxRawDefs.MassLynxScanItem.SET_MASS.value
        )
        if set_mass_str:
            set_mass = float(set_mass_str)
        else:
            set_mass = 0.0
        if set_mass == 0:
            lower_bound, upper_bound = self.info_reader.GetAcquisitionMassRange(
                data.function)
            set_mass = (lower_bound + upper_bound) / 2.
            lower_bound_offset = upper_bound_offset = upper_bound - set_mass
        else:
            lower_bound_offset = upper_bound_offset = 0
        return IsolationWindow(
            lower_bound_offset, set_mass, upper_bound_offset)

    def _is_profile(self, data):
        return self.info_reader.IsContinuum(data.function)

    def _acquisition_information(self, data):
        scan_window = ScanWindow(
            *self.info_reader.GetAcquisitionMassRange(data.function))
        scan_time = self._scan_time(data)
        drift_time = self._drift_time(data)
        event = ScanEventInformation(scan_time, [scan_window], traits={
            'preset scan configuration': data.function + 1,
        })
        if drift_time is not None:
            event.ion_mobility.add_ion_mobility(
                ion_mobility_drift_time, drift_time)
        return ScanAcquisitionInformation('no combination', [event])

    def _polarity(self, data):
        mode = self.info_reader.GetIonMode(data.function)
        s = self.info_reader.GetIonModeString(mode)
        if s.endswith('+'):
            return 1
        elif s.endswith('-'):
            return -1
        raise ValueError("Unknown Ion Mode %r" % (s, ))
        # return 1

    def _activation(self, data):
        energy_str = self.info_reader.GetScanItem(
            data.function, data.block if data.block >= 0 else data.scan, MassLynxRawDefs.MassLynxScanItem.COLLISION_ENERGY)
        if energy_str:
            energy = float(energy_str)
            return ActivationInformation(HCD, energy)

    def get_scan_by_index(self, index):
        return self._make_scan(self.index[index])

    def get_scan_by_id(self, scan_id):
        match = waters_id_pattern.search(scan_id)
        if not match:
            raise KeyError(scan_id)
        fnum1, proc, scan = tuple(map(int, match.groups()))
        fnum = fnum1 - 1
        intervals = self.function_blocks[fnum]
        for interval in intervals:
            ie = self.index[interval[1] - 1]
            snum = int(waters_id_pattern.search(ie.id).group(3))
            if scan <= snum:
                ie = self.index[interval[0]]
                snum = int(waters_id_pattern.search(ie.id).group(3))
                if snum <= scan:
                    for i in range(interval[0], interval[1]):
                        ie = self.index[i]
                        if ie.id == scan_id:
                            return self._make_scan(ie)
                    else:
                        raise KeyError(scan_id)
        raise KeyError(scan_id)