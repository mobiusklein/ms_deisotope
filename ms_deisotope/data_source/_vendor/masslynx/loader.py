import os
import glob
import re
import sys

from collections import defaultdict, OrderedDict

import numpy as np

from ms_deisotope.utils import Base
from ms_deisotope.data_source.common import (
    Scan, ActivationInformation, PrecursorInformation,
    ChargeNotProvided, IsolationWindow)
from ms_deisotope.data_source.scan.loader import ScanDataSource, RandomAccessScanSource
from ms_deisotope.data_source.scan.mobility_frame import IonMobilityFrame, IonMobilitySourceRandomAccessFrameSource
from ms_deisotope.data_source.metadata.scan_traits import (
    ScanAcquisitionInformation, ScanWindow, ScanEventInformation,
    ion_mobility_drift_time)
from ms_deisotope.data_source.metadata.activation import HCD

from . import (MassLynxRawInfoReader,
               MassLynxRawScanReader,
               MassLynxRawChromatogramReader,
               MassLynxRawDefs,
               MassLynxParameters,
               MassLynxLockMassProcessor,
               libload
              )



waters_id_pattern = re.compile(r"function=(\d+) process=(\d+) scan=(\d+)")
waters_frame_id_pattern = re.compile(r"function=(\d+) process=(\d+) startScan=(\d+) endScan=(\d+)")


def is_waters_raw_dir(path):
    '''Detect whether or not the file referenced by ``path``
    is a Waters RAW directory.

    Parameters
    ----------
    path: :class:`str`
        The path to test

    Returns
    -------
    :class:`bool`:
        Whether or not the file is a Waters RAW directory.
    '''
    if not libload.proxy._is_loaded():
        try:
            libload.register_dll()
        except ImportError:
            return False
    try:
        MassLynxRawLoader(path)
        return True
    except Exception:   # pylint: disable=broad-except
        return False


def infer_reader(path):
    '''If the file referenced by ``path`` is a Waters RAW
    directory, return the callable (:class:`MassLynxRawLoader`) to
    open it, otherwise raise an exception.

    Parameters
    ----------
    path: :class:`str`
        The path to test

    Returns
    -------
    :class:`type`:
        The type to use to open the file

    Raises
    ------
    :class:`ValueError`:
        If the file is not a Waters RAW file
    '''
    if is_waters_raw_dir(path):
        return MassLynxRawLoader
    raise ValueError("Not Waters Raw Directory")


determine_if_available = libload.determine_if_available


class Cycle(Base):
    __slots__ = ("function", "block", "start_scan", "end_scan", "index", "id")

    def __init__(self, function=None, block=None, start_scan=None, end_scan=None, index=None, id=None):
        self.function = function
        self.block = block
        self.start_scan = start_scan
        self.end_scan = end_scan
        self.index = index
        self.id = id

    def __str__(self):
        return self.id


class IndexEntry(Base):
    __slots__ = ('function', 'process', 'block', 'scan', 'index', 'id')

    def __init__(self, function=None, process=None, block=None, scan=None, index=None, id=None):
        self.function = function
        self.process = process
        self.block = block
        self.scan = scan
        self.index = index
        self.id = id

    def __str__(self):
        return self.id


class WatersMSECycleSourceMixin(IonMobilitySourceRandomAccessFrameSource):
    def _frame_id(self, data):
        return data.id

    def _frame_index(self, data):
        return data.index

    def _frame_time(self, data):
        scan = self.get_scan_by_index(data.start_scan)
        return scan.scan_time

    def _frame_ms_level(self, data):
        scan = self.get_scan_by_index(data.start_scan)
        return scan.ms_level

    def _frame_start_scan_index(self, data):
        return data.start_scan

    def _frame_end_scan_index(self, data):
        return data.end_scan

    def _frame_precursor_information(self, data):
        scan = self.get_scan_by_index(data.start_scan)
        return scan.precursor_information

    def _frame_activation(self, data):
        scan = self.get_scan_by_index(data.start_scan)
        return scan.activation

    def _frame_isolation_window(self, data):
        scan = self.get_scan_by_index(data.start_scan)
        return scan.isolation_window

    def _frame_polarity(self, data):
        scan = self.get_scan_by_index(data.start_scan)
        return scan.polarity

    def get_frame_by_index(self, index):
        cycle = self.cycle_index[index]
        return self._make_frame(cycle)

    def get_frame_by_time(self, time):
        lo = 0
        hi = len(self.cycle_index)

        best_match = None
        best_error = float('inf')

        if time == float('inf'):
            return self.get_frame_by_index(self.cycle_index[-1].index)

        while hi != lo:
            mid = (hi + lo) // 2
            sid = self.cycle_index[mid]
            frame = self.get_frame_by_index(sid.index)
            if not self._validate_frame(frame):
                sid = self.cycle_index[mid - 1]
                frame = self.get_frame_by_index(sid.index)
                if not self._validate_frame(frame):
                    sid = self.cycle_index[mid - 2]
                    frame = self.get_frame_by_index(sid.index)

            frame_time = frame.time
            err = abs(frame_time - time)
            if err < best_error:
                best_error = err
                best_match = frame
            if frame_time == time:
                # Rewind here
                i = frame.index - 1
                while i >= 0:
                    prev_frame = self.get_frame_by_index(i)
                    if prev_frame.time == frame.time:
                        frame = prev_frame
                        i -= 1
                    else:
                        break
                return frame
            elif (hi - lo) == 1:
                # Rewind here
                frame = best_match
                i = frame.index - 1
                while i >= 0:
                    prev_frame = self.get_frame_by_index(i)
                    if prev_frame.time == frame.time:
                        frame = prev_frame
                        i -= 1
                    else:
                        break
                return frame
            elif frame_time > time:
                hi = mid
            else:
                lo = mid

    def _validate_frame(self, data):
        return True

    def _make_frame(self, data):
        return IonMobilityFrame(data, self)

    def _cache_frame(self, frame):
        pass

    def _default_frame_iterator(self, start_index=None):
        if start_index is None:
            start_index = 0
        for i in range(start_index, len(self.cycle_index)):
            yield self.cycle_index[i]

    def make_frame_iterator(self, iterator=None, grouped=False):
        from ms_deisotope.data_source.scan.scan_iterator import (
            _SingleScanIteratorImpl, _GroupedScanIteratorImpl, MSEIterator)
        if iterator is None:
            iterator = self._default_frame_iterator()

        if grouped == 'mse':
            strategy = MSEIterator(
                iterator, self._make_frame, self.low_energy_function, self.lockmass_function,
                self._validate_frame, self._cache_frame)
        elif grouped:
            strategy = _GroupedScanIteratorImpl(
                iterator, self._make_frame, self._validate_frame, self._cache_frame)
        else:
            strategy = _SingleScanIteratorImpl(
                iterator, self._make_frame, self._validate_frame, self._cache_frame)
        return strategy

    def start_from_frame(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        if scan_id is not None:
            frame = self.get_frame_by_id(scan_id)
            start_index = frame.index
        elif index is not None:
            start_index = index
        elif rt is not None:
            frame = self.get_frame_by_time(rt)
            start_index = frame.index
        if require_ms1 and grouped != 'mse':
            while start_index != 0:
                frame = self.get_frame_by_index(start_index)
                if frame.ms_level > 1:
                    start_index -= 1
                else:
                    break
        self._producer = self.make_frame_iterator(
            self._default_frame_iterator(start_index), grouped=grouped)
        return self


class WatersMassLynxScanSource(ScanDataSource):
    def _make_scan(self, data):
        return Scan(data, self)

    def _scan_time(self, scan):
        if scan.block >= 0:
            return self.info_reader.GetRetentionTime(scan.function, scan.block)
        else:
            return self.info_reader.GetRetentionTime(scan.function, scan.scan)

    def _drift_time(self, scan):
        if scan.block >= 0:
            return self.info_reader.GetDriftTime(scan.function, scan.scan)
        else:
            return None

    def _scan_id(self, scan):
        return scan.id

    def _scan_index(self, scan):
        return scan.index

    def _ms_level(self, scan):
        ms_level, _scan_type = self._translate_function_type(scan.function)
        return ms_level

    def _scan_arrays(self, scan):
        if scan.block >= 0:
            mz, inten = self.scan_reader.ReadDriftScan(
                scan.function, scan.block, scan.scan)
        else:
            mz, inten = self.scan_reader.ReadScan(scan.function, scan.scan)
        return np.array(mz), np.array(inten)

    def _precursor_information(self, scan):
        if self._ms_level(scan) == 1:
            return None
        set_mass_str = self.info_reader.GetScanItem(
            scan.function,
            scan.block if scan.block >= 0 else scan.scan,
            MassLynxRawDefs.MassLynxScanItem.SET_MASS.value
        )
        if set_mass_str:
            set_mass = float(set_mass_str)
        else:
            set_mass = 0.0
        if set_mass == 0:
            lower_bound, upper_bound = self.info_reader.GetAcquisitionMassRange(
                scan.function)
            set_mass = (lower_bound + upper_bound) / 2.
        pinfo = PrecursorInformation(
            set_mass, 0, ChargeNotProvided, source=self, product_scan_id=scan.id)
        return pinfo

    def _isolation_window(self, scan):
        if self._ms_level(scan) == 1:
            return None
        set_mass_str = self.info_reader.GetScanItem(
            scan.function,
            scan.block if scan.block >= 0 else scan.scan,
            MassLynxRawDefs.MassLynxScanItem.SET_MASS.value
        )
        if set_mass_str:
            set_mass = float(set_mass_str)
        else:
            set_mass = 0.0
        if set_mass == 0:
            lower_bound, upper_bound = self.info_reader.GetAcquisitionMassRange(
                scan.function)
            set_mass = (lower_bound + upper_bound) / 2.
            lower_bound_offset = upper_bound_offset = upper_bound - set_mass
        else:
            lower_bound_offset = upper_bound_offset = 0
        return IsolationWindow(
            lower_bound_offset, set_mass, upper_bound_offset)

    def _is_profile(self, scan):
        return self.info_reader.IsContinuum(scan.function)

    def _scan_title(self, scan):
        return self._scan_id(scan)

    def _acquisition_information(self, scan):
        scan_window = ScanWindow(
            *self.info_reader.GetAcquisitionMassRange(scan.function))
        scan_time = self._scan_time(scan)
        drift_time = self._drift_time(scan)
        event = ScanEventInformation(scan_time, [scan_window], traits={
            'preset scan configuration': scan.function + 1,
        })
        if drift_time is not None:
            event._ion_mobility.add_ion_mobility(
                ion_mobility_drift_time, drift_time)
        return ScanAcquisitionInformation('no combination', [event])

    def _polarity(self, scan):
        mode = self.info_reader.GetIonMode(scan.function)
        s = self.info_reader.GetIonModeString(mode)
        if s.endswith('+'):
            return 1
        elif s.endswith('-'):
            return -1
        raise ValueError("Unknown Ion Mode %r" % (s, ))
        # return 1

    def _activation(self, scan):
        energy_str = self.info_reader.GetScanItem(
            scan.function, scan.block if scan.block >= 0 else scan.scan, MassLynxRawDefs.MassLynxScanItem.COLLISION_ENERGY)
        if energy_str:
            energy = float(energy_str)
            return ActivationInformation(HCD, energy)


class MassLynxRawLoader(RandomAccessScanSource, WatersMassLynxScanSource, WatersMSECycleSourceMixin):
    def __init__(self, raw_path, lockmass_config=None):
        if sys.version_info.major == 2:
            if not isinstance(raw_path, str):
                raw_path = str(raw_path)
        self.source_file = raw_path

        self.info_reader = MassLynxRawInfoReader.MassLynxRawInfoReader(
            self.source_file)
        self.chrom_reader = MassLynxRawChromatogramReader.MassLynxRawChromatogramReader(
            self.source_file)
        self.scan_reader = MassLynxRawScanReader.MassLynxRawScanReader(
            self.source_file)
        self.lockmass_processor = MassLynxLockMassProcessor.MassLynxLockMassProcessor()
        self.lockmass_processor.SetRawData(self.scan_reader)
        self.configure_lockmass(lockmass_config)
        self.index = []
        self.cycle_index = []
        self._producer = None
        self.initialize_scan_cache()

        self._build_function_index()
        self._build_scan_index()
        self.lockmass_function = self.function_index_list[-1] + 1
        self.low_energy_function = self.function_index_list[0] + 1

        self._producer = self.make_frame_iterator()

    def configure_lockmass(self, lockmass_config=None):
        if lockmass_config is None:
            return
        if isinstance(lockmass_config, (list, tuple)):
            if len(lockmass_config) == 1:
                mass = float(lockmass_config[0])
                tolerance = 0.25
            else:
                mass = float(lockmass_config[0])
                tolerance = float(lockmass_config[1])
        else:
            mass = float(lockmass_config)
            tolerance = 0.25

        params = MassLynxParameters.MassLynxParameters()
        params.Set(MassLynxRawDefs.LockMassParameter.MASS, mass)
        params.Set(MassLynxRawDefs.LockMassParameter.TOLERANCE, tolerance)
        self.lockmass_processor.SetParameters(params)
        if self.lockmass_processor.CanLockMassCorrect():
            self.lockmass_processor.LockMassCorrect()

    def __repr__(self):
        return "MassLynxRawLoader(%r)" % (self.source_file)

    # Vendor data access support methods
    def _read_header_properties(self):
        path = os.path.join(self.source_file, "_HEADER.TXT")
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
        for dat in glob.glob(os.path.join(self.source_file, "_FUNC*.DAT")):
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
        # self.scan_time_to_function_block_map = defaultdict(list)
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
                    # self.scan_time_to_function_block_map[key *
                    #                                      60].append((fnum, i))
                    function_and_scan_by_rt.append((key, (fnum, i)))
            else:
                for i in range(scan_count):
                    function_and_scan_by_rt.append(
                        (self.info_reader.GetRetentionTime(fnum, i), (fnum, i)))
        function_and_scan_by_rt.sort(key=lambda x: x[0])
        self.index = []
        self.cycle_index = []
        self.function_blocks = defaultdict(list)
        for _rt, (fnum, i) in function_and_scan_by_rt:
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
                cyc = Cycle(
                    fnum, i, block_start, block_end,
                    len(self.cycle_index),
                    id="function=%d process=0 startScan=%d endScan=%d" % (
                        fnum + 1,
                        num_scans_in_block * i + 1,
                        num_scans_in_block * i + num_scans_in_block,
                    ))
                self.cycle_index.append(cyc)
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
                cyc = Cycle(
                    ie.function, ie.block, block_start, block_end,
                    len(self.cycle_index),
                    id="function=%d process=0 startScan=%d endScan=%d" % (
                        fnum + 1,
                        i + 1,
                        i + 2,
                    ))
                self.cycle_index.append(cyc)

    # RandomAccessScanSource methods
    def __len__(self):
        return len(self.index)

    def get_scan_by_index(self, index):
        ie = self.index[index]
        if ie.id in self._scan_cache:
            return self._scan_cache[ie.id]
        scan = self._make_scan(ie)
        self._scan_cache[ie.id] = scan
        return scan

    def get_scan_by_id(self, scan_id):
        match = waters_id_pattern.search(scan_id)
        if not match:
            raise KeyError(scan_id)
        fnum1, _proc, scan = tuple(map(int, match.groups()))
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
                            return self.get_scan_by_index(ie.index)
                    raise KeyError(scan_id)

    def get_scan_by_time(self, time):
        lo = 0
        hi = len(self.index)

        best_match = None
        best_error = float('inf')

        if time == float('inf'):
            return self.get_scan_by_index(self.index[-1].index)

        while hi != lo:
            mid = (hi + lo) // 2
            sid = self.index[mid]
            scan = self.get_scan_by_index(sid.index)
            if not self._validate(scan):
                sid = self.index[mid - 1]
                scan = self.get_scan_by_index(sid.index)
                if not self._validate(scan):
                    sid = self.index[mid - 2]
                    scan = self.get_scan_by_index(sid.index)

            scan_time = scan.scan_time
            err = abs(scan_time - time)
            if err < best_error:
                best_error = err
                best_match = scan
            if scan_time == time:
                # Rewind here
                i = scan.index - 1
                while i >= 0:
                    prev_scan = self.get_scan_by_index(i)
                    if prev_scan.scan_time == scan.scan_time:
                        scan = prev_scan
                        i -= 1
                    else:
                        break
                return scan
            elif (hi - lo) == 1:
                # Rewind here
                scan = best_match
                i = scan.index - 1
                while i >= 0:
                    prev_scan = self.get_scan_by_index(i)
                    if prev_scan.scan_time == scan.scan_time:
                        scan = prev_scan
                        i -= 1
                    else:
                        break
                return scan
            elif scan_time > time:
                hi = mid
            else:
                lo = mid

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        if scan_id is not None:
            scan = self.get_scan_by_id(scan_id)
            start_index = scan.index
        elif index is not None:
            start_index = index
        elif rt is not None:
            scan = self.get_scan_by_time(rt)
            start_index = scan.index
        if require_ms1:
            while start_index != 0:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            current_block = self.index[start_index].block
            while start_index > 0:
                ie = self.index[start_index - 1]
                if ie.block != current_block:
                    break
                else:
                    start_index -= 1
        self._producer = self._make_pointer_iterator(start_index)
        return self

    # ScanIterator methods
    def next(self):
        return next(self._producer)

    def _make_default_iterator(self):
        return self._make_pointer_iterator()

    def _make_scan_index_producer(self, start_index=None, start_time=None):
        if start_index is not None:
            return range(start_index, len(self.index))
        elif start_time is not None:
            raise NotImplementedError()
            # start_index = self._scan_time_to_scan_number(start_time)
            # while start_index != 0:
            #     scan = self.get_scan_by_index(start_index)
            #     if scan.ms_level > 1:
            #         start_index -= 1
            #     else:
            #         break
            # return range(start_index, len(self.index))
        else:
            return range(0, len(self.index))

    def _make_pointer_iterator(self, start_index=None, start_time=None):
        iterator = self._make_scan_index_producer(start_index, start_time)
        for i in iterator:
            yield self._make_scan(self.index[i])

    def reset(self):
        self.initialize_scan_cache()
        self._producer = self.make_frame_iterator()
