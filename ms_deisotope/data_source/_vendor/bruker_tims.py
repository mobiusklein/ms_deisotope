'''This interface to the Bruker TDF file format is based upon the TimsData library
that is part of ProteoWizard (Apache II), and a common kernel of Python ctypes functions
found in dia-pasef and other Python wrappers of libtimsdata.
'''
import os
import re
import sqlite3
import sys
import warnings

from ctypes import cdll, c_char_p, c_uint32, c_int32, c_uint64, c_int64, c_double, c_void_p, POINTER, create_string_buffer
from weakref import WeakValueDictionary

import numpy as np

from ms_peak_picker import reprofile, pick_peaks, PeakIndex, PeakSet

from ms_deisotope.utils import Base
from ms_deisotope.data_source.metadata import software
from ms_deisotope.data_source.scan.loader import ScanDataSource, RandomAccessScanSource
from ms_deisotope.data_source.scan import PrecursorInformation, RawDataArrays
from ms_deisotope.data_source.scan.mobility_frame import IonMobilityFrame, IonMobilitySourceRandomAccessFrameSource

from ms_deisotope.data_source.metadata.activation import ActivationInformation, dissociation_methods
from ms_deisotope.data_source.metadata.scan_traits import (
    IsolationWindow, inverse_reduced_ion_mobility,
    ScanAcquisitionInformation, ScanWindow, ScanEventInformation)
from ms_deisotope.peak_dependency_network import Interval
from ms_deisotope.averagine import neutral_mass as calculate_neutral_mass

if sys.platform[:5] == "win32" or sys.platform[:5] == "win64":
    libname = "timsdata.dll"
elif sys.platform[:5] == "linux":
    libname = "libtimsdata.so"
else:
    raise Exception("Unsupported platform.")


dll = None


def load_library(search_paths=None):
    if search_paths is None:
        search_paths = []
    elif isinstance(search_paths, str):
        search_paths = [search_paths]
    global dll
    if dll is None:
        for lib_path in search_paths:
            try:
                dll = _load_library(lib_path)
            except (Exception) as err:
                continue
            if dll is not None:
                break
    return dll


def _load_library(lib_path):
    dll = cdll.LoadLibrary(os.path.realpath(lib_path))
    dll.tims_open.argtypes = [c_char_p, c_uint32]
    dll.tims_open.restype = c_uint64

    dll.tims_close.argtypes = [c_uint64]
    dll.tims_close.restype = None

    dll.tims_get_last_error_string.argtypes = [c_char_p, c_uint32]
    dll.tims_get_last_error_string.restype = c_uint32

    dll.tims_has_recalibrated_state.argtypes = [c_uint64]
    dll.tims_has_recalibrated_state.restype = c_uint32

    dll.tims_read_scans_v2.argtypes = [
        c_uint64, c_int64, c_uint32, c_uint32, c_void_p, c_uint32]
    dll.tims_read_scans_v2.restype = c_uint32

    # dll.tims_oneoverk0_to_ccs_for_mz.argtypes = [c_double, c_int32, c_double]
    # dll.tims_oneoverk0_to_ccs_for_mz.restype = c_double

    # dll.tims_ccs_to_oneoverk0_for_mz.argtypes = [c_double, c_int32, c_double]
    # dll.tims_ccs_to_oneoverk0_for_mz.restype = c_double


    convfunc_argtypes = [c_uint64, c_int64, POINTER(
        c_double), POINTER(c_double), c_uint32]

    for fn in [dll.tims_index_to_mz,
               dll.tims_mz_to_index,
               dll.tims_scannum_to_oneoverk0,
               dll.tims_oneoverk0_to_scannum,
               dll.tims_scannum_to_voltage,
               dll.tims_voltage_to_scannum]:
        fn.argtypes = convfunc_argtypes
        fn.restype = c_uint32
    return dll


def throw_tims_error(dll_handle):
    """Raise the last thronw timsdata error as a
    :exc:`RuntimeError`
    """

    size = dll_handle.tims_get_last_error_string(None, 0)
    buff = create_string_buffer(size)
    dll_handle.tims_get_last_error_string(buff, size)
    raise RuntimeError(buff.value)


def msms_type_to_ms_level(enum):
    if enum == 0:
        return 1
    else:
        return 2


msms_type_to_label = {
    0: "MS1 Scan",
    2: "MS2 Scan",
    8: "PASEF MS2 Scan",
    9: "DIA-PASEF",
}


msms_type_to_metadata_table = {
    2: "FrameMsMsInfo",
    8: "PasefFrameMsMsInfo",
    9: "DiaFrameMsMsWindows"
}


class TIMSMetadata(object):
    def _read_global_metadata(self):
        self._acquisition_parameters = acquisition_parameters = {}
        self._instrument_configuration = instrument_configuration = {}
        self._software = software_map = {
            "acquisition_software": {},
            "control_software": {},
        }
        q = self.connection.execute("SELECT Key, Value FROM GlobalMetadata;")
        for key, value in q:
            if key == "AcquistionSoftware":
                software_map['acquisition_software']['name'] = value
            elif key == "AcquisitionSoftwareVersion":
                software_map['acquisition_software']['version'] = value
            elif key == "InstrmentFamily":
                instrument_configuration['instrument_family'] = value
            elif key == "InstrmentRevision":
                instrument_configuration['instrument_revision'] = value
            elif key == "InstrumentSerialNumber":
                instrument_configuration['serial_number'] = value
            elif key == "AcquisitionDateTime":
                acquisition_parameters['acquisition_date'] = value
            elif key == "OperatorName":
                acquisition_parameters['operator_name'] = value
            elif key == "MzAcqRangeLower":
                acquisition_parameters['scan_window_lower'] = float(value)
            elif key == "MzAcqRangeUpper":
                acquisition_parameters['scan_window_upper'] = float(value)

    def _build_frame_index(self):
        self._frame_counts = {}
        q = self.connection.execute(
            "SELECT MsMsType, Count(*) FROM Frames GROUP BY MsMsType;")
        total = 0
        for scan_type_enum, count in q:
            self._frame_counts[msms_type_to_label[scan_type_enum]] = count
            total += count
        self._frame_counts['Total'] = total
        frame_id_scan_count = self.connection.execute("SELECT Id, NumScans, Time FROM Frames ORDER BY Id;")
        count_to_id = []
        time_to_id = []
        total = 0
        for f_c in frame_id_scan_count:
            count_to_id.append([total, f_c['Id']])
            total += f_c['NumScans']
            time_to_id.append([f_c['Time'], f_c['Id']])
        self._total_scans = total
        self._scan_count_to_frame_id = np.array(count_to_id)
        self._scan_time_to_frame_id = np.array(time_to_id)

    def _read_metadata(self):
        self._read_global_metadata()
        self._build_frame_index()

    def frame_count(self):
        return self._frame_counts['Total']

    def scan_count(self):
        return self._total_scans


class TIMSFrameInformation(Base):
    def __init__(self, source, id, accumulation_time, max_intensity, msms_type, mz_calibration, num_peaks, num_scans,
                 polarity, property_group, ramp_time, scan_mode, summed_intensities, t1, t2, time, tims_calibration,
                 tims_id, pasef_precursors=None):
        if pasef_precursors is None:
            pasef_precursors = []
        self.source = source
        self.id = id
        self.accumulation_time = accumulation_time
        self.max_intensity = max_intensity
        self.msms_type = msms_type
        self.mz_calibration = mz_calibration
        self.num_peaks = num_peaks
        self.num_scans = num_scans
        self.polarity = polarity
        self.property_group = property_group
        self.ramp_time = ramp_time
        self.scan_mode = scan_mode
        self.summed_intensities = summed_intensities
        self.t1 = t1
        self.t2 = t2
        self.time = time
        self.tims_calibration = tims_calibration
        self.tims_id = tims_id
        self.pasef_precursors = pasef_precursors

    def get_scan(self, start, end=None):
        if end is None:
            return self.source.get_scan_by_id("frame=%d scan=%d" % (self.id, start))
        return self.source.get_scan_by_id("frame=%d startScan=%d endScan=%d" % (self.id, start, end))

    @classmethod
    def from_query(cls, source, rowdict):
        frame = cls(
            source,
            rowdict['Id'],
            rowdict['AccumulationTime'],
            rowdict['MaxIntensity'],
            rowdict['MsMsType'],
            rowdict['MzCalibration'],
            rowdict['NumPeaks'],
            rowdict['NumScans'],
            rowdict['Polarity'],
            rowdict['PropertyGroup'],
            rowdict['RampTime'],
            rowdict['ScanMode'],
            rowdict['SummedIntensities'],
            rowdict['T1'],
            rowdict['T2'],
            rowdict['Time'],
            rowdict['TimsCalibration'],
            rowdict['TimsId'],
        )
        return frame


class PASEFPrecursorInformation(Base):
    def __init__(self, frame_id, start_scan, end_scan, isolation_mz, isolation_width, collision_energy, monoisotopic_mz,
                 charge, average_scan_number, intensity, parent):
        self.frame_id = frame_id
        self.start_scan = start_scan
        self.end_scan = end_scan
        self.isolation_mz = isolation_mz
        self.isolation_width = isolation_width
        self.collision_energy = collision_energy
        self.monoisotopic_mz = monoisotopic_mz
        self.charge = charge
        self.average_scan_number = average_scan_number
        self.intensity = intensity
        self.parent = parent

    @property
    def neutral_mass(self):
        return calculate_neutral_mass(self.monoisotopic_mz, self.charge)

    @classmethod
    def from_query(cls, rowdict):
        pasef_pinfo = cls(
            rowdict['Frame'],
            rowdict['ScanNumBegin'],
            rowdict['ScanNumEnd'],
            rowdict['IsolationMz'],
            rowdict['IsolationWidth'],
            rowdict['CollisionEnergy'],
            rowdict['MonoisotopicMz'],
            rowdict['Charge'],
            rowdict['ScanNumber'],
            rowdict['Intensity'],
            rowdict['Parent']
        )
        return pasef_pinfo


class TIMSSpectrumDataBase(Base):
    def __init__(self, frame, start_scan, end_scan=None):
        self.frame = frame
        self.start_scan = start_scan
        self._end_scan = end_scan

    @property
    def end_scan(self):
        if self._end_scan is None:
            return self.start_scan + 1
        else:
            return self._end_scan

    def is_combined(self):
        if self._end_scan is not None and self._end_scan - self.start_scan > 1:
            return True
        return False

    def make_id_string(self):
        if self._end_scan is None:
            return "frame=%d scan=%d" % (self.frame.id, self.start_scan + 1)
        else:
            return "frame=%d startScan=%d endScan=%d" % (self.frame.id, self.start_scan + 1, self.end_scan + 1)


class TIMSPASEFSpectrumData(TIMSSpectrumDataBase):
    def __init__(self, frame, start_scan, pasef_precursor, end_scan=None):
        super(TIMSPASEFSpectrumData, self).__init__(
            frame, start_scan, end_scan)
        self.pasef_precursor = pasef_precursor


default_scan_merging_parameters = {
    "fwhm": 0.04,
    "dx": 0.001
}


class TIMSAPI(object):
    initial_frame_buffer_size = 128

    def _convert_callback(self, frame_id, input_data, func):
        if isinstance(input_data, np.ndarray) and input_data.dtype == np.float64:
            # already supports buffer protocol, no extra copy
            in_array = input_data
        else:
            # convert data to appropriate float data buffer
            in_array = np.array(input_data, dtype=np.float64)
        cnt = len(in_array)
        out = np.empty(shape=cnt, dtype=np.float64)
        success = func(self.handle, frame_id,
                       in_array.ctypes.data_as(POINTER(c_double)),
                       out.ctypes.data_as(POINTER(c_double)),
                       cnt)
        if success == 0:
            throw_tims_error(self.dll)
        return out

    def index_to_mz(self, frame_id, mzs):
        return self._convert_callback(frame_id, mzs, self.dll.tims_index_to_mz)

    def mz_to_index(self, frame_id, mzs):
        return self._convert_callback(frame_id, mzs, self.dll.tims_mz_to_index)

    def scan_number_to_one_over_K0(self, frame_id, mzs):
        return self._convert_callback(frame_id, mzs, self.dll.tims_scannum_to_oneoverk0)

    def one_over_K0_to_scan_number(self, frame_id, mzs):
        return self._convert_callback(frame_id, mzs, self.dll.tims_oneoverk0_to_scannum)

    def scan_number_to_voltage(self, frame_id, mzs):
        return self._convert_callback(frame_id, mzs, self.dll.tims_scannum_to_voltage)

    def voltage_to_scan_number(self, frame_id, mzs):
        return self._convert_callback(frame_id, mzs, self.dll.tims_voltage_to_scannum)
    # Output: list of tuples (indices, intensities)

    def read_scans(self, frame_id, scan_begin, scan_end):
        # buffer-growing loop
        while True:
            cnt = int(self.initial_frame_buffer_size)
            buf = np.empty(shape=cnt, dtype=np.uint32)
            buffer_size_in_bytes = 4 * cnt

            required_len = self.dll.tims_read_scans_v2(self.handle, frame_id, scan_begin, scan_end,
                                                       buf.ctypes.data_as(
                                                           POINTER(c_uint32)),
                                                       buffer_size_in_bytes)
            if required_len == 0:
                throw_tims_error(self.dll)

            if required_len > buffer_size_in_bytes:
                if required_len > 16777216:
                    # arbitrary limit for now...
                    raise RuntimeError("Maximum expected frame size exceeded.")
                self.initial_frame_buffer_size = required_len / 4 + 1  # grow buffer
            else:
                break

        result = []
        d = scan_end - scan_begin
        for i in range(scan_begin, scan_end):
            npeaks = buf[i-scan_begin]
            indices = buf[d:d + npeaks]
            d += npeaks
            intensities = buf[d:d + npeaks]
            d += npeaks
            result.append((indices, intensities))

        return result

    def read_spectrum(self, frame_id, scan_begin, scan_end):
        scans = self.read_scans(frame_id, scan_begin, scan_end)
        # Summarize on a grid
        allind = []
        allint = np.array([], dtype=float)
        for scan in scans:
            indices = np.array(scan[0])
            if len(indices) > 0:
                intens = scan[1]
                allind = np.concatenate((allind, indices))
                allint = np.concatenate((allint, intens))
        allmz = self.index_to_mz(frame_id, allind)
        return allmz, allint

    def read_spectrum_grid(self, frame_id, scan_begin, scan_end):
        scans = self.read_scans(frame_id, scan_begin, scan_end)
        # Summarize on a grid
        scan_numbers = []
        scan_mzs = []
        scan_intensities = []
        for i, scan in enumerate(scans, scan_begin):
            indices = np.array(scan[0])
            if len(indices) > 0:
                intens = scan[1]
                mzs = self.index_to_mz(frame_id, indices)
                scan_numbers.append(i)
                scan_mzs.append(mzs)
                scan_intensities.append(intens)
            else:
                scan_numbers.append(i)
                scan_mzs.append(np.array([], float))
                scan_intensities.append(np.array([], float))
        return scan_numbers, scan_mzs, scan_intensities

    def _decode_activation_method(self, frame):
        mode = frame.scan_mode
        if mode in (2, 8, 9):
            method = dissociation_methods["collision-induced dissociation"]
        elif mode in (3, 4, 5):
            method = dissociation_methods['in-source collision-induced dissociation']
        else:
            print("Unknown Scan Mode %d, Unknown Dissociation. Returning CID" % (mode, ))
            method = dissociation_methods["collision-induced dissociation"]
        return method

    def _query_frame_by_id_number(self, frame_id):
        cursor = self.connection.execute(
            "SELECT * FROM Frames WHERE Id={0};".format(frame_id))
        frame = TIMSFrameInformation.from_query(self, dict(cursor.fetchone()))
        # MS1
        if frame.msms_type == 0:
            pass
        # PASEF MS2
        elif frame.msms_type == 8:
            pasef_cursor = self.connection.execute(
                """SELECT Frame, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy, MonoisotopicMz,
                          Charge, ScanNumber, Intensity, Parent FROM PasefFrameMsMsInfo f JOIN Precursors p on p.id=f.precursor
                          WHERE Frame = ?
                          ORDER BY ScanNumBegin
                          """, (frame_id, ))
            frame.pasef_precursors.extend(
                map(PASEFPrecursorInformation.from_query, pasef_cursor))
        else:
            warnings.warn("No support for MSMSType %r yet" %
                          (frame.msms_type, ))
        return frame


class BrukerTIMSScanDataSource(RandomAccessScanSource):
    _scan_merging_parameters = default_scan_merging_parameters.copy()

    def _is_profile(self, scan):
        if scan.is_combined():
            return True
        return False

    def _scan_time(self, scan):
        return scan.frame.time

    def _scan_id(self, scan):
        return scan.make_id_string()

    def _scan_title(self, scan):
        return self._scan_id(scan)

    def _ms_level(self, scan):
        if scan.frame.msms_type == 0:
            return 1
        return 2

    def _polarity(self, scan):
        if scan.frame.polarity == "+":
            return 1
        elif scan.frame.polarity == '-':
            return -1
        else:
            return None

    def _activation(self, scan):
        method = self._decode_activation_method(scan.frame)
        precursor = self._locate_pasef_precursor_for(scan)
        if precursor is not None:
            collision_energy = precursor.collision_energy
        return ActivationInformation(method, collision_energy)

    def _locate_pasef_precursor_for(self, scan):
        if scan.is_combined():
            # raise ValueError("Cannot determine precursor for combined spectra yet")
            query_interval = Interval(scan.start_scan, scan.end_scan)
            matches = []
            for precursor in scan.frame.pasef_precursors:
                if query_interval.overlaps(Interval(precursor.start_scan, precursor.end_scan)):
                    matches.append(precursor)
            n_matches = len(matches)
            if n_matches == 0:
                return None
            elif n_matches > 1:
                raise ValueError("Multiple precursors found for scan interval!")
            else:
                return matches[0]

        else:
            scan_number = scan.start_scan
            for precursor in scan.frame.pasef_precursors:
                if precursor.start_scan <= scan_number < precursor.end_scan:
                    return precursor
            return None

    def _precursor_information(self, scan):
        if scan.frame.msms_type == 8:
            precursor = self._locate_pasef_precursor_for(scan)
            if precursor is not None:
                mz = precursor.monoisotopic_mz
                intensity = precursor.intensity
                charge = precursor.charge
                parent_frame_id = precursor.parent

                if not scan.is_combined():
                    current_drift_time = self.scan_number_to_one_over_K0(scan.frame.id, [scan.start_scan])
                else:
                    current_drift_time = self.scan_number_to_one_over_K0(
                        scan.frame.id, [np.ceil(precursor.average_scan_number)])

                parent_frame = self.get_frame_by_id(parent_frame_id)._data
                precursor_scan_id = TIMSPASEFSpectrumData(parent_frame, precursor.start_scan, None, precursor.end_scan).make_id_string()
                product_scan_id = scan.make_id_string()
                pinfo = PrecursorInformation(
                    mz, intensity, charge, precursor_scan_id,
                    product_scan_id=product_scan_id, source=self, annotations={
                        inverse_reduced_ion_mobility: current_drift_time,
                    })
                return pinfo
        return None

    def _isolation_window(self, scan):
        if scan.frame.msms_type == 8:
            precursor = self._locate_pasef_precursor_for(scan)
            if precursor is not None:
                width = precursor.isolation_width / 2
                window = IsolationWindow(width, precursor.isolation_mz, width)
                return window
        return None

    def _scan_index(self, scan):
        cursor = self.connection.execute("SELECT sum(NumScans) FROM Frames WHERE Id < ?", (scan.frame.id, ))
        result = cursor.fetchone()[0]
        if result is None:
            result = 0
        return result + scan.start_scan

    def _acquisition_information(self, scan):
        idx = self._scan_index(scan) - \
            self._scan_count_to_frame_id[scan.frame.id - 1, 0]
        drift_time = self.scan_number_to_one_over_K0(scan.frame.id, [idx])[0]
        scan_time = self._scan_time(scan)
        evt = ScanEventInformation(scan_time, [
            ScanWindow(self._acquisition_parameters['scan_window_lower'],
                       self._acquisition_parameters['scan_window_upper'])
                       ])
        evt._ion_mobility.add_ion_mobility(inverse_reduced_ion_mobility, drift_time)
        return ScanAcquisitionInformation('none', [evt])

    def _get_centroids(self, scan):
        mzs, intensities = self.read_spectrum(
            scan.frame.id, scan.start_scan, scan.end_scan)
        sort_mask = np.argsort(mzs)
        mzs = mzs[sort_mask]
        intensities = intensities[sort_mask]
        centroids = pick_peaks(mzs, intensities, peak_mode="centroid")
        if centroids is None:
            centroids = PeakIndex(
                np.array([], float), np.array([], float), [])
        return centroids

    def _scan_arrays(self, scan):
        if scan.is_combined():
            mzs, intensities = self.read_spectrum(
                scan.frame.id, scan.start_scan, scan.end_scan)
            if len(mzs) == 0:
                return np.array([], dtype=float), np.array([], dtype=float)
            sort_mask = np.argsort(mzs)
            mzs = mzs[sort_mask]
            intensities = intensities[sort_mask]
            centroids = pick_peaks(mzs, intensities, peak_mode="centroid")
            if centroids is None:
                return np.array([], dtype=float), np.array([], dtype=float)
            mzs, intensities = reprofile(
                centroids, dx=self._scan_merging_parameters['dx'],
                override_fwhm=self._scan_merging_parameters['fwhm'])
            return mzs, intensities
        else:
            mzs, intensities = self.read_spectrum(scan.frame.id, scan.start_scan, scan.end_scan)
            return mzs, intensities


class BrukerTIMSFrameSource(IonMobilitySourceRandomAccessFrameSource):
    def _frame_id(self, data):
        return "frame=%d startScan=%d endScan=%d" % (data.id, 1, data.num_scans)

    def _frame_index(self, data):
        return data.id - 1

    def _frame_time(self, data):
        return data.time

    def _frame_ms_level(self, data):
        if data.msms_type == 0:
            return 1
        return 2

    def _frame_start_scan_index(self, data):
        start_i, frame_id = self._scan_count_to_frame_id[data.id - 1, :]
        assert frame_id == data.id
        return start_i

    def _frame_end_scan_index(self, data):
        start_i, frame_id = self._scan_count_to_frame_id[data.id - 1, :]
        assert frame_id == data.id
        return start_i + data.num_scans

    def _frame_precursor_information(self, data):
        pinfos = []
        for pasef in data.pasef_precursors:
            mz = pasef.monoisotopic_mz
            inten = pasef.intensity
            charge = pasef.charge
            parent_merged_id = 'frame=%d startScan=%d endScan=%d' % (pasef.parent, pasef.start_scan, pasef.end_scan)
            product_merged_id = 'frame=%d startScan=%d endScan=%d' % (pasef.frame_id, pasef.start_scan, pasef.end_scan)
            pinfos.append(PrecursorInformation(
                mz, inten, charge, parent_merged_id, self,
                product_scan_id=product_merged_id, annotations={
                    "parent_frame_index": pasef.parent - 1,
                    "start_scan": pasef.start_scan,
                    "end_scan": pasef.end_scan,
                }))
        if not pinfos:
            return None
        return pinfos

    def _frame_activation(self, data):
        activations = []
        method = self._decode_activation_method(data)
        for pasef in data.pasef_precursors:
            collision_energy = pasef.collision_energy
            activations.append(ActivationInformation(method, collision_energy))
        if not activations:
            return None
        return activations

    def _frame_isolation_window(self, data):
        isolations = []
        for pasef in data.pasef_precursors:
            width = pasef.isolation_width / 2.
            isolations.append(IsolationWindow(width, pasef.isolation_mz, width))
        if not isolations:
            return None
        return isolations

    def _frame_polarity(self, data):
        if data.polarity == "+":
            return 1
        elif data.polarity == '-':
            return -1
        else:
            warnings.warn("Unknown polarity %r" % (data.polarity, ))
            return 1

    def get_frame_by_time(self, frame_time):
        i = np.searchsorted(self._scan_time_to_frame_id[:, 0], frame_time) - 1
        _time, fi = self._scan_time_to_frame_id[i, :]
        return self.get_frame_by_id(fi)

    def get_frame_by_index(self, index):
        return self.get_frame_by_id(index + 1)

    def get_frame_by_id(self, frame_id):
        if isinstance(frame_id, str):
            match = multi_scan_id_parser.search(frame_id)
            if match is None:
                raise KeyError(frame_id)
            frame_id = int(match.group(1))
        if frame_id in self._frame_cache:
            return self._frame_cache[frame_id]
        frame = self._query_frame_by_id_number(frame_id)
        frame = self._make_frame(frame)
        self._frame_cache[frame_id] = frame
        return frame

    def _ms2_frames_for_parent_id(self, parent_frame_id):
        product_frame_ids = self.connection.execute("""
            SELECT Frame
            FROM PasefFrameMsMsInfo JOIN Precursors ON Precursors.id = PasefFrameMsMsInfo.precursor
            WHERE Parent = ?;""", (parent_frame_id, )).fetchall()
        product_frames = [self.get_frame_by_id(
            i) for i, in product_frame_ids]
        return product_frames

    def _validate_frame(self, data):
        return True

    def _make_frame(self, data):
        return IonMobilityFrame(data, self)

    def _cache_frame(self, frame):
        pass

    def _default_frame_iterator(self, start_index=None):
        if start_index is None:
            start_index = 0
        for i in range(start_index, self.frame_count()):
            yield self._query_frame_by_id_number(i + 1)

    def make_frame_iterator(self, iterator=None, grouped=False):
        from ms_deisotope.data_source.scan.scan_iterator import (
            _SingleScanIteratorImpl, _InterleavedGroupedScanIteratorImpl)
        if iterator is None:
            iterator = self._default_frame_iterator()
        if grouped:
            strategy = _InterleavedGroupedScanIteratorImpl(
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
        if require_ms1:
            base_frame = self.get_frame_by_index(start_index)
            if base_frame.ms_level > 1:
                precursors = base_frame.precursor_information

            while start_index != 0:
                frame = self.get_frame_by_index(start_index)
                if frame.ms_level > 1:
                    start_index -= 1
                else:
                    break
        iterator = self.make_frame_iterator(
            self._default_frame_iterator(start_index), grouped=grouped)
        return iterator


single_scan_id_parser = re.compile(r"frame=(\d+) scan=(\d+)")
multi_scan_id_parser = re.compile(r"frame=(\d+) startScan=(\d+) endScan=(\d+)")


class BrukerTIMSLoader(BrukerTIMSFrameSource, BrukerTIMSScanDataSource, TIMSMetadata, TIMSAPI):

    def __init__(self, analysis_directory, use_recalibrated_state=False, scan_merging_parameters=None):
        if sys.version_info.major == 2:
            if isinstance(analysis_directory, str):
                analysis_directory = analysis_directory.decode('utf8')
            if not isinstance(analysis_directory, unicode):
                raise ValueError("analysis_directory must be a Unicode string.")
        if sys.version_info.major == 3:
            if isinstance(analysis_directory, bytes):
                analysis_directory = analysis_directory.decode('utf8')
            if not isinstance(analysis_directory, str):
                raise ValueError("analysis_directory must be a string.")
        if scan_merging_parameters is None:
            scan_merging_parameters = default_scan_merging_parameters.copy()
        else:
            for key, value in default_scan_merging_parameters.items():
                scan_merging_parameters.setdefault(key, value)

        self.dll = load_library()

        self.handle = self.dll.tims_open(
            analysis_directory.encode('utf-8'), 1 if use_recalibrated_state else 0)
        if self.handle == 0:
            throw_tims_error(self.dll)

        self._source_file_name = analysis_directory
        self.connection = sqlite3.connect(os.path.join(analysis_directory, "analysis.tdf"))
        self.connection.row_factory = sqlite3.Row

        self.initialize_frame_cache()
        self.initialize_scan_cache()

        self.initial_frame_buffer_size = 128 # may grow in readScans()
        self._read_metadata()
        self._scan_merging_parameters = scan_merging_parameters
        self._producer = self.make_frame_iterator(grouped=True)

    @property
    def source_file_name(self):
        return self._source_file_name

    def __del__(self):
        if hasattr(self, 'handle'):
            self.dll.tims_close(self.handle)

    def _describe_frame(self, frame_id):
        cursor = self.connection.execute("SELECT * FROM Frames WHERE Id={0};".format(frame_id))
        return dict(cursor.fetchone())

    def __len__(self):
        return self.frame_count()

    def next(self):
        return next(self._producer)

    def get_scan_by_index(self, scan_index):
        frame_index = np.searchsorted(
            self._scan_count_to_frame_id[:, 0], scan_index + 1)
        if frame_index > 0:
            frame_index -= 1
        c, i = self._scan_count_to_frame_id[frame_index, :]
        frame_id = i
        remainder = scan_index - c + 1
        scan_number = remainder
        scan_id = "frame=%d scan=%d" % (frame_id, scan_number)
        return self.get_scan_by_id(scan_id)

    def get_scan_by_time(self, scan_time):
        i = np.searchsorted(self._scan_time_to_frame_id[:, 0], scan_time) - 1
        time, fi = self._scan_time_to_frame_id[i, :]
        remainder = scan_time - time
        try:
            next_time = self._scan_time_to_frame_id[i + 1, 0]
        except IndexError:
            # An instance with exactly one frame?
            if i == 0:
                next_time = time * 2.0
            else:
                raise
        duration = next_time - time
        n_scans = self.get_frame_by_id(fi)._data.num_scans
        grid_space = np.linspace(0, duration, n_scans)
        i = np.searchsorted(grid_space, remainder) - 1
        scan_id = 'frame=%d scan=%d' % (fi, i)
        return self.get_scan_by_id(scan_id)

    def get_scan_by_id(self, scan_id):
        match = single_scan_id_parser.match(scan_id)
        if match is None:
            match = multi_scan_id_parser.match(scan_id)
            if match is None:
                raise ValueError("%r does not look like a TIMS nativeID" % (scan_id, ))
            else:
                frame_id, start_scan, end_scan = map(int, match.groups())
                frame = self.get_frame_by_index(frame_id - 1)._data
                if frame.msms_type == 8:
                    pasef_scan = TIMSPASEFSpectrumData(
                        frame, start_scan - 1, None, end_scan)
                    pasef_scan.pasef_precursor = self._locate_pasef_precursor_for(pasef_scan)
                    scan_obj = self._make_scan(pasef_scan)
                else:
                    scan_obj = self._make_scan(TIMSSpectrumDataBase(frame, start_scan - 1, end_scan))
        else:
            frame_id, scan_number = map(int, match.groups())
            frame = self.get_frame_by_index(frame_id - 1)._data
            if frame.msms_type == 8:
                pasef_scan = TIMSPASEFSpectrumData(frame, scan_number - 1, None)
                pasef_scan.pasef_precursor = self._locate_pasef_precursor_for(pasef_scan)
                scan_obj = self._make_scan(pasef_scan)
            else:
                scan_obj = self._make_scan(TIMSSpectrumDataBase(frame, scan_number - 1))
        # Cache scan here
        return scan_obj

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        '''Reconstruct an iterator which will start from the scan matching one of ``scan_id``,
        ``rt``, or ``index``. Only one may be provided.

        After invoking this method, the iterator this object wraps will be changed to begin
        yielding scan bunchs (or single scans if ``grouped`` is ``False``).

        This method will trigger several random-access operations, making it prohibitively
        expensive for normally compressed files.

        Arguments
        ---------
        scan_id: str, optional
            Start from the scan with the specified id.
        rt: float, optional
            Start from the scan nearest to specified time (in minutes) in the run. If no
            exact match is found, the nearest scan time will be found, rounded up.
        index: int, optional
            Start from the scan with the specified index.
        require_ms1: bool, optional
            Whether the iterator must start from an MS1 scan. True by default.
        grouped: bool, optional
            whether the iterator should yield scan bunches or single scans. True by default.
        '''
        raise NotImplementedError()

    def _ms1_frame_iterator(self):
        ms1_frame_ids = self.connection.execute("SELECT Id FROM Frames WHERE MsMsType = 0;").fetchall()
        for frame_id, in ms1_frame_ids:
            yield frame_id

    def _make_default_iterator(self, frame=False):
        if frame:
            return self._ms1_frame_iterator()
        else:
            raise NotImplementedError()
