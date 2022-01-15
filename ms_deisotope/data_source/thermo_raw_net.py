'''Thermo RAW file reading implementation using the pure .NET
RawFileReader library released in 2017.

This module provides :class:`ThermoRawLoader`, a :class:`~.RandomAccessScanSource`
implementation.

Depends upon the ``pythonnet`` project which provides the :mod:`clr`
module, enabling nearly seamless interoperation with the Common Language
Runtime.

The public interface of this module should be identical to
:mod:`ms_deisotope.data_source.thermo_raw`.

.. note::
    This interface was largely based upon the APIs that ProteoWizard used, both
    in order to understand how the Thermo libraries really worked, and to maintain
    parity with it.
'''

import sys
import os

from collections import OrderedDict

import numpy as np

from pyteomics.auxiliary import unitfloat

from six import string_types as basestring


from ms_peak_picker import PeakSet, PeakIndex, simple_peak

from ms_deisotope.data_source.common import (
    PrecursorInformation, ChargeNotProvided, Scan,
    ActivationInformation, MultipleActivationInformation,
    IsolationWindow, ScanDataSource, ScanEventInformation,
    ScanAcquisitionInformation, ScanWindow, RandomAccessScanSource)

from ms_deisotope.data_source._thermo_helper import (
    _InstrumentMethod, ThermoRawScanPtr, FilterString,
    _make_id, _id_template, _RawFileMetadataLoader, analyzer_map)

from ms_deisotope.data_source.metadata.activation import (
    supplemental_term_map, dissociation_methods_map)
from ms_deisotope.data_source.metadata.sample import Sample
from ms_deisotope.data_source.metadata.scan_traits import FAIMS_compensation_voltage

def _try_number(string):
    try:
        x = float(string)
        return x
    except (TypeError, ValueError):
        return string


_DEFAULT_DLL_PATH = os.path.join(
    os.path.dirname(
        os.path.realpath(__file__)),
    "_vendor",
    "ThermoRawFileReader_3_0_41",
    "Libraries")


# late binding imports

Business = None
_RawFileReader = None
clr = None
NullReferenceException = Exception
Marshal = None
IntPtr = None
Int64 = None


def is_thermo_raw_file(path):
    '''Detect whether or not the file referenced by ``path``
    is a Thermo RAW file.

    Parameters
    ----------
    path: :class:`str`
        The path to test

    Returns
    -------
    :class:`bool`:
        Whether or not the file is a Thermo RAW file.
    '''
    if not _test_dll_loaded():
        try:
            register_dll()
        except ImportError:
            return False
    with open(path, 'rb') as fh:
        lead_bytes = fh.read(32)
        decoded = lead_bytes.decode("utf-16")[1:9]
        if decoded == "Finnigan":
            try:
                source = _RawFileReader.RawFileReaderAdapter.FileFactory(path)
                source.SelectInstrument(Business.Device.MS, 1)
                return True
            except NullReferenceException:   # pylint: disable=broad-except
                return False
        return False


def infer_reader(path):
    '''If the file referenced by ``path`` is a Thermo RAW
    file, return the callable (:class:`ThermoRawLoader`) to
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
        If the file is not a Thermo RAW file
    '''
    if is_thermo_raw_file(path):
        return ThermoRawLoader
    raise ValueError("Not Thermo Raw File")


def determine_if_available():
    '''Checks whether or not the .NET-based Thermo
    RAW file reading feature is available.

    Returns
    -------
    :class:`bool`:
        Whether or not the feature is enabled.
    '''
    try:
        return _register_dll([_DEFAULT_DLL_PATH])
    except (OSError, ImportError):
        return False


def _register_dll(search_paths=None):
    '''Start the Common Language Runtime interop service by importing
    the :mod:`clr` module from Pythonnet, and then populate the global
    names referring to .NET entities, and finally attempt to locate the
    ThermoRawFileReader DLLs by searching alogn ``search_paths``.

    Parameters
    ----------
    search_paths: list
        The paths to check along for the ThermoRawFileReader DLL bundle.

    Returns
    -------
    :class:`bool`:
        Whether or not the .NET library successfully loaded
    '''
    from ms_deisotope.config import get_config
    if search_paths is None:
        search_paths = []
    search_paths = list(search_paths)
    search_paths.append(_DEFAULT_DLL_PATH)
    # Take user-specified search paths first.
    search_paths = get_config().get('vendor_readers', {}).get('thermo-net', []) + search_paths
    global _RawFileReader, Business, clr, NullReferenceException   # pylint: disable=global-statement
    global Marshal, IntPtr, Int64   # pylint: disable=global-statement
    if _test_dll_loaded():
        return True
    try:
        import clr  # pylint: disable=redefined-outer-name
        from System import NullReferenceException  # pylint: disable=redefined-outer-name
        clr.AddReference("System.Runtime")
        clr.AddReference("System.Runtime.InteropServices")
        from System import IntPtr, Int64  # pylint: disable=redefined-outer-name
        from System.Runtime.InteropServices import Marshal  # pylint: disable=redefined-outer-name
    except ImportError:
        return False
    for path in search_paths:
        sys.path.append(path)
        try:
            clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
            clr.AddReference('ThermoFisher.CommonCore.Data')
        except OSError:
            continue
        try:
            import ThermoFisher.CommonCore.Data.Business as Business  # pylint: disable=redefined-outer-name
            import ThermoFisher.CommonCore.RawFileReader as _RawFileReader  # pylint: disable=redefined-outer-name
        except ImportError:
            continue
    return _test_dll_loaded()


def register_dll(search_paths=None):
    '''Register the location of the Thermo RawFileReader DLL bundle with
    the Common Language Runtime interop system and load the .NET symbols
    used by this feature.

    Parameters
    ----------
    search_paths: list
        The paths to check along for the ThermoRawFileReader DLL bundle.

    '''
    if search_paths is None:
        search_paths = []
    loaded = _register_dll(search_paths)
    if not loaded:
        msg = '''The ThermoFisher.CommonCore libraries could not be located and loaded.'''
        raise ImportError(msg)


def _test_dll_loaded():
    return _RawFileReader is not None


def _copy_double_array(src):
    '''A quick and dirty implementation of the fourth technique shown in
    https://mail.python.org/pipermail/pythondotnet/2014-May/001525.html for
    copying a .NET Array[Double] to a NumPy ndarray[np.float64] via a raw
    memory copy.

    ``int_ptr_tp`` must be an integer type that can hold a pointer. On Python 2
    this is :class:`long`, and on Python 3 it is :class:`int`.
    '''
    # When the input .NET array pointer is None, return an empty array. On Py2
    # this would happen automatically, but not on Py3, and perhaps not safely on
    # all Py2 because it relies on pythonnet and the .NET runtime properly checking
    # for nulls.
    if src is None:
        return np.array([], dtype=np.float64)
    dest = np.empty(len(src), dtype=np.float64)
    Marshal.Copy(
        src, 0,
        IntPtr.__overloads__[Int64](dest.__array_interface__['data'][0]),
        len(src))
    return dest


class RawReaderInterface(ScanDataSource):
    ''':class:`~.ScanDataSource` implementation for Thermo's RawFileReader API.
    Not intended for direct instantiation.
    '''

    def _scan_arrays(self, scan):
        scan_number = scan.scan_number + 1
        stats = self._source.GetScanStatsForScanNumber(scan_number)
        segscan = self._source.GetSegmentedScanFromScanNumber(scan_number, stats)
        mzs = _copy_double_array(segscan.Positions)
        inten = _copy_double_array(segscan.Intensities)
        return mzs, inten

    def _pick_peaks_vendor(self, scan, *args, **kwargs):
        scan_info = Business.Scan.FromFile(self._source, scan.scan_number + 1)
        if scan_info.HasCentroidStream:
            stream = self._source.GetCentroidStream(scan.scan_number + 1, 0)
            mzs = stream.Masses
            intens = stream.Intensities
            peaks = PeakSet([simple_peak(mzs[i], intens[i], 0.001) for i in range(len(mzs))])
            peaks.reindex()
            arrays = self._scan_arrays(scan)
            return PeakIndex(arrays[0], arrays[1], peaks)
        else:
            raise NotImplementedError()

    def _scan_id(self, scan):
        scan_number = scan.scan_number
        return _make_id(scan_number + 1)

    def _is_profile(self, scan):
        return not self._source.IsCentroidScanFromScanNumber(
            scan.scan_number + 1)

    def _polarity(self, scan):
        filter_string = self._filter_string(scan)
        return filter_string.data['polarity']

    def _scan_title(self, scan):
        return "%s %r" % (self._scan_id(scan), self._filter_string(scan))

    def _filter_string(self, scan):
        if scan.filter_string is None:
            scan_number = scan.scan_number
            scan.filter_string = FilterString(self._source.GetFilterForScanNumber(scan_number + 1).Filter)
        return scan.filter_string

    def _scan_index(self, scan):
        scan_number = scan.scan_number
        return scan_number

    def _scan_time(self, scan):
        scan_number = scan.scan_number
        return self._source.RetentionTimeFromScanNumber(scan_number + 1)

    def _ms_level(self, scan):
        scan_number = scan.scan_number
        f = self._source.GetFilterForScanNumber(scan_number + 1)
        return f.MSOrder

    def _isolation_window(self, scan):
        scan_number = scan.scan_number
        ms_level = self._ms_level(scan)
        width = 0
        trailer = self._trailer_values(scan)
        filt = self._source.GetFilterForScanNumber(scan_number + 1)
        seq_index = filt.MSOrder - 2
        try:
            # Fetch the isolation window width from the old location first, which
            # will be correct on old files, where the new API won't be right.
            width = trailer['MS%d Isolation Width' % ms_level]
        except KeyError:
            # Fall back to the new API, which is akin to our only hope here?
            width = filt.GetIsolationWidth(seq_index)
        width /= 2.0
        offset = filt.GetIsolationWidthOffset(seq_index)
        precursor_mz = filt.GetMass(seq_index)
        return IsolationWindow(width, precursor_mz + offset, width)

    def _trailer_values(self, scan):
        if scan.trailer_values is not None:
            return scan.trailer_values
        scan_number = scan.scan_number
        trailers = self._source.GetTrailerExtraInformation(scan_number + 1)
        scan.trailer_values = OrderedDict(
            zip([label.strip(":") for label in trailers.Labels], map(_try_number, trailers.Values)))
        return scan.trailer_values

    def _precursor_information(self, scan):
        scan_number = scan.scan_number
        filt = self._source.GetFilterForScanNumber(scan_number + 1)
        precursor_mz = filt.GetMass(filt.MSOrder - 2)
        trailers = self._trailer_values(scan)
        _precursor_mz = float(trailers.get("Monoisotopic M/Z", 0))
        if _precursor_mz > 0:
            precursor_mz = _precursor_mz

        # imitate proteowizard's firmware bug correction
        isolation_window = self._isolation_window(scan)
        if (isolation_window.upper + isolation_window.lower) / 2 <= 2.0:
            if (isolation_window.target - 3.0 > precursor_mz) or (isolation_window.target + 2.5 < precursor_mz):
                precursor_mz = isolation_window.target
        elif precursor_mz not in isolation_window:
            precursor_mz = isolation_window.target
        charge = int(trailers.get("Charge State", 0))
        if charge == 0:
            charge = ChargeNotProvided
        inten = 0
        precursor_scan_number = None
        precursor_scan_number = trailers.get('Master Scan Number')
        if precursor_scan_number == 0:
            precursor_scan_number = None
        if precursor_scan_number is not None:
            precursor_scan_number = int(precursor_scan_number) - 1
            if self.get_scan_by_index(precursor_scan_number).ms_level >= self._ms_level(scan):
                precursor_scan_number = None
        if precursor_scan_number is None:
            current_level = self._ms_level(scan)
            # We want to start looking for the most recent spectrum with the next lowest MS
            # level. The expecteation is that we couldn't determine the precursor scan accurately,
            # so we'll probe backwards until the next best candidate.
            current_index = self._scan_index(scan)
            # Use the index of the current scan to check the index of the previous scans at a lower
            # MS level. This may be None, in which case there is no earlier scan recorded.
            lookup = self._get_previous_scan_index_for_ms_level(current_index, current_level - 1)
            if lookup is not None:
                last_index = lookup
            else:
                last_index = current_index - 1
            i = 0
            while last_index >= 0 and i < 100:
                prev_scan = self.get_scan_by_index(last_index)
                if prev_scan.ms_level >= current_level:
                    last_index -= 1
                else:
                    precursor_scan_number = prev_scan._data.scan_number
                    break
                i += 1
        if precursor_scan_number is not None:
            precursor_scan_id = self.get_scan_by_index(precursor_scan_number).id
        else:
            import warnings
            warnings.warn("Could not resolve precursor scan for %s" % (self._scan_id(scan), ))
            precursor_scan_id = None
        return PrecursorInformation(
            precursor_mz, inten, charge, precursor_scan_id,
            source=self, product_scan_id=self._scan_id(scan))

    def _get_scan_segment(self, scan):
        trailer = self._trailer_values(scan)
        try:
            return int(trailer['Scan Segment'])
        except KeyError:
            return 1

    def _get_scan_event(self, scan):
        trailer = self._trailer_values(scan)
        try:
            return int(trailer['Scan Event'])
        except KeyError:
            return 1

    def _activation(self, scan):
        filter_string = self._filter_string(scan)
        tandem_sequence = filter_string.get("tandem_sequence")
        # If the tandem sequence exists, the last entry is the most recent tandem acquisition.
        # It will list contain one or more activation types. Alternatively, multiple activations
        # of the same precursor may exist in the list as separate events in the tandem sequence.
        if tandem_sequence is not None:
            activation_event = tandem_sequence[-1]
            activation_type = list(activation_event.get("activation_type"))
            has_supplemental_activation = filter_string.get("supplemental_activation")

            if activation_type is not None:
                energy = list(activation_event.get("activation_energy"))
                if len(tandem_sequence) > 1:
                    prev_event = tandem_sequence[-2]
                    # Merge previous tandem sequences of the same precursor
                    if abs(prev_event['isolation_mz'] - activation_event['isolation_mz']) < 1e-3:
                        activation_type = list(prev_event.get("activation_type")) + activation_type
                        energy = list(prev_event.get("activation_energy")) + energy
                        has_supplemental_activation = True

                if has_supplemental_activation and len(activation_type) > 1:
                    activation_type.append(supplemental_term_map[
                        dissociation_methods_map[activation_type[-1]]])
                if len(activation_type) == 1:
                    return ActivationInformation(activation_type[0], energy[0])
                else:
                    return MultipleActivationInformation(activation_type, energy)
        return None

    def _acquisition_information(self, scan):
        fline = self._filter_string(scan)
        event = self._get_scan_event(scan)
        trailer_extras = self._trailer_values(scan)
        traits = {
            'preset scan configuration': event,
            'filter string': fline,
        }
        cv = fline.get("compensation_voltage")
        if cv is not None:
            traits[FAIMS_compensation_voltage] = cv
        event = ScanEventInformation(
            self._scan_time(scan),
            injection_time=unitfloat(trailer_extras.get('Ion Injection Time (ms)', 0.0), 'millisecond'),
            window_list=[ScanWindow(
                fline.get("scan_window")[0], fline.get("scan_window")[1])], traits=traits)
        return ScanAcquisitionInformation("no combination", [event])

    def _instrument_configuration(self, scan):
        fline = self._filter_string(scan)
        try:
            confid = self._analyzer_to_configuration_index[analyzer_map[fline.data.get("analyzer")]]
            return self._instrument_config[confid]
        except KeyError:
            return None

    def _annotations(self, scan):
        fline = self._filter_string(scan)
        trailer_extras = self._trailer_values(scan)
        annots = {
            "filter string": fline,
        }
        microscans = trailer_extras.get("Micro Scan Count")
        if microscans is not None:
            annots['[Thermo Trailer Extra]Micro Scan Count'] = float(microscans)
        scan_segment = trailer_extras.get("Scan Segment")
        if scan_segment is not None:
            annots['[Thermo Trailer Extra]Scan Segment'] = int(scan_segment)
        scan_event = trailer_extras.get("Scan Event")
        if scan_event is not None:
            annots['[Thermo Trailer Extra]Scan Event'] = int(scan_event)
        mono_mz = float(trailer_extras.get("Monoisotopic M/Z", 0))
        if mono_mz is not None and mono_mz > 0:
            annots['[Thermo Trailer Extra]Monoisotopic M/Z'] = mono_mz
        hcd_ev = _trailer_float(trailer_extras.get('HCD Energy eV'))
        if hcd_ev is not None and hcd_ev > 0:
            annots['[Thermo Trailer Extra]HCD Energy eV'] = float(hcd_ev)
        hcd_energies = trailer_extras.get('HCD Energy')
        if hcd_energies is not None:
            if isinstance(hcd_energies, basestring) and not hcd_energies.strip():
                pass
            else:
                annots['[Thermo Trailer Extra]HCD Energy'] = hcd_energies
        return annots


def _trailer_float(value):
    if value is None:
        return None
    try:
        value = value.strip()
    except AttributeError:
        return value
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


class ThermoRawLoader(RawReaderInterface, RandomAccessScanSource, _RawFileMetadataLoader):
    '''Reads scans from Thermo Fisher RAW files directly. Provides both iterative and
    random access.
    '''

    def __init__(self, source_file, _load_metadata=True, **kwargs):
        if not _test_dll_loaded():
            register_dll()
        self._source = _RawFileReader.RawFileReaderAdapter.FileFactory(source_file)
        self._source.SelectInstrument(Business.Device.MS, 1)
        self.source_file = source_file
        self._producer = None
        self._scan_type_index = dict()

        self._method = None
        self._analyzer_to_configuration_index = {}
        self._instrument_config = {}

        self.make_iterator()
        self.initialize_scan_cache()
        self._first_scan_time = self.get_scan_by_index(0).scan_time
        self._last_scan_time = self.get_scan_by_index(
            self._source.RunHeaderEx.LastSpectrum - 1).scan_time
        self._index = self._pack_index()

        if _load_metadata:
            self._method = self._parse_method()
            self._build_scan_type_index()
            self._get_instrument_info()

    def _has_ms1_scans(self):
        if self._scan_type_index:
            return 1 in self._scan_type_index
        else:
            # metadata has not been loaded so best to assume there is
            return True

    def _has_msn_scans(self):
        if self._scan_type_index:
            return max(self._scan_type_index) > 1
        else:
            # metadata has not been loaded so best to assume there is
            return True

    def has_msn_scans(self):
        return self._has_msn_scans()

    def has_ms1_scans(self):
        return self._has_ms1_scans()

    def __reduce__(self):
        return self.__class__, (self.source_file, False), self.__getstate__()

    def __getstate__(self):
        state = {
            "method": self._method,
            "scan_type_index": self._scan_type_index,
            "analyzer_to_configuration_index": self._analyzer_to_configuration_index,
            "instrument_config": self._instrument_config,
            "previous_ms_levels": self._previous_ms_levels,
        }
        return state

    def __setstate__(self, state):
        self._method = state['method']
        self._scan_type_index = state['scan_type_index']
        self._analyzer_to_configuration_index = state['analyzer_to_configuration_index']
        self._instrument_config = state['instrument_config']
        self._previous_ms_levels = state['previous_ms_levels']

    @property
    def index(self):
        '''Accesses the scan index

        Returns
        -------
        :class:`collections.OrderedDict`
            Maps scan ID to index
        '''
        return self._index

    def __len__(self):
        return len(self.index)

    def __repr__(self):
        return "ThermoRawLoader(%r)" % (self.source_file)

    def close(self):
        '''Close the underlying file reader.
        '''
        if self._source is not None:
            self._source.Close()
            self._source = None
        self._dispose()

    def __del__(self):
        self.close()

    def reset(self):
        self.make_iterator(None)
        self.initialize_scan_cache()

    def _pack_index(self):
        index = OrderedDict()
        for sn in range(self._source.RunHeaderEx.FirstSpectrum,
                        self._source.RunHeaderEx.LastSpectrum + 1):
            index[_make_id(sn)] = sn
        return index

    def _get_instrument_model_name(self):
        return self._source.GetInstrumentData().Model

    def _get_instrument_serial_number(self):
        return self._source.GetInstrumentData().SerialNumber

    def _parse_method(self):
        try:
            method_count = self._source.InstrumentMethodsCount
            if method_count == 0:
                return _InstrumentMethod('')
            # the data acquisition method should be the last method
            return _InstrumentMethod(self._source.GetInstrumentMethod(method_count - 1))
        except NullReferenceException:   # pylint: disable=broad-except
            return _InstrumentMethod('')

    def _scan_time_to_scan_number(self, scan_time):
        scan_number = self._source.ScanNumberFromRetentionTime(scan_time) - 1
        return scan_number

    def samples(self):
        """Describe the sample(s) used to generate the mass spectrometry
        data contained in this file.

        Returns
        -------
        :class:`list` of :class:`~.Sample`
        """
        result = []
        si = self._source.SampleInformation
        sample = Sample(si.SampleId or 'sample_1')
        sample.name = si.SampleName or si.SampleId
        if si.SampleVolume:
            sample.parameters['sample volume'] = si.SampleVolume
        if si.SampleWeight:
            sample.parameters['sample mass'] = si.SampleWeight
        if si.Vial:
            sample.parameters['sample vial'] = si.Vial
        if si.Barcode:
            sample.parameters['sample barcode'] = si.Barcode
        result.append(sample)
        return result

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
        scan_number = int(index)
        try:
            return self._scan_cache[scan_number]
        except KeyError:
            package = ThermoRawScanPtr(scan_number)
            if not package.validate(self):
                raise IndexError(index)
            scan = Scan(package, self)
            self._scan_cache[scan_number] = scan
            return scan

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
        scan_number = int(str(scan_id).replace(_id_template, '')) - 1
        try:
            return self._scan_cache[scan_number]
        except KeyError:
            package = ThermoRawScanPtr(scan_number)
            if not package.validate(self):
                raise KeyError(str(scan_id))
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
        scan_number = self._scan_time_to_scan_number(time)
        try:
            return self._scan_cache[scan_number]
        except KeyError:
            package = ThermoRawScanPtr(scan_number)
            scan = Scan(package, self)
            self._scan_cache[scan_number] = scan
            return scan

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        '''Reconstruct an iterator which will start from the scan matching one of ``scan_id``,
        ``rt``, or ``index``. Only one may be provided.

        After invoking this method, the iterator this object wraps will be changed to begin
        yielding scan bunchs (or single scans if ``grouped`` is ``False``).

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
        if scan_id is not None:
            scan_number = int(str(scan_id).replace(_id_template, '')) - 1
        elif index is not None:
            scan_number = int(index)
        elif rt is not None:
            scan_number = self._scan_time_to_scan_number(rt)
        if require_ms1:
            start_index = scan_number
            while start_index != 0:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            scan_number = start_index
        iterator = self._make_pointer_iterator(start_index=scan_number)
        if grouped:
            self._producer = self._scan_group_iterator(iterator)
        else:
            self._producer = self._single_scan_iterator(iterator)
        return self

    def _make_scan_index_producer(self, start_index=None, start_time=None):
        if start_index is not None:
            return range(start_index, self._source.RunHeaderEx.LastSpectrum)
        elif start_time is not None:
            start_index = self._scan_time_to_scan_number(start_time)
            while start_index != 0:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            return range(start_index, self._source.RunHeaderEx.LastSpectrum)
        else:
            return range(0, self._source.RunHeaderEx.LastSpectrum)

    def _make_pointer_iterator(self, start_index=None, start_time=None):
        iterator = self._make_scan_index_producer(start_index, start_time)
        for i in iterator:
            yield ThermoRawScanPtr(i)

    def _make_default_iterator(self):
        return self._make_pointer_iterator()

    def _make_cache_key(self, scan):
        return scan._data.scan_number

    def next(self):
        return next(self._producer)
