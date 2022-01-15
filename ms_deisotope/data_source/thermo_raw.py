'''Thermo RAW file reading implementation using the Windows-specific COM
interface provided by MSFileReader.dll library provided by Thermo.

This module provides :class:`ThermoRawLoader`, a :class:`~.RandomAccessScanSource`
implementation.

Depends upon :mod:`comtypes` to interact with the operating system's
COM API and to register and create COM wrapper modules on demand.

The public interface of this module should be identical to
:mod:`ms_deisotope.data_source.thermo_raw_net`.

.. note::
    This interface was largely based upon the APIs that ProteoWizard used, both
    in order to understand how the Thermo libraries really worked, and to maintain
    parity with it.

'''
# pragma: no cover

import warnings
from collections import OrderedDict

import logging

import numpy as np

from pyteomics.auxiliary import unitfloat

from ms_deisotope.data_source.common import (
    ScanDataSource, RandomAccessScanSource,
    Scan, PrecursorInformation, ChargeNotProvided,
    ActivationInformation, IsolationWindow,
    ScanAcquisitionInformation, ScanEventInformation, ScanWindow,
    MultipleActivationInformation)

from ms_deisotope.data_source._thermo_helper import (
    _RawFileMetadataLoader, analyzer_map,
    FilterString, _id_template, _InstrumentMethod,
    _make_id, ThermoRawScanPtr)

from ms_deisotope.data_source.metadata.activation import (
    supplemental_term_map, dissociation_methods_map)
from ms_deisotope.data_source.metadata.sample import Sample
from ms_deisotope.data_source.metadata.scan_traits import FAIMS_compensation_voltage

try:
    from ms_deisotope.data_source._vendor.MSFileReader import (  # pylint: disable=unused-import
        ThermoRawfile as _ThermoRawFileAPI, register_dll,
        log as _api_logger, safearray_as_ndarray)

    comtypes_logger = logging.getLogger("comtypes")
    comtypes_logger.setLevel("INFO")
    _api_logger.setLevel("INFO")

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
        try:
            _ThermoRawFileAPI(path)
            return True
        except (WindowsError, IOError, ImportError):
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
        '''Checks whether or not the COM-based Thermo
        RAW file reading feature is available.

        This is done by attempting to instantiate the
        COM-provided object, which queries the Windows
        registry for the MSFileReader.dll.

        Returns
        -------
        :class:`bool`:
            Whether or not the feature is enabled.
        '''
        try:
            _ThermoRawFileAPI.create_com_object()
            return True
        except ImportError:
            return False
except ImportError as e:  # pragma: no cover
    message = str(e)

    def is_thermo_raw_file(*args, **kwargs):
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
        return False


    def infer_reader(*args, **kwargs):
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
        raise ValueError(message)


    def register_dll(*args, **kwargs):
        '''Register the location of Thermo's MSFileReader.dll
        with the COM interop system.
        '''
        warnings.warn("no-op: %s" % (message,))
        return False


    def determine_if_available():
        '''Checks whether or not the COM-based Thermo
        RAW file reading feature is available.

        This is done by attempting to instantiate the
        COM-provided object, which queries the Windows
        registry for the MSFileReader.dll.

        Returns
        -------
        :class:`bool`:
            Whether or not the feature is enabled.
        '''
        warnings.warn("no-op: %s" % (message,))
        return False

try:
    range = xrange
except NameError:
    pass


class ThermoRawDataInterface(ScanDataSource):
    ''':class:`~.ScanDataSource` implementation for Thermo's MSFileReader API.
    '''
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
        filter_string = self._filter_string(scan)
        return filter_string.data['polarity']

    def _filter_string(self, scan):
        if scan.filter_string is None:
            scan.filter_string = FilterString(self._source.GetFilterForScanNum(scan.scan_number))
        return scan.filter_string

    def _scan_title(self, scan):
        return "%s %r" % (self._scan_id(scan), self._filter_string(scan))

    def _scan_arrays(self, scan):
        with safearray_as_ndarray:
            arrays, _ = self._source.GetMassListFromScanNum(
                scan.scan_number)
        mz, intensity = arrays
        if isinstance(mz, np.ndarray):
            return mz.copy(), intensity.copy()
        return np.array(mz), np.array(intensity)

    def _precursor_information(self, scan):
        if self._ms_level(scan) == 1:
            return None
        scan_number = scan.scan_number
        pinfo_struct = self._source.GetPrecursorInfoFromScanNum(scan_number)
        trailer = self._trailer_values(scan)
        precursor_scan_number = trailer.get('Master Scan Number')
        if precursor_scan_number is not None:
            precursor_scan_number = int(precursor_scan_number) - 1
            if self.get_scan_by_index(precursor_scan_number).ms_level >= self._ms_level(scan):
                precursor_scan_number = None

        labels, _, _ = self._source.GetAllMSOrderData(scan_number)
        if pinfo_struct:
            # this struct field is unreliable and may fall outside the
            # isolation window
            mz = pinfo_struct.monoIsoMass
            charge = pinfo_struct.chargeState
            intensity = float(labels.intensity[0])
            # this struct field is unreliable, and simple to infer
            # precursor_scan_number = pinfo_struct.scanNumber + 1
        else:
            mz = labels.mass[0]
            intensity = float(labels.intensity[0])
            charge = labels.charge[0]
        if not charge:
            charge = ChargeNotProvided
        _mz = trailer.get('Monoisotopic M/Z', 0.0)
        # prefer the trailer m/z if available?
        if _mz > 0:
            mz = _mz

        # imitate proteowizard's firmware bug correction
        isolation_window = self._isolation_window(scan)
        if (isolation_window.upper + isolation_window.lower) / 2 <= 2.0:
            if (isolation_window.target - 3.0 > mz) or (isolation_window.target + 2.5 < mz):
                mz = isolation_window.target
        elif mz not in isolation_window:
            mz = isolation_window.target
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

    def _isolation_window(self, scan):
        ms_level = self._ms_level(scan)
        if ms_level == 1:
            return None
        isolation_width = 0
        trailer = self._trailer_values(scan)
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
        fline = self._filter_string(scan)
        try:
            confid = self._analyzer_to_configuration_index[analyzer_map[fline.data.get("analyzer")]]
            return self._instrument_config[confid]
        except KeyError:
            return None

    def _trailer_values(self, scan):
        if scan.trailer_values is not None:
            return scan.trailer_values
        trailer_extras = self._source.GetTrailerExtraForScanNum(scan.scan_number)
        scan.trailer_values = trailer_extras
        return trailer_extras

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
            injection_time=unitfloat(trailer_extras.get(
                'Ion Injection Time (ms)', 0.0), 'millisecond'),
            window_list=[ScanWindow(
                fline.get("scan_window")[0], fline.get("scan_window")[1])],
            traits=traits)
        return ScanAcquisitionInformation("no combination", [event])

    def _annotations(self, scan):
        fline = self._filter_string(scan)
        trailer_extras = self._trailer_values(scan)
        annots = {
            "filter string": fline,
        }
        microscans = trailer_extras.get("Micro Scan Count")
        if microscans is not None:
            annots['[Thermo Trailer Extra]Micro Scan Count'] = microscans
        scan_segment = trailer_extras.get("Scan Segment")
        if scan_segment is not None:
            annots['[Thermo Trailer Extra]Scan Segment'] = scan_segment
        scan_event = trailer_extras.get("Scan Event")
        if scan_event is not None:
            annots['[Thermo Trailer Extra]Scan Event'] = scan_event
        mono_mz = trailer_extras.get("Monoisotopic M/Z")
        if mono_mz is not None and mono_mz > 0:
            annots['[Thermo Trailer Extra]Monoisotopic M/Z'] = mono_mz
        hcd_ev = trailer_extras.get('HCD Energy eV')
        if hcd_ev is not None and hcd_ev > 0:
            annots['[Thermo Trailer Extra]HCD Energy eV'] = hcd_ev
        return annots


class ThermoRawLoader(ThermoRawDataInterface, RandomAccessScanSource, _RawFileMetadataLoader):
    '''Reads scans from Thermo Fisher RAW files directly. Provides both iterative and
    random access.
    '''

    def __init__(self, source_file, _load_metadata=True, **kwargs):
        self.source_file = source_file
        self._source = _ThermoRawFileAPI(self.source_file)
        self._producer = None
        self._scan_type_index = dict()
        self.make_iterator()
        self.initialize_scan_cache()
        self._index = self._pack_index()
        self._first_scan_time = self.get_scan_by_index(0).scan_time
        self._last_scan_time = self.get_scan_by_id(self._source.LastSpectrumNumber).scan_time
        if _load_metadata:
            self._method = self._parse_method()
            self._build_scan_type_index()
            self._get_instrument_info()

    def _get_instrument_model_name(self):
        return self._source.GetInstModel()

    def _get_instrument_serial_number(self):
        return self._source.GetInstSerialNumber()

    def samples(self):
        """Describe the sample(s) used to generate the mass spectrometry
        data contained in this file.

        Returns
        -------
        :class:`list` of :class:`~.Sample`
        """
        result = []
        si = self._source
        sample = Sample(si.GetSeqRowSampleID() or 'sample_1')
        sample.name = si.GetSeqRowSampleName() or si.GetSeqRowSampleID()
        if si.GetSeqRowSampleVolume():
            sample.parameters['sample volume'] = si.GetSeqRowSampleVolume()
        if si.GetSeqRowSampleWeight():
            sample.parameters['sample mass'] = si.GetSeqRowSampleWeight()
        if si.GetSeqRowVial():
            sample.parameters['sample vial'] = si.GetSeqRowVial()
        if si.GetSeqRowBarcode():
            sample.parameters['sample barcode'] = si.GetSeqRowBarcode()
        result.append(sample)
        return result

    def _parse_method(self):
        return _InstrumentMethod(self._source.GetInstMethod())

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

    def __repr__(self):
        return "ThermoRawLoader(%r)" % (self.source_file)

    def _pack_index(self):
        index = OrderedDict()
        for sn in range(1, self._source.NumSpectra + 1):
            index[_make_id(sn)] = sn
        return index

    def __len__(self):
        return len(self.index)

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

    def _scan_time_to_scan_number(self, scan_time):
        if scan_time < self._first_scan_time:
            scan_time = self._first_scan_time
        elif scan_time > self._last_scan_time:
            scan_time = self._last_scan_time
        scan_number = self._source.ScanNumFromRT(scan_time)
        return scan_number

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
            if not package.validate(self):
                raise KeyError(str(scan_id))
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
            if not package.validate(self):
                raise IndexError(index)
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
            return range(start_index + 1, self._source.NumSpectra + 1)
        elif start_time is not None:
            start_index = self._scan_time_to_scan_number(start_time)
            while start_index != 1:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            return range(start_index, self._source.NumSpectra + 1)
        else:
            return range(1, self._source.NumSpectra + 1)

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
