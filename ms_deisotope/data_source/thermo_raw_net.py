import sys
import os
import numpy as np

from collections import OrderedDict


from ms_deisotope.data_source.common import (
    PrecursorInformation, ChargeNotProvided, Scan,
    ActivationInformation, MultipleActivationInformation,
    IsolationWindow, ScanDataSource, ScanEventInformation,
    ScanAcquisitionInformation, ScanWindow, RandomAccessScanSource)

from ms_deisotope.data_source.metadata.activation import (
    supplemental_term_map, dissociation_methods_map)

from ms_deisotope.data_source._thermo_helper import (
    _InstrumentMethod, ThermoRawScanPtr, FilterString,
    _make_id, _id_template, _RawFileMetadataLoader, analyzer_map)


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


def is_thermo_raw_file(path):
    if not _test_dll_loaded():
        try:
            register_dll()
        except ImportError:
            return False
    try:
        source = _RawFileReader.RawFileReaderAdapter.FileFactory(path)
        source.SelectInstrument(Business.Device.MS, 1)
        return True
    except (NullReferenceException):
        return False


def infer_reader(path):
    if is_thermo_raw_file(path):
        return ThermoRawLoader
    raise ValueError("Not Thermo Raw File")


def determine_if_available():
    try:
        return _register_dll([_DEFAULT_DLL_PATH])
    except (OSError, ImportError):
        return False


def _register_dll(search_paths=None):
    if search_paths is None:
        search_paths = []
    global _RawFileReader, Business, clr, NullReferenceException
    if _test_dll_loaded():
        return True
    try:
        import clr
        from System import NullReferenceException
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
            import ThermoFisher.CommonCore.Data.Business as Business
            import ThermoFisher.CommonCore.RawFileReader as _RawFileReader
        except ImportError:
            continue
    return _test_dll_loaded()


def register_dll(search_paths=None):
    if search_paths is None:
        search_paths = []
    search_paths = list(search_paths)
    search_paths.append(_DEFAULT_DLL_PATH)
    loaded = _register_dll(search_paths)
    if not loaded:
        msg = '''The ThermoFisher.CommonCore libraries could not be located and loaded.'''
        raise ImportError(msg)


def _test_dll_loaded():
    return _RawFileReader is not None


class RawReaderInterface(ScanDataSource):

    def _scan_arrays(self, scan):
        scan_number = scan.scan_number + 1
        stats = self._source.GetScanStatsForScanNumber(scan_number)
        segscan = self._source.GetSegmentedScanFromScanNumber(scan_number, stats)
        mzs = np.array(list(segscan.Positions), dtype=np.float64)
        inten = np.array(list(segscan.Intensities), dtype=np.float64)
        return mzs, inten

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
        scan_number = scan.scan_number
        return FilterString(self._source.GetFilterForScanNumber(scan_number + 1).Filter)

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
        filt = self._source.GetFilterForScanNumber(scan_number + 1)
        seq_index = filt.MSOrder - 2
        width = filt.GetIsolationWidth(seq_index)
        offset = filt.GetIsolationWidthOffset(seq_index)
        precursor_mz = filt.GetMass(seq_index)
        return IsolationWindow(width, precursor_mz + offset, width)

    def _trailer_values(self, scan):
        scan_number = scan.scan_number
        trailers = self._source.GetTrailerExtraInformation(scan_number + 1)
        return OrderedDict(zip(map(lambda x: x.strip(":"), trailers.Labels), trailers.Values))

    def _infer_precursor_scan_number(self, scan):
        precursor_scan_number = None
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
        return precursor_scan_number

    def _precursor_information(self, scan):
        scan_number = scan.scan_number
        filt = self._source.GetFilterForScanNumber(scan_number + 1)
        precursor_mz = filt.GetMass(filt.MSOrder - 2)
        trailers = self._trailer_values(scan)
        charge = int(trailers.get("Charge State", 0))
        if charge == 0:
            charge = ChargeNotProvided
        inten = 0
        precursor_scan_number = None
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
        if precursor_scan_number is not None:
            precursor_scan_id = self.get_scan_by_index(precursor_scan_number).id
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
        traits = {
            'preset scan configuration': event,
            'filter string': fline,
        }
        event = ScanEventInformation(
            self._scan_time(scan),
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
            "filter_string": fline,
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
        hcd_ev = trailer_extras.get('HCD Energy eV')
        if hcd_ev is not None and hcd_ev > 0:
            annots['[Thermo Trailer Extra]HCD Energy eV'] = hcd_ev
        hcd_energies = trailer_extras.get('HCD Energy')
        if hcd_energies is not None and hcd_energies:
            annots['[Thermo Trailer Extra]HCD Energy'] = hcd_energies
        return annots


class ThermoRawLoader(RawReaderInterface, RandomAccessScanSource, _RawFileMetadataLoader):
    def __init__(self, source_file, _load_metadata=True, **kwargs):
        if not _test_dll_loaded():
            register_dll()
        self._source = _RawFileReader.RawFileReaderAdapter.FileFactory(source_file)
        self._source.SelectInstrument(Business.Device.MS, 1)
        self.source_file = source_file
        self._producer = None
        self._scan_type_index = dict()
        self.make_iterator()
        self.initialize_scan_cache()
        self._first_scan_time = self.get_scan_by_index(0).scan_time
        self._last_scan_time = self.get_scan_by_index(self._source.RunHeaderEx.LastSpectrum - 1).scan_time
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
        return self.__class__, (self.source_file, False)

    @property
    def index(self):
        return self._index

    def __len__(self):
        return len(self.index)

    def __repr__(self):
        return "ThermoRawLoader(%r)" % (self.source_file)

    def close(self):
        if self._source is not None:
            self._source.Close()
            self._source = None

    def __del__(self):
        self.close()

    def reset(self):
        self.make_iterator(None)
        self.initialize_scan_cache()

    def _pack_index(self):
        index = OrderedDict()
        for sn in range(self._source.RunHeaderEx.FirstSpectrum - 1, self._source.RunHeaderEx.LastSpectrum):
            index[_make_id(sn)] = sn
        return index

    def _parse_method(self):
        try:
            return _InstrumentMethod(self._source.GetInstrumentMethod(0))
        except NullReferenceException:
            return _InstrumentMethod('')

    def _scan_time_to_scan_number(self, scan_time):
        scan_number = self._source.ScanNumberFromRetentionTime(scan_time) - 1
        return scan_number

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
                raise KeyError(index)
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
            return range(start_index, self._source.RunHeaderEx.LastSpectrum - 1)
        elif start_time is not None:
            start_index = self._scan_time_to_scan_number(start_time)
            while start_index != 0:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            return range(start_index, self._source.RunHeaderEx.LastSpectrum - 1)
        else:
            return range(0, self._source.RunHeaderEx.LastSpectrum - 1)

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
