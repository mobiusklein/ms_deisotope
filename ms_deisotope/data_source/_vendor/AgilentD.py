import os
import glob
import warnings
import logging

from collections import deque

from six import string_types as basestring
from lxml import etree

try:
    log = logging.getLogger(os.path.basename(__file__))
except Exception:
    log = None
from collections import OrderedDict, defaultdict
from weakref import WeakValueDictionary

try:
    WindowsError
except NameError:
    raise ImportError("Platform Not Supported")


try:
    import comtypes
    from comtypes.client import GetModule, CreateObject
except (ImportError, NameError) as e:
    raise ImportError("Could not import comtypes")

import numpy as np

from ms_deisotope.data_source.common import (
    ScanDataSource,
    RandomAccessScanSource,
    Scan,
    ScanBunch,
    PrecursorInformation,
    ActivationInformation,
    IsolationWindow,
    InstrumentInformation,
    ComponentGroup,
    component,
    FileInformation,
    SourceFile,
    ScanAcquisitionInformation,
    ScanEventInformation,
    ScanWindow)

try:
    # Load previously built COM wrapper
    from comtypes.gen import (
        MassSpecDataReader,
        BaseCommon,
        BaseDataAccess)
    DLL_IS_LOADED = True
except (ImportError, TypeError):
    DLL_IS_LOADED = False


_default_paths = []


def _register_dll_dir(search_paths=None):
    from ms_deisotope.config import get_config
    if search_paths is None:
        search_paths = []
    global DLL_IS_LOADED
    if DLL_IS_LOADED:
        return True
    search_paths = list(search_paths)
    search_paths.extend(_default_paths)
    search_paths.extend(get_config().get('vendor_readers', {}).get('agilent-com', []))
    for dll_dir in search_paths:
        try:
            GetModule(os.path.join(dll_dir, 'MassSpecDataReader.tlb'))
            GetModule(os.path.join(dll_dir, 'BaseCommon.tlb'))
            GetModule(os.path.join(dll_dir, 'BaseDataAccess.tlb'))
            DLL_IS_LOADED = True
            return True
        except Exception:
            continue
    else:
        return False


def register_dll_dir(search_paths=None):
    if search_paths is None:
        search_paths = []
    if isinstance(search_paths, basestring):
        search_paths = [search_paths]
    loaded = _register_dll_dir(search_paths)
    if not loaded:
        log.debug("Could not resolve Agilent-related DLL")
        search_paths.extend(_default_paths)
        msg = '''
        1) The MassSpecDataReader, BaseCommon, BaseDataAccess DLLs/TLBs may not be installed and
           therefore not registered to the COM server.
        2) The MassSpecDataReader, BaseCommon, BaseDataAccess DLLs/TLBs may not be on these paths:
        %s
        ''' % ('\n'.join(search_paths))
        raise ImportError(msg)


class CaseInsensitiveDict(dict):
    def __init__(self, template=None):
        if isinstance(template, dict):
            template = {k.lower(): v for k, v in template.items()}
        dict.__init__(self, template)

    def __getitem__(self, key):
        key = key.lower()
        return dict.__getitem__(self, key)

    def __delitem__(self, key):
        return super(CaseInsensitiveDict, self).__delitem__(key.lower())

    def __setitem__(self, key, value):
        key = key.lower()
        return dict.__setitem__(self, key, value)

    def __contains__(self, key):
        return super(CaseInsensitiveDict, self).__contains__(key.lower())


device_to_component_group_map = CaseInsensitiveDict({
    "QTOF": [
        ComponentGroup("analyzer", [component("quadrupole")], 2),
        ComponentGroup("analyzer", [component("quadrupole")], 3),
        ComponentGroup("analyzer", [component("time-of-flight")], 4)
    ],
    "Quadrupole": [
        ComponentGroup("analyzer", [component("quadrupole")], 2),
    ],
    "TandemQuadrupole": [
        ComponentGroup("analyzer", [component("quadrupole")], 2),
        ComponentGroup("analyzer", [component("quadrupole")], 3),
        ComponentGroup("analyzer", [component("quadrupole")], 4)
    ],
    "IonTrap": [
        ComponentGroup("analyzer", [component("iontrap")], 2)
    ],
    "TOF": [
        ComponentGroup("analyzer", [component("time-of-flight")], 2)
    ]

})


polarity_map = {
    1: -1,
    0: 1,
    3: 0,
    2: None
}


ion_mode_map = {
    0: 'Unspecified',
    1: 'Mixed',
    2: 'EI',
    4: 'CI',
    8: 'Maldi',
    16: 'Appi',
    32: 'Apci',
    64: 'ESI',
    128: 'NanoEsi',
    512: 'MsChip',
    1024: 'ICP',
    2048: 'Jetstream'
}

ionization_map = CaseInsensitiveDict({
    "EI": component("electron ionization"),
    "CI": component("chemical ionization"),
    "ESI": component("electrospray ionization"),
    "NanoEsi": component("nanoelectrospray"),
    "Appi": component('atmospheric pressure photoionization'),
    "Apci": component("atmospheric pressure chemical ionization"),
    "Maldi": component("matrix assisted laser desorption ionization"),
    "MsChip": component("nanoelectrospray"),
    "ICP": component("plasma desorption ionization"),
    "Jetstream": component("nanoelectrospray")
})


inlet_map = CaseInsensitiveDict({
    "EI": component("direct inlet"),
    "CI": component("direct inlet"),
    "Maldi": component("particle beam"),
    "Appi": component("direct inlet"),
    "Apci": component("direct inlet"),
    "Esi": component("electrospray inlet"),
    "NanoEsi": component("nanospray inlet"),
    "MsChip": component("nanospray inlet"),
    "ICP": component("component(inductively coupled plasma"),
    "JetStream": component("nanospray inlet"),
})


peak_mode_map = {
    'profile': 0,
    'centroid': 1,
    'profilepreferred': 2,
    'centroidpreferred': 3
}

device_type_map = {
    0: 'Unknown',
    1: 'Mixed',
    2: 'Quadrupole',
    3: 'IsocraticPump',
    4: 'TOF',
    5: 'TandemQuadrupole',
    6: 'QTOF',
    10: 'FlourescenceDetector',
    11: 'ThermalConductivityDetector',
    12: 'RefractiveIndexDetector',
    13: 'MultiWavelengthDetector',
    14: 'ElectronCaptureDetector',
    15: 'VariableWavelengthDetector',
    16: 'AnalogDigitalConverter',
    17: 'EvaporativeLightScatteringDetector',
    18: 'GCDetector',
    19: 'FlameIonizationDetector',
    20: 'ALS',
    21: 'WellPlateSampler',
    22: 'MicroWellPlateSampler',
    23: 'DiodeArrayDetector',
    31: 'CANValves',
    32: 'QuaternaryPump',
    33: 'ChipCube',
    34: 'Nanopump',
    40: 'ThermostattedColumnCompartment',
    41: 'CTC',
    42: 'CapillaryPump',
    50: 'IonTrap'
}


scan_type_map = CaseInsensitiveDict({
    "Unspecified": 0,
    "All": 7951,
    "AllMS": 15,
    "AllMSN": 7936,
    "Scan": 1,
    "SelectedIon": 2,
    "HighResolutionScan": 4,
    "TotalIon": 8,
    "MultipleReaction": 256,
    "ProductIon": 512,
    "PrecursorIon": 1024,
    "NeutralLoss": 2048,
    "NeutralGain": 4096
})


PEAK_MODE = 0


def make_scan_id_string(scan_id):
    return "scanId=%s" % (scan_id,)


class AgilentDScanPtr(object):
    def __init__(self, index):
        self.index = index

    def __repr__(self):
        return "AgilentDScanPtr(%d)" % (self.index,)


class AgilentDDataInterface(ScanDataSource):
    def _get_spectrum_obj(self, scan, peak_mode=PEAK_MODE):
        index = scan.index
        spectrum = self.source.GetSpectrum_8(rowNumber=index, storageType=peak_mode)
        return spectrum

    def _get_scan_record(self, scan):
        index = scan.index
        record = self.source.GetScanRecord(index)
        return record

    def _scan_index(self, scan):
        return scan.index

    def _scan_id(self, scan):
        record = self._get_scan_record(scan)
        return make_scan_id_string(record.ScanId)

    def _scan_title(self, scan):
        return self._scan_id(scan)

    def _scan_arrays(self, scan):
        spectrum = self._get_spectrum_obj(scan)
        return (np.array(spectrum.XArray, dtype=float),
                np.array(spectrum.YArray, dtype=float))

    def _polarity(self, scan):
        record = self._get_scan_record(scan)
        polarity_enum = record.IonPolarity
        polarity = polarity_map.get(polarity_enum)
        if polarity in (0, None):
            warnings.warn("Unknown Scan Polarity: %r" % (polarity,))
        return polarity

    def _scan_time(self, scan):
        record = self._get_scan_record(scan)
        return record.retentionTime

    def _is_profile(self, scan):
        spectrum_obj = self._get_spectrum_obj(scan)
        mode = spectrum_obj.MSStorageMode
        return mode in (0, 2, 3)

    def _ms_level(self, scan):
        record = self._get_scan_record(scan)
        return record.MSLevel

    def _precursor_information(self, scan):
        if self._ms_level(scan) < 2:
            return None
        spectrum_obj = self._get_spectrum_obj(scan)
        precursor_scan_id = make_scan_id_string(spectrum_obj.ParentScanId)
        n, ions = spectrum_obj.GetPrecursorIon()
        if n < 1:
            return None
        mz = ions[0]
        charge, _ = spectrum_obj.GetPrecursorCharge()
        intensity, _ = spectrum_obj.GetPrecursorIntensity()
        return PrecursorInformation(mz, intensity, charge, precursor_scan_id, self)

    def _acquisition_information(self, scan):
        spectrum_obj = self._get_spectrum_obj(scan)
        try:
            low = spectrum_obj.MeasuredMassRange.Start
            high = spectrum_obj.MeasuredMassRange.End
        except Exception:
            arrays = self._scan_arrays(scan)
            mz_array = arrays[0]
            if len(mz_array) != 0:
                low = mz_array.min()
                high = mz_array.max()
            else:
                low = high = 0
        window = ScanWindow(low, high)
        event = ScanEventInformation(
            self._scan_time(scan),
            window_list=[window])
        return ScanAcquisitionInformation("no combination", [event])

    def _activation(self, scan):
        record = self._get_scan_record(scan)
        return ActivationInformation('cid', record.CollisionEnergy)

    def _isolation_window(self, scan):
        if self._ms_level(scan) < 2:
            return None
        spectrum_obj = self._get_spectrum_obj(scan)
        n, ions = spectrum_obj.GetPrecursorIon()
        if n < 1:
            return None
        return IsolationWindow(0, ions[0], 0)

    def _instrument_configuration(self, scan):
        return self._instrument_config[1]


class _AgilentDDirectory(object):
    @staticmethod
    def create_com_object():
        if not DLL_IS_LOADED:
            raise WindowsError("Could not locate Agilent DLLs")
        reader = CreateObject('Agilent.MassSpectrometry.DataAnalysis.MassSpecDataReader')
        return reader

    @staticmethod
    def create_com_object_filter():
        if not DLL_IS_LOADED:
            raise WindowsError("Could not locate Agilent DLLs")
        no_filter = CreateObject('Agilent.MassSpectrometry.DataAnalysis.MsdrPeakFilter')
        return no_filter

    @staticmethod
    def is_valid(path):
        if os.path.exists(path):
            if os.path.isdir(path):
                return os.path.exists(os.path.join(path, "AcqData", "Contents.xml"))
        return False


class _AgilentMethod(object):
    def __init__(self, method_parameters):
        self.parameters = list(method_parameters)

    def __getitem__(self, i):
        return self.parameters[i]

    def __len__(self):
        return len(self.parameters)

    def __iter__(self):
        return iter(self.parameters)

    def __repr__(self):
        return "_AgilentMethod(%d)" % (len(self),)

    def search_by_name(self, name):
        for param in self:
            try:
                if param['Name'].lower() == name.lower():
                    return param
            except (AttributeError, KeyError):
                continue


class _AgilentDMetadataLoader(object):
    def _has_ms1_scans(self):
        return bool(self._scan_types_flags & scan_type_map['Scan'])

    def _has_msn_scans(self):
        return bool(self._scan_types_flags & scan_type_map['ProductIon'])

    def has_msn_scans(self):
        return self._has_msn_scans()

    def has_ms1_scans(self):
        return self._has_ms1_scans()

    def file_description(self):
        fi = FileInformation(contents={}, source_files=[])
        if self._has_ms1_scans():
            fi.add_content("MS1 spectrum")
        if self._has_msn_scans():
            fi.add_content("MSn spectrum")
        basename = os.path.basename
        dirname = os.path.dirname
        file_queue = deque()
        file_queue.extend(glob.glob(os.path.join(self.dirpath, "AcqData", "*")))
        # for source_file in file_queue:
        while file_queue:
            source_file = file_queue.popleft()
            if os.path.isdir(source_file):
                file_queue.extendleft(glob.glob(os.path.join(source_file, "*")))
            else:
                sf = SourceFile(
                    basename(source_file), dirname(source_file),
                    None, *("Agilent MassHunter nativeID format", "Agilent MassHunter format"))
                sf.add_checksum("sha1")
                fi.add_file(sf, check=False)
        return fi

    def _get_instrument_info(self):
        ion_modes_flags = self.source.MSScanFileInformation.IonModes
        ionization = []
        for bit, label in ion_mode_map.items():
            if ion_modes_flags & bit:
                ionization.append(label)
        configs = []
        i = 1
        for ionizer in ionization:
            groups = [ComponentGroup("source", [ionization_map[ionizer], inlet_map[ionizer]], 1)]
            groups.extend(device_to_component_group_map[self.device])
            config = InstrumentInformation(i, groups)
            i += 1
            configs.append(config)
        self._instrument_config = {
            c.id: c for c in configs
        }
        return configs

    def instrument_configuration(self):
        return sorted(self._instrument_config.values(), key=lambda x: x.id)

    def data_processing(self):
        return []

    def _acquisition_method_xml_path(self):
        return os.path.join(self.dirpath, "AcqData", "AcqMethod.xml")

    def _parse_method_xml(self):
        try:
            path = self._acquisition_method_xml_path()
            tree = etree.parse(path)
            nsmap = {"ns": "http://tempuri.org/DataFileReport.xsd"}
            elt = tree.find(".//ns:SCICDevicesXml", namespaces=nsmap)
            method_xml = etree.fromstring(elt.text)
        except (IOError, OSError, ValueError, TypeError) as e:
            print(e)
            self._method = []
            return self._method
        method = list()
        for section in method_xml.iterfind(".//SectionInfo"):
            section_dict = {}
            for child in section:
                name = child.tag
                value = child.text
                section_dict[name] = value
            method.append(section_dict)
        method = _AgilentMethod(method)
        self._method = method
        return method


_ADM = _AgilentDMetadataLoader
_ADD = _AgilentDDirectory


class AgilentDLoader(AgilentDDataInterface, _ADD, RandomAccessScanSource, _ADM):

    def __init__(self, dirpath, **kwargs):
        self.dirpath = dirpath
        self.dirpath = os.path.abspath(self.dirpath)
        self.dirpath = os.path.normpath(self.dirpath)

        self.source = self.create_com_object()
        self.filter = self.create_com_object_filter()

        try:
            self.source.OpenDataFile(self.dirpath)
        except comtypes.COMError as err:
            raise IOError(str(err))
        self._TIC = self.source.GetTIC()
        self.device = self._TIC.DeviceName
        self._n_spectra = self._TIC.TotalDataPoints
        self._scan_types_flags = self.source.MSScanFileInformation.ScanTypes

        self._producer = self._scan_group_iterator()
        self.initialize_scan_cache()
        self._index = self._pack_index()
        self._get_instrument_info()

    def __reduce__(self):
        return self.__class__, (self.dirpath,)

    @property
    def index(self):
        return self._index

    def __len__(self):
        return len(self.index)

    def __repr__(self):
        return "AgilentDLoader(%r)" % (self.dirpath)

    def reset(self):
        self.make_iterator(None)
        self.initialize_scan_cache()

    def close(self):
        # seems to make attempting to re-open the same datafile cause a segfault
        # self.source.CloseDataFile()
        self._dispose()

    def _pack_index(self):
        index = OrderedDict()
        for sn in range(self._n_spectra):
            rec = self._get_scan_record(AgilentDScanPtr(sn))
            index[make_scan_id_string(rec.ScanId)] = sn
        return index

    def _make_pointer_iterator(self, start_index=None, start_time=None):
        iterator = self._make_scan_index_producer(start_index, start_time)
        for i in iterator:
            yield AgilentDScanPtr(i)

    def _make_default_iterator(self):
        return self._make_pointer_iterator()

    def _make_scan_index_producer(self, start_index=None, start_time=None):
        if start_index is not None:
            return range(start_index, self._n_spectra)
        elif start_time is not None:
            start_index = self._source.ScanNumFromRT(start_time)
            while start_index != 0:
                scan = self.get_scan_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            return range(start_index, self._n_spectra)
        else:
            return range(0, self._n_spectra)

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
        index = self._index[scan_id]
        return self.get_scan_by_index(index)

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
            package = AgilentDScanPtr(scan_number)
            scan = Scan(package, self)
            self._cache_scan(scan)
            return scan

    def get_scan_by_time(self, time):
        time_array = self._TIC.XArray
        lo = 0
        hi = self._n_spectra

        if time == float('inf'):
            return self.get_scan_by_index(len(self) - 1)

        best_match = None
        best_error = float('inf')

        while hi != lo:
            mid = (hi + lo) // 2
            scan_time = time_array[mid]
            err = abs(scan_time - time)
            if err < best_error:
                best_error = err
                best_match = mid
            if scan_time == time:
                return self.get_scan_by_index(mid)
            elif (hi - lo) == 1:
                return self.get_scan_by_index(best_match)
            elif scan_time > time:
                hi = mid
            else:
                lo = mid

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
            scan_number = self.get_scan_by_id(scan_id).index
        elif index is not None:
            scan_number = int(index)
        elif rt is not None:
            scan_number = self.get_scan_by_time(rt).index
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

    def _make_cache_key(self, scan):
        return scan._data.index

    def _single_scan_iterator(self, iterator=None, mode=None):
        if iterator is None:
            iterator = self._make_scan_index_producer()
        for ix in iterator:
            packed = self.get_scan_by_index(ix)
            self._cache_scan(packed)
            yield packed

    def _scan_group_iterator(self, iterator=None, mode=None):
        if iterator is None:
            iterator = self._make_scan_index_producer()

        precursor_scan = None
        product_scans = []

        current_level = 1

        for ix in iterator:
            packed = self.get_scan_by_index(ix)
            self._cache_scan(packed)
            if packed.ms_level > 1:
                # inceasing ms level
                if current_level < packed.ms_level:
                    current_level = packed.ms_level
                # decreasing ms level
                elif current_level > packed.ms_level:
                    current_level = packed.ms_level
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
