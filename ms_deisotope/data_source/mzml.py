import numpy as np
from pyteomics import mzml
from .common import (
    PrecursorInformation, ScanDataSource,
    ChargeNotProvided, ActivationInformation)
from weakref import WeakValueDictionary
from .xml_reader import XMLReaderBase


class MzMLDataInterface(ScanDataSource):
    """Provides implementations of all of the methods needed to implement the
    :class:`ScanDataSource` for mzML files. Not intended for direct instantiation.
    """
    def _scan_arrays(self, scan):
        try:
            return scan['m/z array'], scan["intensity array"]
        except KeyError:
            return np.array([]), np.array([])

    def _precursor_information(self, scan):
        pinfo_dict = scan["precursorList"]['precursor'][0]["selectedIonList"]['selectedIon'][0]
        try:
            precursor_scan_id = scan["precursorList"]['precursor'][0]['spectrumRef']
        except KeyError:
            precursor_scan_id = None
        pinfo = PrecursorInformation(
            mz=pinfo_dict['selected ion m/z'],
            intensity=pinfo_dict.get('peak intensity', 0.0),
            charge=pinfo_dict.get('charge state', ChargeNotProvided),
            precursor_scan_id=precursor_scan_id,
            source=self,
            product_scan_id=self._scan_id(scan))
        return pinfo

    # Optional
    def _scan_title(self, scan):
        try:
            return scan["spectrum title"]
        except KeyError:
            return scan["id"]

    def _scan_id(self, scan):
        return scan["id"]

    def _scan_index(self, scan):
        return scan['index']

    def _ms_level(self, scan):
        return scan['ms level']

    def _scan_time(self, scan):
        try:
            return scan['scanList']['scan'][0]['scan start time']
        except KeyError:
            return None

    def _is_profile(self, scan):
        return "profile spectrum" in scan

    def _polarity(self, scan):
        if "positive scan" in scan:
            return 1
        elif "negative scan" in scan:
            return -1

    def _activation(self, scan):
        try:
            struct = dict(scan['precursorList']['precursor'][0]['activation'])
            activation = None
            for key in tuple(struct):
                if key in ActivationInformation.dissociation_methods:
                    activation = key
                    struct.pop(key)
                    break
            else:
                activation = "unknown dissociation method"
            energy = struct.pop("collision energy", -1)
            return ActivationInformation(activation, energy, struct)
        except KeyError:
            return None


def _find_section(source, section):
    value = next(source.iterfind(section))
    source.reset()
    return value


class MzMLLoader(MzMLDataInterface, XMLReaderBase):
    """Reads scans from PSI-HUPO mzML XML files. Provides both iterative and
    random access.

    Attributes
    ----------
    source_file: str
        Path to file to read from.
    source: pyteomics.mzml.MzML
        Underlying scan data source
    """
    __data_interface__ = MzMLDataInterface

    def __init__(self, source_file, use_index=True):
        self.source_file = source_file
        self._source = mzml.MzML(source_file, read_schema=True, iterative=True, use_index=use_index)
        self._producer = self._scan_group_iterator()
        self._scan_cache = WeakValueDictionary()
        self._use_index = use_index

    def __reduce__(self):
        return self.__class__, (self.source_file, self._use_index)

    def file_description(self):
        return _find_section(self._source, "fileDescription")

    def instrument_configuration(self):
        return _find_section(self._source, "instrumentConfigurationList")

    def data_processing(self):
        return _find_section(self._source, "dataProcessingList")

    def samples(self):
        return _find_section(self._source, "sampleList")

    def _validate(self, scan):
        return "m/z array" in scan._data

    def _yield_from_index(self, scan_source, start):
        offset_provider = scan_source._offset_index.offsets
        keys = offset_provider.keys()
        if start is not None:
            if isinstance(start, basestring):
                start = keys.index(start)
            elif isinstance(start, int):
                start = start
            else:
                raise TypeError("Cannot start from object %r" % start)
        else:
            start = 0
        for key in keys[start:]:
            yield scan_source.get_by_id(key)


# PyteomicsMzMLLoader = MzMLLoader
