from collections import namedtuple
import abc

from ms_peak_picker import pick_peaks
from ..averagine import neutral_mass, mass_charge_ratio
from ..utils import Constant

ScanBunch = namedtuple("ScanBunch", ["precursor", "products"])


ChargeNotProvided = Constant("ChargeNotProvided")
IsolationWindowNotProvided = Constant("IsolationWindowNotProvided")


class ScanDataSourceBase(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def scan_arrays(self, scan):
        raise NotImplementedError()

    @abc.abstractmethod
    def precursor_information(self, scan):
        raise NotImplementedError()

    @abc.abstractmethod
    def scan_title(self, scan):
        raise NotImplementedError()

    @abc.abstractmethod
    def scan_id(self, scan):
        raise NotImplementedError()

    @abc.abstractmethod
    def scan_index(self, scan):
        raise NotImplementedError()

    @abc.abstractmethod
    def ms_level(self, scan):
        raise NotImplementedError()

    @abc.abstractmethod
    def scan_time(self, scan):
        raise NotImplementedError()


class ScanIteratorBase(ScanDataSourceBase):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def next(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def get_scan_by_id(self, scan_id):
        raise NotImplementedError()

    def _make_scan(self, data):
        return Scan(data, self)

    def __next__(self):
        return self.next()

    def __iter__(self):
        return self


class Scan(object):

    def __init__(self, data, source, peak_set=None, deconvoluted_peak_set=None, product_scans=None):
        if product_scans is None:
            product_scans = []
        self.source = source
        self.peak_set = peak_set
        self.deconvoluted_peak_set = deconvoluted_peak_set

        self._data = data

        self._arrays = None
        self._id = None
        self._title = None
        self._ms_level = None
        self._scan_time = None
        self._precursor_information = None
        self._index = None

        self.product_scans = product_scans

    def __getitem__(self, key):
        return self._data[key]

    def __iter__(self):
        return iter(self.peak_set)

    def has_peak(self, *args, **kwargs):
        return self.peak_set.has_peak(*args, **kwargs)

    @property
    def ms_level(self):
        if self._ms_level is None:
            self._ms_level = self.source.ms_level(self._data)
        return self._ms_level

    @property
    def scan_time(self):
        if self._scan_time is None:
            self._scan_time = self.source.scan_time(self._data)
        return self._scan_time

    @property
    def arrays(self):
        return self.source.scan_arrays(self._data)

    @property
    def title(self):
        return self.source.scan_title(self._data)

    @property
    def id(self):
        return self.source.scan_id(self._data)

    @property
    def index(self):
        if self._index is None:
            self._index = self.source.scan_index(self._data)
        return self._index

    @property
    def precursor_information(self):
        if self.ms_level < 2:
            return None
        if self._precursor_information is None:
            self._precursor_information = self.source.precursor_information(self._data)
        return self._precursor_information

    def __repr__(self):
        return "Scan(%s %s)" % (
            self.id,
            self.precursor_information if self.precursor_information else '')

    def pick_peaks(self, *args, **kwargs):
        mzs, intensities = self.arrays
        self.peak_set = pick_peaks(mzs, intensities, *args, **kwargs)
        return self

    def pack(self):
        return ProcessedScan(
            self.id, self.title, self.precursor_information,
            self.ms_level, self.scan_time,
            self.deconvoluted_peak_set)


class PrecursorInformation(object):
    def __init__(self, mz, intensity, charge, precursor_scan_id=None, _source=None,
                 extracted_neutral_mass=0, extracted_charge=0, extracted_peak_height=0,
                 peak=None, extracted_peak=None):
        self.mz = mz
        self.intensity = intensity
        self.charge = charge
        self.precursor_scan_id = precursor_scan_id
        self._source = _source
        self.extracted_neutral_mass = extracted_neutral_mass
        self.extracted_charge = extracted_charge
        self.extracted_peak_height = extracted_peak_height
        self.peak = peak
        self.extracted_peak = extracted_peak

    def __repr__(self):
        return "PrecursorInformation(mz=%0.4f/%0.4f, intensity=%0.4f/%0.4f, charge=%d/%d, scan_id=%s)" % (
            self.mz,
            mass_charge_ratio(self.extracted_neutral_mass, self.extracted_charge if self.extracted_charge != 0 else 1)
            if self.extracted_neutral_mass != 0. else 0.,
            self.intensity, self.extracted_peak_height or 0., self.charge,
            self.extracted_charge or 0., self.precursor_scan_id)

    @property
    def neutral_mass(self):
        return neutral_mass(self.mz, self.charge)

    @property
    def extracted_mz(self):
        return mass_charge_ratio(self.extracted_neutral_mass, self.extracted_charge)

    @property
    def precursor(self):
        return self._source.get_scan_by_id(self.precursor_scan_id)


class ProcessedScan(object):
    def __init__(self, id, title, precursor_information, ms_level, scan_time, peaks):
        self.id = id
        self.title = title
        self.precursor_information = precursor_information
        self.ms_level = ms_level
        self.scan_time = scan_time
        self.peaks = peaks

    def __iter__(self):
        return iter(self.peaks)

    def __getitem__(self, index):
        return self.peaks[index]

    def __repr__(self):
        return "ProcessedScan(id=%s, ms_level=%d, %d peaks)" % (self.id, self.ms_level, len(self.peaks))
