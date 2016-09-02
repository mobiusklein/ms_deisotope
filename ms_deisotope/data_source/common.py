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

    @abc.abstractmethod
    def is_profile(self, scan):
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


class DetachedAccessError(Exception):
    pass


class DataAccessProxy(object):
    def attach(self, source):
        self.source = source

    def detach(self):
        self.source = None

    def __getstate__(self):
        return ()

    def __setstate__(self, state):
        self.source = None

    def raise_if_detached(self):
        if self.source is None:
            raise DetachedAccessError("Cannot perform operation. Instance is detached.")

    def get_scan_by_id(self, scan_id):
        self.raise_if_detached()
        return self.source.get_scan_by_id(scan_id)


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
        self._is_profile = None
        self._polarity = None

        self.product_scans = product_scans

    def _load(self):
        self.arrays
        self.id
        self.title
        self.ms_level
        self.scan_time
        self.index
        self.precursor_information

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
    def is_profile(self):
        if self._is_profile is None:
            self._is_profile = self.source.is_profile(self._data)
        return self._is_profile

    @property
    def polarity(self):
        if self._polarity is None:
            self._polarity = self.source.polarity(self._data)
        return self._polarity

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
        return "Scan(%s, index=%d, time=%s%s)" % (
            self.id, self.index, self.scan_time,
            " " + str(self.precursor_information) if self.precursor_information else '')

    def pick_peaks(self, *args, **kwargs):
        mzs, intensities = self.arrays
        self.peak_set = pick_peaks(mzs, intensities, *args, **kwargs)
        return self

    def pack(self):
        precursor_info = self.precursor_information
        return ProcessedScan(
            self.id, self.title, precursor_info,
            self.ms_level, self.scan_time, self.index,
            self.peak_set.pack(),
            self.deconvoluted_peak_set)


class PrecursorInformation(object):
    def __init__(self, mz, intensity, charge, precursor_scan_id=None, source=None,
                 extracted_neutral_mass=0, extracted_charge=0, extracted_intensity=0,
                 peak=None, extracted_peak=None, defaulted=False, orphan=False):
        self.mz = mz
        self.intensity = intensity
        self.charge = charge
        self.precursor_scan_id = precursor_scan_id
        self.source = source
        self.extracted_neutral_mass = extracted_neutral_mass
        self.extracted_charge = extracted_charge
        self.extracted_intensity = extracted_intensity
        self.peak = peak
        self.extracted_peak = extracted_peak
        self.defaulted = defaulted
        self.orphan = orphan

    def __repr__(self):
        return "PrecursorInformation(mz=%0.4f/%0.4f, intensity=%0.4f/%0.4f, charge=%d/%d, scan_id=%s)" % (
            self.mz,
            mass_charge_ratio(self.extracted_neutral_mass, self.extracted_charge if self.extracted_charge != 0 else 1)
            if self.extracted_neutral_mass != 0. else 0.,
            self.intensity, self.extracted_intensity or 0., self.charge,
            self.extracted_charge or 0., self.precursor_scan_id)

    def __getstate__(self):
        return (self.mz, self.intensity, self.charge, self.precursor_scan_id, None, self.extracted_neutral_mass,
                self.extracted_charge, self.extracted_intensity, self.peak, self.extracted_peak)

    def __setstate__(self, state):
        (self.mz, self.intensity, self.charge, self.precursor_scan_id, self.source, self.extracted_neutral_mass,
         self.extracted_charge, self.extracted_intensity, self.peak, self.extracted_peak) = state

    def extract(self, peak, override_charge=None):
        self.extracted_neutral_mass = peak.neutral_mass
        self.extracted_charge = int(peak.charge) if override_charge is None else override_charge
        self.extracted_intensity = peak.intensity
        self.extracted_peak = peak

    def default(self):
        self.extracted_neutral_mass = self.neutral_mass
        self.extracted_charge = int(self.charge)
        self.extracted_intensity = self.intensity
        self.defaulted = True

    @property
    def neutral_mass(self):
        return neutral_mass(self.mz, self.charge)

    @property
    def extracted_mz(self):
        return mass_charge_ratio(self.extracted_neutral_mass, self.extracted_charge)

    @property
    def precursor(self):
        return self.source.get_scan_by_id(self.precursor_scan_id)


class ProcessedScan(object):
    def __init__(self, id, title, precursor_information, ms_level, scan_time, index, peak_set, deconvoluted_peak_set):
        self.id = id
        self.title = title
        self.precursor_information = precursor_information
        self.ms_level = ms_level
        self.scan_time = scan_time
        self.index = index
        self.peak_set = peak_set
        self.deconvoluted_peak_set = deconvoluted_peak_set

    def __iter__(self):
        return iter(self.deconvoluted_peak_set)

    def __getitem__(self, index):
        return self.deconvoluted_peak_set[index]

    def __repr__(self):
        return "ProcessedScan(id=%s, ms_level=%d, %d peaks)" % (self.id, self.ms_level, len(self.deconvoluted_peak_set))
