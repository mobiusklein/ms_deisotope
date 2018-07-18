import json

from collections import OrderedDict

from ms_deisotope.data_source.common import PrecursorInformation, ChargeNotProvided
from ms_deisotope.utils import Base


class MSRecordBase(Base):
    def __init__(self, scan_time, drift_time=None, **kwargs):
        self.scan_time = scan_time
        self.drift_time = drift_time
        self._extra_keys = kwargs.keys()
        self.__dict__.update(kwargs)

    def __getitem__(self, key):
        return getattr(self, key)

    def to_dict(self):
        package = {
            "scan_time": self.scan_time,
        }
        if self.drift_time is not None:
            package['drift_time'] = self.drift_time
        for key in self._extra_keys:
            package[key] = self[key]
        return package

    def get(self, key, default=None):
        try:
            return self[key]
        except Exception:
            return default


class MS1Record(MSRecordBase):
    def __init__(self, scan_time=None, product_scan_ids=None, msms_peaks=None, drift_time=None, **kwargs):
        super(MS1Record, self).__init__(scan_time, drift_time, **kwargs)
        self.product_scan_ids = product_scan_ids or []
        self.msms_peaks = msms_peaks or []

    def to_dict(self):
        package = super(MS1Record, self).to_dict()
        package['product_scan_ids'] = self.product_scan_ids
        package['msms_peaks'] = self.msms_peaks
        return package


class MSnRecord(MSRecordBase):
    def __init__(self, scan_time=None, neutral_mass=None, mz=None, intensity=None, charge=None,
                 precursor_scan_id=None, product_scan_id=None, defaulted=None, orphan=None,
                 drift_time=None, **kwargs):
        super(MSnRecord, self).__init__(scan_time, drift_time, **kwargs)
        self.neutral_mass = neutral_mass
        self.mz = mz
        self.intensity = intensity
        self.charge = charge
        self.precursor_scan_id = precursor_scan_id
        self.product_scan_id = product_scan_id
        self.defaulted = defaulted
        self.orphan = orphan

    def to_dict(self):
        package = super(MSnRecord, self).to_dict()
        package['neutral_mass'] = self.neutral_mass
        package['mz'] = self.mz
        package['intensity'] = self.intensity
        package['charge'] = self.charge
        package['precursor_scan_id'] = self.precursor_scan_id
        package['product_scan_id'] = self.product_scan_id
        package['defaulted'] = self.defaulted
        package['orphan'] = self.orphan
        return package


class ExtendedScanIndex(object):
    SCHEMA_VERSION = "1.1"

    def __init__(self, ms1_ids=None, msn_ids=None, schema_version=None):
        if schema_version is None:
            schema_version = self.SCHEMA_VERSION
        if ms1_ids is None:
            ms1_ids = {}
        if msn_ids is None:
            msn_ids = {}
        self.ms1_ids = OrderedDict(ms1_ids)
        self.msn_ids = OrderedDict(msn_ids)
        self.schema_version = schema_version

    def get_scan_dict(self, key):
        try:
            return self.ms1_ids[key]
        except KeyError:
            return self.msn_ids[key]

    def __getitem__(self, key):
        return self.get_scan_dict(key)

    def _package_precursor_information(self, product):
        precursor_information = product.precursor_information
        if precursor_information.extracted_neutral_mass != 0:
            charge = precursor_information.extracted_charge
            if charge == ChargeNotProvided:
                charge = 'ChargeNotProvided'
            package = {
                "neutral_mass": precursor_information.extracted_neutral_mass,
                "mz": precursor_information.extracted_mz,
                "intensity": precursor_information.extracted_intensity,
                "charge": charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "product_scan_id": product.id,
                "scan_time": product.scan_time,
                "defaulted": precursor_information.defaulted,
                "orphan": precursor_information.orphan
            }
            if product.has_ion_mobility():
                package['drift_time'] = product.drift_time
        else:
            charge = precursor_information.charge
            if charge == ChargeNotProvided:
                charge = 'ChargeNotProvided'
            package = {
                "neutral_mass": precursor_information.neutral_mass,
                "mz": precursor_information.mz,
                "intensity": precursor_information.intensity,
                "charge": charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "product_scan_id": product.id,
                "scan_time": product.scan_time,
                "defaulted": precursor_information.defaulted,
                "orphan": precursor_information.orphan
            }
            if product.has_ion_mobility():
                package['drift_time'] = product.drift_time
        return package

    def add_scan_bunch(self, bunch):
        package = {
            "scan_time": bunch.precursor.scan_time,
            "product_scan_ids": [
                product.id for product in bunch.products
            ],
            "msms_peaks": [
                p.index.neutral_mass for p in bunch.precursor.deconvoluted_peak_set
                if p.chosen_for_msms
            ] if bunch.precursor.deconvoluted_peak_set is not None else [],
        }
        if bunch.precursor.has_ion_mobility():
            package['drift_time'] = bunch.precursor.drift_time
        self.ms1_ids[bunch.precursor.id] = MS1Record(**package)
        for product in bunch.products:
            self.msn_ids[product.id] = MSnRecord(**self._package_precursor_information(product))

    def serialize(self, handle):
        mapping = {
            "ms1_ids": [(k, v.to_dict()) for k, v in self.ms1_ids.items()],
            "msn_ids": [(k, v.to_dict()) for k, v in self.msn_ids.items()],
            "schema_version": self.schema_version,
        }
        json.dump(mapping, handle)

    def merge(self, other):
        dup = ExtendedScanIndex(self.ms1_ids, self.msn_ids)
        dup.ms1_ids.update(other.ms1_ids)
        dup.msn_ids.update(other.msn_ids)
        return dup

    @staticmethod
    def index_file_name(name):
        return name + '-idx.json'

    @classmethod
    def deserialize(cls, handle):
        mapping = json.load(handle)
        ms1_ids = mapping.get("ms1_ids", [])
        mapping['ms1_ids'] = [(k, MS1Record(**v)) for k, v in ms1_ids]
        msn_ids = mapping.get("msn_ids", [])
        mapping['msn_ids'] = [(k, MSnRecord(**v)) for k, v in msn_ids]
        return cls(**mapping)

    def get_precursor_information(self, bind=None):
        out = []
        for key, info in self.msn_ids.items():
            mz = info['mz']
            neutral_mass = info['neutral_mass']
            charge = info['charge']
            if charge == "ChargeNotProvided":
                charge = ChargeNotProvided
            intensity = info['intensity']
            precursor_scan_id = info['precursor_scan_id']
            product_scan_id = info['product_scan_id']
            pinfo = PrecursorInformation(
                mz, intensity, charge, precursor_scan_id,
                bind, neutral_mass, charge, intensity,
                product_scan_id=product_scan_id)
            out.append(pinfo)
        return out

    def find_msms_by_precursor_mass(self, neutral_mass, mass_error_tolerance=1e-5, bind=None):
        m = neutral_mass
        w = neutral_mass * mass_error_tolerance
        lo = m - w
        hi = m + w
        out = []
        for pinfo in self.get_precursor_information(bind):
            if lo <= pinfo.neutral_mass <= hi:
                out.append(pinfo)
        return out
