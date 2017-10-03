import json
from collections import OrderedDict

from ms_deisotope.averagine import neutral_mass
from ms_deisotope.data_source.common import PrecursorInformation


class ExtendedScanIndex(object):
    SCHEMA_VERSION = "1.0"

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
            package = {
                "neutral_mass": precursor_information.extracted_neutral_mass,
                "mz": precursor_information.extracted_mz,
                "intensity": precursor_information.extracted_intensity,
                "charge": precursor_information.extracted_charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "product_scan_id": product.id,
                "scan_time": product.scan_time,
                "defaulted": precursor_information.defaulted,
                "orphan": precursor_information.orphan
            }
        else:
            package = {
                "neutral_mass": neutral_mass(
                    precursor_information.mz, precursor_information.charge),
                "mz": precursor_information.mz,
                "intensity": precursor_information.intensity,
                "charge": precursor_information.charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "product_scan_id": product.id,
                "scan_time": product.scan_time
            }
        return package

    def add_scan_bunch(self, bunch):
        self.ms1_ids[bunch.precursor.id] = {
            "scan_time": bunch.precursor.scan_time,
            "product_scan_ids": [
                product.id for product in bunch.products
            ],
            "msms_peaks": [
                p.index.neutral_mass for p in bunch.precursor.deconvoluted_peak_set
                if p.chosen_for_msms
            ] if bunch.precursor.deconvoluted_peak_set is not None else [],
        }
        for product in bunch.products:
            self.msn_ids[product.id] = self._package_precursor_information(product)

    def serialize(self, handle):
        mapping = {
            "ms1_ids": list(self.ms1_ids.items()),
            "msn_ids": list(self.msn_ids.items()),
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
        return cls(**mapping)

    def get_precursor_information(self, bind=None):
        out = []
        for key, info in self.msn_ids.items():
            mz = info['mz']
            neutral_mass = info['neutral_mass']
            charge = info['charge']
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
