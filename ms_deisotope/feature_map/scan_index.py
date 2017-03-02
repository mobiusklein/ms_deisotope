import json
from collections import OrderedDict

from ms_deisotope.averagine import neutral_mass


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
                "scan_time": product.scan_time
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
