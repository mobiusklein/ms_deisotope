import io
import json

from typing import Any, ClassVar, Dict, List, Optional, OrderedDict

from ms_deisotope.data_source.common import (
    PrecursorInformation, ChargeNotProvided,
    ActivationInformation, ScanBase)
from ms_deisotope.data_source.scan.base import ScanBunch
from ms_deisotope.data_source.scan.loader import ScanDataSource, ScanIterator

from .feature_map import NeutralMassIndex

from ms_deisotope.utils import Base
from ms_deisotope.qc.isolation import CoIsolation


class MSRecordBase(Base):
    scan_time: float
    drift_time: float
    _extra_keys: List[str]

    def __init__(self, scan_time, drift_time=None, **kwargs):
        self.scan_time = scan_time
        self.drift_time = drift_time
        self._extra_keys = list(kwargs.keys())
        self.__dict__.update(kwargs)

    def __getitem__(self, key):
        return getattr(self, key)

    def __eq__(self, other):
        return self.to_dict() == other.to_dict()

    def __ne__(self, other):
        return not self == other

    def to_dict(self) -> Dict[str, Any]:
        package = {
            "scan_time": self.scan_time,
        }
        if self.drift_time is not None:
            package['drift_time'] = self.drift_time
        for key in self._extra_keys:
            package[key] = self[key]
        return package

    def get(self, key: str, default=None):
        try:
            return self[key]
        except Exception:
            return default


class MS1Record(MSRecordBase):
    product_scan_ids: List[str]
    msms_peaks: List

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
    mz: float
    neutral_mass: float
    intensity: float
    charge: int
    precursor_scan_id: str
    product_scan_id: str
    defaulted: bool
    orphan: bool
    coisolation: List[CoIsolation]
    activation: ActivationInformation

    def __init__(self, scan_time=None, neutral_mass=None, mz=None, intensity=None, charge=None,
                 precursor_scan_id=None, product_scan_id=None, defaulted=None, orphan=None,
                 drift_time=None, coisolation=None, activation=None, **kwargs):
        super(MSnRecord, self).__init__(scan_time, drift_time, **kwargs)
        self.neutral_mass = neutral_mass
        self.mz = mz
        self.intensity = intensity
        self.charge = charge
        self.precursor_scan_id = precursor_scan_id
        self.product_scan_id = product_scan_id
        self.defaulted = defaulted
        self.orphan = orphan
        self.coisolation = [CoIsolation(*c) if c else c for c in (coisolation or [])]
        self.activation = ActivationInformation.from_dict(
            activation) if isinstance(activation, dict) else activation

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
        package['coisolation'] = self.coisolation
        package['activation'] = self.activation.to_dict() if self.activation is not None else None
        return package


class ExtendedScanIndex(object):
    """
    An extra index that holds scan-level metadata in memory independent of
    the :class:`~.ScanBase` object itself.

    Attributes
    ----------
    ms1_ids : :class:`~.OrderedDict`
        An ordered mapping from :term:`scan_id` to :class:`MS1Record`
    msn_ids : :class:`~.OrderedDict`
        An ordered mapping from :term:`scan_id` to :class:`MSNRecord`
    schema_version : str
        The version string for the schema of the serialization format
    _index_bind : :class:`~.RandomAccessScanSource`
        The data source to bind when calling :meth:`get_precursor_information` and
        :meth:`find_msms_by_precursor_mass`
    _mass_search_index : :class:`~.NeutralMassIndex`
        A fast-to-search collection to make :meth:`find_msms_by_precursor_mass` faster
    """

    SCHEMA_VERSION: ClassVar[str] = "1.1"

    ms1_ids: OrderedDict[str, MS1Record]
    msn_ids: OrderedDict[str, MSnRecord]
    schema_version: str

    _index_bind: Optional[ScanDataSource]
    _mass_search_index: Optional[NeutralMassIndex]

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

        self._index_bind = None
        self._mass_search_index = None

    def clear(self):
        """Discard all information held in memory, analogous to :meth:`dict.clear`."""
        self.ms1_ids.clear()
        self.msn_ids.clear()
        self._mass_search_index = None

    def get_scan_dict(self, key):
        try:
            return self.ms1_ids[key]
        except KeyError:
            return self.msn_ids[key]

    def __getitem__(self, key):
        return self.get_scan_dict(key)

    def _package_precursor_information(self, product: ScanBase) -> Dict[str, Any]:
        precursor_information = product.precursor_information
        package = {
            "product_scan_id": product.id,
            "scan_time": product.scan_time,
            "ms_level": product.ms_level,
        }
        if precursor_information is None:
            pass
        elif precursor_information.extracted_neutral_mass != 0:
            charge = precursor_information.extracted_charge
            if charge == ChargeNotProvided:
                charge = 'ChargeNotProvided'
            package.update({
                "neutral_mass": precursor_information.extracted_neutral_mass,
                "mz": precursor_information.extracted_mz,
                "intensity": precursor_information.extracted_intensity,
                "charge": charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "defaulted": precursor_information.defaulted,
                "orphan": precursor_information.orphan,
                "coisolation": precursor_information.coisolation,
                "activation": product.activation,
            })
        else:
            charge = precursor_information.charge
            if charge == ChargeNotProvided:
                charge = 'ChargeNotProvided'
            package.update({
                "neutral_mass": precursor_information.neutral_mass,
                "mz": precursor_information.mz,
                "intensity": precursor_information.intensity,
                "charge": charge,
                "precursor_scan_id": precursor_information.precursor_scan_id,
                "defaulted": precursor_information.defaulted,
                "orphan": precursor_information.orphan,
                "coisolation": precursor_information.coisolation,
                "activation": product.activation,
            })
        try:
            if product.has_ion_mobility():
                package['drift_time'] = product.drift_time
        except AttributeError:
            pass
        return package

    def add_scan(self, scan: ScanBase):
        """
        Add ``scan`` to the index.

        Parameters
        ----------
        scan: :class:`~.ScanBase`
            The scan to add to the index
        """
        if scan.ms_level == 1:
            package = {
                "scan_time": scan.scan_time,
                # would be populated if add_scan_bunch were used
                "product_scan_ids": [],
            }
            try:
                if scan.has_ion_mobility():
                    package['drift_time'] = scan.drift_time
            except AttributeError:
                pass
            self.ms1_ids[scan.id] = MS1Record(**package)
        else:
            self.msn_ids[scan.id] = MSnRecord(**self._package_precursor_information(scan))

    def add_scan_bunch(self, bunch: ScanBunch):
        """
        Add each scan object in ``bunch`` to the index.

        Parameters
        ----------
        bunch : :class:`~.ScanBunch`
            The bunch of scans to add.
        """
        if bunch.precursor is not None:
            package = {
                "scan_time": bunch.precursor.scan_time,
                "product_scan_ids": [
                    product.id for product in bunch.products
                ],
            }
            try:
                if bunch.precursor.has_ion_mobility():
                    package['drift_time'] = bunch.precursor.drift_time
            except AttributeError:
                pass
            self.ms1_ids[bunch.precursor.id] = MS1Record(**package)
        for product in bunch.products:
            self.msn_ids[product.id] = MSnRecord(**self._package_precursor_information(product))

    def update_from_reader(self, reader: ScanIterator):
        """
        Iterate over ``reader``, accumulating scans in the index.

        Parameters
        ----------
        reader: :class:`~.ScanIterator`
            The scan iterator to consume scans from
        """
        for bunch in reader:
            if isinstance(bunch, ScanBase):
                self.add_scan(bunch)
            else:
                self.add_scan_bunch(bunch)

    def dump(self, handle: io.TextIOBase):
        """
        Serialize the index to JSON.

        Parameters
        ----------
        handle: file-like
            The file-like object to write the index to
        """
        mapping = {
            "ms1_ids": [(k, v.to_dict()) for k, v in self.ms1_ids.items()],
            "msn_ids": [(k, v.to_dict()) for k, v in self.msn_ids.items()],
            "schema_version": self.schema_version,
        }
        json.dump(mapping, handle)

    serialize = dump

    def merge(self, other: 'ExtendedScanIndex'):
        """
        Combine the indices in ``other`` with those in ``self``,
        return a copy containing both collections' data.

        Parameters
        ----------
        other: :class:`ExtendedScanIndex`
            Another scan index to combine with this one

        Returns
        -------
        :class:`ExtendedScanIndex`
        """
        dup = ExtendedScanIndex(self.ms1_ids, self.msn_ids)
        dup.ms1_ids.update(other.ms1_ids)
        dup.msn_ids.update(other.msn_ids)
        return dup

    @staticmethod
    def index_file_name(name: str) -> str:
        """
        Create a standard file name based on source file name ``name``
        for storing the index

        Parameters
        ----------
        name: str
            The path to the source file to create an adjacent index file
            name for.

        Returns
        -------
        str
        """
        return name + '-idx.json'

    @classmethod
    def load(cls, handle: io.TextIOBase):
        """
        Construct a :class:`ExtendedScanIndex` instance from a file object

        Parameters
        ----------
        handle : file-like
            The file to read from

        Returns
        -------
        :class:`ExtendedScanIndex`
        """
        mapping = json.load(handle)
        ms1_ids = mapping.get("ms1_ids", [])
        mapping['ms1_ids'] = [(k, MS1Record(**v)) for k, v in ms1_ids]
        msn_ids = mapping.get("msn_ids", [])
        mapping['msn_ids'] = [(k, MSnRecord(**v)) for k, v in msn_ids]
        return cls(**mapping)

    deserialize = load

    def get_precursor_information(self, bind: Optional[ScanDataSource]=None) -> List[PrecursorInformation]:
        """
        Create a list of :class:`~.PrecursorInformation` objects
        from :attr:`msn_ids`'s records.

        Returns
        -------
        list
        """
        out = []
        for _, info in self.msn_ids.items():
            mz = info['mz']
            neutral_mass = info['neutral_mass']
            if neutral_mass is None:
                continue
            charge = info['charge']
            if charge == "ChargeNotProvided":
                charge = ChargeNotProvided
            intensity = info['intensity']
            precursor_scan_id = info['precursor_scan_id']
            product_scan_id = info['product_scan_id']
            pinfo = PrecursorInformation(
                mz, intensity, charge, precursor_scan_id,
                bind, neutral_mass, charge, intensity,
                product_scan_id=product_scan_id,
                coisolation=info.get('coisolation'))
            out.append(pinfo)
        return out

    def _get_mass_search_index(self, bind: Optional[ScanDataSource]=None) -> NeutralMassIndex:
        if self._mass_search_index is not None and (bind is self._index_bind or bind is None):
            return self._mass_search_index
        pinfos = self.get_precursor_information(bind)
        index = NeutralMassIndex(pinfos)
        self._index_bind = bind
        self._mass_search_index = index
        return index

    def find_msms_by_precursor_mass(self, neutral_mass: float, mass_error_tolerance: float=1e-5, bind: Optional[ScanDataSource]=None) -> List[PrecursorInformation]:
        """
        Find all entries in :attr:`msn_ids` which are within
        ``mass_error_tolerance`` of ``neutral_mass``.

        Returns
        -------
        list
        """
        index = self._get_mass_search_index(bind)
        out = index.find_all(neutral_mass, mass_error_tolerance)
        return out
