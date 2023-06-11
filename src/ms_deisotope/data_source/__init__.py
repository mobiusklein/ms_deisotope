try:
    from ms_deisotope._c.units import patch_pyteomics
    patch_pyteomics()
except ImportError:
    pass

from .scan import (
    ScanProxyContext,
    ScanBase,
    Scan,
    ProcessedScan,
    ScanBunch,
    PrecursorInformation,
    ScanDataSource,
    ScanIterator,
    RandomAccessScanSource,
    IonMobilityFrame,
    ProcessedIonMobilityFrame,
    Generic3DIonMobilityFrameSource,
    FramedIonMobilityFrameSource,
    IonMobilitySource,
    IonMobilitySourceRandomAccessFrameSource
)

from typing import Any
from .infer_type import MSFileLoader
from .mzml import MzMLLoader
from .mzxml import MzXMLLoader
from .mgf import MGFLoader
from .common import (
    IsolationWindow,
    ChargeNotProvided
)

from .metadata.file_information import (
    FileInformation,
    FileFormat, file_formats,
    FileContent, content_keys,
    IDFormat, id_formats)

from .metadata.instrument_components import (
    Component, ComponentGroup, InstrumentInformation,
    components)

from .metadata.scan_traits import (
    ScanAcquisitionInformation, ScanEventInformation,
    ScanWindow, scan_attributes)

from .metadata.activation import (
    ActivationInformation, MultipleActivationInformation,
    DissociationMethod, dissociation_methods)


from ._compression import get_opener

from .text import scan_from_csv
from .memory import make_scan, ScanCollection, make_scan_collection


ProcessedRandomAccessScanSource = RandomAccessScanSource[Any, ProcessedScan]


__all__ = [
    "MSFileLoader", "MzMLLoader",
    "MzXMLLoader", "MGFLoader",

    "Scan", "ProcessedScan", "ScanBase",

    "PrecursorInformation",

    "ActivationInformation", "MultipleActivationInformation",
    "DissociationMethod", "dissociation_methods",

    "ScanAcquisitionInformation", "ScanEventInformation",
    "IsolationWindow", "ScanWindow", "scan_attributes",

    "FileInformation", "FileFormat", "file_formats",
    "FileContent", "content_keys",
    "IDFormat", "id_formats",

    "Component", "ComponentGroup", "InstrumentInformation",
    "components",

    "ScanDataSource", "ScanIterator", "ScanBunch",
    "ScanWindow", "RandomAccessScanSource", "ChargeNotProvided",
    "get_opener", "ScanProxyContext", 'scan_from_csv',
    "make_scan", "ScanCollection", "make_scan_collection",
    "ProcessedRandomAccessScanSource",

    "IonMobilityFrame", "ProcessedIonMobilityFrame",
    "Generic3DIonMobilityFrameSource", "FramedIonMobilityFrameSource",
    "IonMobilitySource", "IonMobilitySourceRandomAccessFrameSource",
]
