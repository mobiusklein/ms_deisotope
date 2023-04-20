import os
import re
import time
import logging
import hashlib
import json

import dataclasses
from dataclasses import dataclass, field as dfield
from collections import defaultdict
from threading import RLock, Thread, local
from typing import Callable, Dict, List, Optional

import idzip

import ms_deisotope
from ms_deisotope.data_source.scan.base import ScanBase
from ms_deisotope.data_source.scan.loader import RandomAccessScanSource

from ms_deisotope.utils import LRUDict
from ms_deisotope.data_source._compression import get_opener
from ms_deisotope.data_source.usi import USI


BLOCK_SIZE = 2 ** 20
logger = logging.getLogger(__name__)

@dataclass(frozen=True)
class FileType:
    accession: str
    name: str

    def __str__(self):
        return self.name

    def to_dict(self) -> Dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, state):
        return cls(**state)


PeakListFile = FileType("MS:1002850", "Peak list file URI")
AssociatedRawFile = FileType("MS:1002846", "Associated raw file URI")


@dataclass(frozen=True)
class PROXIDataFile:
    uri: str
    file_type: FileType

    def to_dict(self) -> Dict:
        d = self.file_type.to_dict()
        d['value'] = self.uri
        return d

    @classmethod
    def from_dict(cls, d):
        ft = FileType(d['accession'], d['name'])
        return cls(d['value'], ft)


@dataclass
class PROXIDataset:
    accession: str = dfield(default=None)
    title: str = dfield(default=None)
    summary: str = dfield(default=None)
    instruments: List[str] = dfield(default_factory=list)
    identifiers: List[str] = dfield(default_factory=list)
    species: List[str] = dfield(default_factory=list)
    modifications: List[str] = dfield(default_factory=list)
    contacts: List[str] = dfield(default_factory=list)
    publications: List[str] = dfield(default_factory=list)
    keywords: List[str] = dfield(default_factory=list)
    dataset_link: str = dfield(default=None)
    dataset_files: List[PROXIDataFile] = dfield(default_factory=list)
    links: List[str] = dfield(default_factory=list)

    def to_dict(self):
        d = dataclasses.asdict(self)
        d['datasetLinks'] = d.pop("dataset_links", None)
        # TODO: Make this consistent whether or not it is a dict or an instance
        d['dataFiles'] = [PROXIDataFile(f['uri'], FileType.from_dict(f['file_type'])).to_dict() for f in d.pop('dataset_files', [])]
        return d

    @classmethod
    def from_dict(cls, d):
        d = d.copy()
        d['dataset_links'] = d.pop('datasetLinks', None)
        d['dataset_files'] = [
            PROXIDataFile.from_dict(f) for f in d.pop('datasetFiles', [])]
        return cls(**d)

    @classmethod
    def new(cls, dataset_id):
        self = cls()
        self.accession = dataset_id
        self.title = dataset_id
        self.summary = dataset_id
        return self


class PROXIDatasetAPIMixinBase(object):

    def _datafile_uri(self, dataset_id, data_file):
        raise NotImplementedError()

    def _is_peaklist(self, d):
        base = os.path.basename(d)
        if base.endswith(".gz"):
            base = base[:-3]
        base = base.lower()
        _, ext = base.rsplit(".", 1)
        if ext in tuple(map(str.lower, self.peaklist_file_types)):
            return True
        return False

    @property
    def peaklist_file_types(self):
        """
        Alias to allow :attr:`valid_file_types` to vary.

        Returns
        -------
        list[str]
        """
        return self.valid_file_types

    def describe_dataset(self, dataset_id: str) -> PROXIDataset:
        record = PROXIDataset.new(dataset_id)
        data_files = self._find_mz_files(self.index[dataset_id])
        for d in data_files:
            uri = self._datafile_uri(dataset_id, d)
            if self._is_peaklist(d):
                f = PROXIDataFile(uri, PeakListFile)
            else:
                f = PROXIDataFile(uri, AssociatedRawFile)
            record.dataset_files.append(f)
        return record

    def enumerate_datasets(self, result_type=None):
        if result_type == 'compact':
            result = list(self.index.keys())
            return result
        return [
            self.describe_dataset(d).to_dict() for d in self.index.keys()]


class IndexManagingMixin(object):

    def __init__(self, **kwargs):
        self.index = self.build_index()
        self.monitor_thread = self.monitor_index(
            kwargs.get("update_frequency"))

    def _monitor_index(self, update_frequency=None):
        if update_frequency is None:
            # Default to one hour, 60^2 seconds
            update_frequency = 60.0 ** 2
        logger.info(
            "Monitoring archive for changes every %0.2f seconds.", update_frequency)
        while True:
            time.sleep(update_frequency)
            logger.info("Rebuilding Index!")
            self.index = self.build_index(use_cache=False)

    def monitor_index(self, update_frequency=None):
        update_thread = Thread(
            target=self._monitor_index,
            args=(update_frequency, ))
        update_thread.daemon = True
        update_thread.start()
        return update_thread

    def _checksum(self, string):
        return hashlib.new("md5", string).hexdigest()

    def _find_index_cache(self, use_cache=True):
        cache_name = hashlib.md5(
            self.base_uri.encode('utf8')).hexdigest() + '.json'
        if os.path.exists(cache_name) and use_cache:
            logger.info("Reading cache from %r", cache_name)
            try:
                with open(cache_name, 'rt') as fh:
                    index = defaultdict(list, json.load(fh))
                    return cache_name, index
            except IOError:
                logger.error("Failed to read index from disk", exc_info=True)
        return cache_name, None


class SpectrumSerializerServiceMixin(object):
    scan_filters: List[Callable[[
        ScanBase, Optional[RandomAccessScanSource]], ScanBase]]

    def __init__(self, scan_filters: Optional[List]=None):
        self.scan_filters = scan_filters or []

    def process_scan(self, scan: ScanBase, reader: Optional[RandomAccessScanSource]):
        """
        Process ``scan`` with each callable in :attr:`scan_filters` prior to serialization.

        Arguments
        ---------
        scan : :class:`~.ScanBase`
            The scan being processed
        reader : :class:`~.RandomAccessScanSource`
            The data source for ``scan``

        Returns
        -------
        :class:`~.ScanBase`
        """
        for filt in self.scan_filters:
            scan = filt(scan, reader)
        return scan

    def add_scan_filter(self, filt: Callable[[ScanBase, Optional[RandomAccessScanSource]], ScanBase]):
        """
        Add a new scan filtering callable to the chain of filters.

        This new filter is applied after all previous filters.

        Parameters
        ----------
        filt : Callable[[ScanBase, Optional[RandomAccessScanSource]], ScanBase]
            A callable object that transforms a :class:`~.ScanBase` object with optional access to its
            source.
        """
        self.scan_filters.append(filt)

    def build_response(self, usi: USI, scan: ScanBase):
        """
        Build a response :class:`dict` representing ``scan`` for the requested ``usi``.

        Parameters
        ----------
        usi : :class:`~pyteomics.usi.USI`
            The parsed USI requested.
        scan : :class:`~.ScanBase`
            The scan retrieved for this request.

        Returns
        -------
        dict
        """
        payload = {
            "usi": str(usi),
            "status": "READABLE",
            "attributes": [
                {"accession": "MS:1000511", "name": "ms level",
                    "value": scan.ms_level},
                {"accession": "MS:1008025", "name": "scan number",
                    "value": scan.index + 1},
                {"accession": "MS:1003061", "name": "spectrum name", "value": scan.id},
            ]
        }
        payload = self._write_peak_arrays(usi, scan, payload)
        if scan.ms_level > 1:
            mz = None
            charge = None
            if scan.precursor_information:
                mz = scan.precursor_information.mz
                charge = scan.precursor_information.charge
            elif scan.isolation_window:
                mz = scan.isolation_window.target
            if mz:
                payload['attributes'].append({
                    "accession": "MS:1000827", "name": "isolation window target m/z",
                    "value": mz,
                })
            if charge:
                payload['attributes'].append({
                    "accession": "MS:1000041", "name": "charge state",
                    "value": charge,
                })
        return payload

    def _write_peak_arrays(self, usi, scan, payload):
        mzs = []
        intensities = []
        charges = None
        base_peak = None
        nb_peaks = None
        if scan.deconvoluted_peak_set is not None:
            mzs = [p.mz for p in scan.deconvoluted_peak_set]
            intensities = [p.intensity for p in scan.deconvoluted_peak_set]
            charges = [p.charge for p in scan.deconvoluted_peak_set]
            is_profile = False
            nb_peaks = len(scan.deconvoluted_peak_set)
            base_peak = scan.base_peak.deconvoluted().intensity
        elif scan.peak_set is not None:
            mzs = [p.mz for p in scan.peak_set]
            intensities = [p.intensity for p in scan.peak_set]
            charges = None
            is_profile = False
            nb_peaks = len(scan.peak_set)
            base_peak = scan.base_peak.centroided().intensity
        else:
            mzs = scan.arrays.mz.tolist()
            intensities = scan.arrays.intensity.tolist()
            is_profile = scan.is_profile
            nb_peaks = len(mzs)
            base_peak = scan.base_peak.raw().intensity

        payload['m/z array'] = mzs
        payload['intensity array'] = intensities
        if charges is not None:
            payload['charge array'] = charges
        payload['attributes'].append(
            {"accession": "MS:1000128", "name": "profile spectrum"} if is_profile else
            {"accession": "MS:1000127", "name": "centroid spectrum"})
        payload['attributes'].append({
            "accession": 'MS:1000505', 'name': 'base peak intensity', 'value': base_peak
        })
        payload['attributes'].append({
            "accession": 'MS:1003059', 'name': 'number of peaks', 'value': nb_peaks
        })
        return payload


class DatasetNotAvailable(ValueError):
    """An error code indicating that the requested dataset isn't available."""


class UnrecognizedIndexFlag(ValueError):
    """An error code indicating that the USI index flag was not recognized."""


class MassSpectraDataArchiveBase(SpectrumSerializerServiceMixin, PROXIDatasetAPIMixinBase, IndexManagingMixin):
    valid_file_types = ['mzML', 'mzXML', 'mgf', 'mzMLb']
    special_openers = set()

    def __init__(self, base_uri, cache_size=None, scan_filters: Optional[List] = None, **kwargs):
        if cache_size is None:
            cache_size = 1
        self.base_uri = base_uri
        self.cache = LRUDict(cache_size=cache_size)
        IndexManagingMixin.__init__(self, **kwargs)
        SpectrumSerializerServiceMixin.__init__(self, scan_filters=scan_filters)

    def build_index(self, use_cache=True):
        file_types = self.valid_file_types + ['json']

        cache_name, index = self._find_index_cache(use_cache=use_cache)
        if index is not None and use_cache:
            return index

        index = defaultdict(list)

        self._traverse_files(file_types, index)

        reduced = {
            k: self._find_mz_files(v) for k, v in index.items()
        }

        n_datasets = len(reduced)
        n_ms_files = 0
        for _, v in reduced.items():
            n_ms_files += len(v)

        logger.info(
            "Collected %d Datasets with %d MS Data Files in Index", n_datasets, n_ms_files)
        logger.info("Saving index to %r", cache_name)
        try:
            with open(cache_name, 'wt') as fh:
                json.dump(index, fh, indent=2, sort_keys=True)
        except IOError:
            logger.error("Failed to write index to disk", exc_info=True)
        return index

    def _traverse_files(self, file_types, index):
        raise NotImplementedError()

    def find_ms_files_for(self, dataset, run_name):
        result = []
        members = self.index.get(dataset)
        if run_name is None:
            return members

        if members is None:
            raise DatasetNotAvailable((dataset, run_name))

        for member in members:
            base = member.split("/")[-1]
            if base.startswith(run_name):
                result.append(member)
            elif base.endswith("json"):
                if base.startswith(run_name.rsplit(".", 1)[0]):
                    result.append(member)
        # Validate that there is only one MS file for the name
        i = 0
        for res in sorted(result):
            if re.search(r"-byte-offsets.json$", res):
                continue
            if re.search(r"-idx.json$", res):
                continue
            i += 1
        if i > 1:
            raise ValueError("Multiple MS Data Files Found")
        return self._find_mz_files(result, first=True), self._find_index_files(result, first=True)

    def _find_mz_files(self, namelist, first=False):
        pattern = re.compile(r"(?:%s)(?:.gz)?$" %
                             "|".join(self.valid_file_types))
        names = []
        for name in namelist:
            if pattern.search(name):
                if first:
                    return name
                names.append(name)
        if first and not names:
            return None
        return names

    def _find_index_files(self, namelist, first=False):
        names = []
        for name in namelist:
            if re.search(r"-byte-offsets.json$", name):
                if first:
                    return name
                names.append(name)
        if first and not names:
            return None
        return names

    def _opener(self, uri, block_size=None, **kwargs):
        raise NotImplementedError()

    def open_ms_file(self, ms_file_uri, index_uri=None):
        if ms_file_uri in self.cache:
            logger.info("Found %r in the cache", ms_file_uri)
            return self.cache[ms_file_uri]
        logger.info("Opening %r", ms_file_uri)
        base_name =  ms_file_uri.split("/")[-1]
        if base_name.split('.')[-1] in self.special_openers:
            reader = ms_deisotope.MSFileLoader(ms_file_uri)
        else:
            mzml_fh = self._opener(ms_file_uri)
            if ms_file_uri.endswith("gz"):
                mzml_fh = idzip.IdzipFile(fileobj=mzml_fh)
            index_fh = None
            if index_uri is not None:
                index_fh = self._opener(index_uri, block_size=BLOCK_SIZE)
            else:
                logger.warning("Byte offset index is missing for %r", ms_file_uri)
            reader = ms_deisotope.MSFileLoader(mzml_fh, index_file=index_fh)
        logger.info("Finished opening %r", ms_file_uri)
        lock = RLock()
        self.cache[ms_file_uri] = (reader, lock)
        return reader, lock

    def list_collections(self):
        return list(self.index.keys())

    def list_datasets(self):
        return {
            k: [f.split("/")[-1] for f in self._find_mz_files(v)] for k, v in self.index.items()
        }

    def get_scan_from_usi(self, usi, reader):
        if usi.scan_identifier_type == 'scan':
            index = int(usi.scan_identifier) - 1
        elif usi.scan_identifier_type == "index":
            index = int(usi.scan_identifier)
        else:
            raise UnrecognizedIndexFlag(usi.scan_identifier_type)
        scan = reader.get_scan_by_index(index)
        return scan

    def enumerate_spectra_for(self, dataset, ms_run=None):
        raise NotImplementedError()

    def handle_usi(self, usi, convert_json=True):
        """
        Receive a USI, parse it into its constituent pieces, locate
        the relevant spectrum, and return it formatted message.

        Parameters
        ----------
        usi : str
            The identifier of the spectrum to be retrieved
        convert_json : bool, optional
            Whether to format the spectrum as JSON or not. Defaults to :const:`True`.

        Returns
        -------
        spectrum : dict or Scan
            The spectrum, rendered according to the `convert_json` or not.
        """
        usi = USI.parse(str(usi))
        mzml_uri, index_uri = self.find_ms_files_for(usi.dataset, usi.datafile)
        reader, lock = self.open_ms_file(mzml_uri, index_uri)
        with lock:
            scan = self.get_scan_from_usi(usi, reader)
            scan = self.process_scan(scan, reader)
            if convert_json:
                payload = self.build_response(usi, scan)
                return payload
            else:
                return scan

    def __call__(self, usi):
        return self.handle_usi(usi)


class FSMassSpectraDataArchive(MassSpectraDataArchiveBase):
    def __init__(self, base_uri, cache_size=None, **kwargs):
        self.valid_file_types = self.valid_file_types[:]
        from ms_deisotope.data_source.thermo_raw_net import determine_if_available
        if determine_if_available():
            self.valid_file_types.append("raw")
            self.special_openers = set(self.special_openers) | {"raw"}
        super(FSMassSpectraDataArchive, self).__init__(
            base_uri, cache_size=cache_size, **kwargs)

    def _opener(self, uri, block_size=None, **kwargs):
        return open(uri, 'rb')

    def _datafile_uri(self, dataset_id, data_file):
        return os.path.normpath('/'.join([self.base_uri, data_file])).replace(os.sep, '/')

    def _traverse_files(self, file_types, index):
        for prefix, dirs, files in os.walk(self.base_uri):
            if os.sep in prefix:
                tokens = prefix.split(os.sep)
            else:
                tokens = prefix.split("/")
            if tokens[-1] == '':
                base = tokens[-2]
            else:
                base = tokens[-1]

            for f in files:
                try:
                    ext = f.rsplit(".", 1)[1]
                except IndexError:
                    continue
                if ext == 'gz':
                    ext = f.rsplit(".", 2)[1]
                if ext in file_types:
                    index[base].append(
                        os.path.normpath(
                            os.path.join(prefix, f)).replace(
                                os.sep, "/"))


class S3MassSpectraDataArchive(MassSpectraDataArchiveBase):
    def __init__(self, base_uri, cache_size=None, s3_args=None, **kwargs):
        if s3_args is None:
            s3_args = {}
        self._tls = local()
        self.s3_args = s3_args
        super(S3MassSpectraDataArchive, self).__init__(
            base_uri=base_uri, cache_size=cache_size, **kwargs)

    @property
    def s3_service(self):
        try:
            return self._tls.s3_service
        except (AttributeError, KeyError):
            from s3fs import S3FileSystem
            self._tls.s3_service = S3FileSystem(**self.s3_args)
            return self._tls.s3_service

    def _opener(self, uri, block_size=None, **kwargs):
        if block_size is None:
            block_size = BLOCK_SIZE
        return self.s3_service.open(uri, block_size=block_size)

    def _datafile_uri(self, dataset_id, data_file):
        return "s3://" + data_file

    def _traverse_files(self, file_types, index):
        for prefix, _dirs, files in self.s3_service.walk(self.base_uri):
            for f in files:
                query = f
                if query.endswith(".gz"):
                    query = query[:-3]
                ext = query.rsplit(".", 1)[-1]
                for t in file_types:
                    if t == ext:
                        index[prefix.replace(self.base_uri[5:], '')].append(
                            prefix + '/' + f)
                        break
        for group, members in list(index.items()):
            if "/" in group:
                base, _rest = group.split("/", 1)
                index[base].extend(members)
                index.pop(group)


class PeakPickingScanFilter(object):
    def __init__(self, peak_picking_args=None, *args, **kwargs):
        if peak_picking_args is None:
            peak_picking_args = {}
        self.peak_picking_args = peak_picking_args

    def __call__(self, scan, reader):
        scan.pick_peaks(**self.peak_picking_args)
        return scan


class DeconvolutingScanFilter(PeakPickingScanFilter):
    def __init__(self, *args, averagine=ms_deisotope.peptide, scorer=ms_deisotope.MSDeconVFitter(0), truncate_after: float=0.8, **kwargs):
        self.deconvolution_args = kwargs.get('deconvolution_args', {})
        self.deconvolution_args.setdefault("averagine", averagine)
        self.deconvolution_args.setdefault("scorer", scorer)
        self.deconvolution_args.setdefault("truncate_after", truncate_after)
        super(DeconvolutingScanFilter, self).__init__(*args, **kwargs)

    def __call__(self, scan, reader):
        scan = super().__call__(scan, reader)
        scan.deconvolute(**self.deconvolution_args)
        return scan
