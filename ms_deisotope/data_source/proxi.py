import os
import re
import time
import logging
import hashlib
import json

from collections import defaultdict
from threading import RLock, Thread

import idzip

import ms_deisotope
from ms_deisotope.utils import LRUDict
from ms_deisotope.data_source.usi import USI


BLOCK_SIZE = 2 ** 20
logger = logging.getLogger("proxi_backend")


class PROXIDatasetAPIMixinBase(object):
    def _datafile_uri(self, dataset_id, data_file):
        raise NotImplementedError()

    def _is_peaklist(self, d):
        base = os.path.basename(d)
        if base.endswith(".gz"):
            base = base[:-3]
        base = base.lower()
        _, ext = base.rsplit(".", 1)
        if ext in ["mzml", 'mgf', 'mzxml', 'mzmlb', 'ms2']:
            return True
        return False

    def describe_dataset(self, dataset_id):
        record = {
            "accession": dataset_id,
            "title": dataset_id,
            "summary": dataset_id,
            "instruments": [],
            "identifiers": [],
            "species": [],
            "modifications": [],
            "contacts": [],
            "publications": [],
            "keywords": [],
            "datasetLink": None,
            "dataFiles": [],
            "links": [],
        }
        data_files = self._find_mz_files(self.index[dataset_id])
        for d in data_files:
            if self._is_peaklist(d):
                rep = {"accession": "MS:1002850", "name": "Peak list file URI",
                       "value": self._datafile_uri(dataset_id, d)}
            else:
                rep = {"accession": "MS:1002846", "name": "Associated raw file URI",
                       "value": self._datafile_uri(dataset_id, d)}
            record['dataFiles'].append(rep)
        return record

    def enumerate_datasets(self, result_type=None):
        if result_type == 'compact':
            result = list(self.index.keys())
            return result
        return [
            self.describe_dataset(d) for d in self.index.keys()]


class MassSpectraDataArchiveBase(PROXIDatasetAPIMixinBase):
    valid_file_types = ['mzML', 'mzXML', 'mgf', 'mzMLb']

    def __init__(self, cache_size=None, **kwargs):
        if cache_size is None:
            cache_size = 1
        self.cache = LRUDict(cache_size=cache_size)
        self.index = self.build_index()
        self.monitor_thread = self.monitor_index(
            kwargs.get("update_frequency"))

    def _monitor_index(self, update_frequency: float = None):
        if update_frequency is None:
            # Default to one hour, 60^2 seconds
            update_frequency = 60.0 ** 2
        logger.info(
            "Monitoring archive for changes every %0.2f seconds.", update_frequency)
        while True:
            time.sleep(update_frequency)
            logger.info("Rebuilding Index!")
            self.index = self.build_index(use_cache=False)

    def monitor_index(self, update_frequency: float = None):
        update_thread = Thread(
            target=self._monitor_index,
            args=(update_frequency, ), daemon=True)
        update_thread.start()
        return update_thread

    def _checksum(self, string):
        return hashlib.new("md5", string).hexdigest()

    def build_index(self, use_cache=True):
        raise NotImplementedError()

    def _find_index_cache(self, path, use_cache=True):
        cache_name = hashlib.md5(
            path.encode('utf8')).hexdigest() + '.json'
        if os.path.exists(cache_name) and use_cache:
            logger.info("Reading cache from %r", cache_name)
            try:
                with open(cache_name, 'rt') as fh:
                    index = defaultdict(list, json.load(fh))
                    return cache_name, index
            except IOError:
                logger.error("Failed to read index from disk", exc_info=True)
        return cache_name, None

    def _build_index(self, path, use_cache=True):
        file_types = self.valid_file_types + ['json']

        cache_name, index = self._find_index_cache(path, use_cache=use_cache)
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

    def _traverse_files(self, file_types, index):
        raise NotImplementedError()

    def find_ms_files_for(self, dataset, run_name):
        result = []
        members = self.index.get(dataset)

        if members is None:
            raise ValueError("DatasetNotAvailable")

        for member in members:
            base = member.split("/")[-1]
            if base.startswith(run_name):
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

    def open_ms_file(self, mzml_uri, index_uri=None):
        logger.info("Opening %r", mzml_uri)
        mzml_fh = self._opener(mzml_uri)
        if mzml_uri.endswith("gz"):
            mzml_fh = idzip.IdzipFile(fileobj=mzml_fh)
        index_fh = None
        if index_uri is not None:
            index_fh = self._opener(index_uri, block_size=BLOCK_SIZE)
        else:
            logger.warning("Byte offset index is missing for %r", mzml_uri)
        reader = ms_deisotope.MSFileLoader(mzml_fh, index_file=index_fh)
        lock = RLock()
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
            raise ValueError("UnrecognizedIndexFlag", usi.scan_identifier_type)
        scan = reader.get_scan_by_index(index)
        return scan

    def process_scan(self, scan, reader):
        # No processing
        return scan

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

        payload['mzs'] = mzs
        payload['intensities'] = intensities
        if charges is not None:
            payload['charges'] = charges
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

    def build_response(self, usi, scan):
        payload = {
            "attributes": [
                {"accession": "MS:1000511", "name": "ms level",
                    "value": str(scan.ms_level)},
                {"accession": "MS:1008025", "name": "scan number",
                    "value": str(scan.index + 1)},
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
                    "value": str(mz),
                })
            if charge:
                payload['attributes'].append({
                    "accession": "MS:1000041", "name": "charge state",
                    "value": str(charge),
                })
        return payload

    def handle_usi(self, usi, convert_json=True):
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


class ScanProcessingMixin(object):
    def __init__(self, *args, **kwargs):
        self.peak_picking_args = kwargs.get('peak_picking_args', {})
        self.deconvolution_args = kwargs.get('deconvolution_args', {})
        self.deconvolution_args.setdefault("averagine", ms_deisotope.peptide)
        self.deconvolution_args.setdefault(
            "scorer", ms_deisotope.MSDeconVFitter(0))
        self.deconvolution_args.setdefault("truncate_after", 0.8)
        super(ScanProcessingMixin, self).__init__(*args, **kwargs)

    def process_scan(self, scan, reader):
        scan.pick_peaks(**self.peak_picking_args)
        scan.deconvolute(**self.deconvolution_args)
        return scan
