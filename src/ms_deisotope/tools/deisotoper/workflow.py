from os import PathLike
from typing import IO, Callable, ClassVar, Dict, Any, Optional, Type

import ms_deisotope
from ms_deisotope.data_source.scan.loader import ScanDataSource
from ms_deisotope.task import TaskBase

from .output import ThreadedMzMLScanStorageHandler, NullScanStorageHandler, ScanStorageHandlerBase
from .scan_generator import ScanGenerator


class ScanSink(object):
    scan_generator: ScanGenerator
    scan_store: ScanStorageHandlerBase
    _scan_store_type: Type[ScanStorageHandlerBase]

    def __init__(self, scan_generator: ScanGenerator, storage_type: Type[ScanStorageHandlerBase] = NullScanStorageHandler):
        self.scan_generator = scan_generator
        self.scan_store = None
        self._scan_store_type = storage_type

    @property
    def scan_source(self) -> Optional[PathLike]:
        try:
            return self.scan_generator.scan_source
        except AttributeError:
            return None

    @property
    def sample_run(self):
        try:
            return self.scan_store.sample_run
        except AttributeError:
            return None

    def configure_storage(self, storage_path: PathLike=None,
                          name: Optional[str]=None,
                          source: Optional[ScanDataSource]=None,
                          stream_cls: Optional[Callable[[str], IO]] = None):
        self.scan_store = self._scan_store_type.configure_storage(
            storage_path, name, source, stream_cls=stream_cls)
        packer_type = self.scan_store.get_scan_packer_type()
        if source is not None:
            source._scan_packer_type = packer_type

    def configure_iteration(self, *args, **kwargs):
        self.scan_generator.configure_iteration(*args, **kwargs)

    def store_scan(self, scan):
        if self.scan_store is not None:
            self.scan_store.accumulate(scan)

    def commit(self):
        if self.scan_store is not None:
            self.scan_store.commit()

    def complete(self):
        if self.scan_store is not None:
            self.scan_store.complete()
        self.scan_generator.close()

    def next_scan(self):
        scan = next(self.scan_generator)
        self.store_scan(scan)
        return scan

    def __iter__(self):
        return self

    def __next__(self):
        return self.next_scan()

    def next(self):
        return self.next_scan()


class SampleConsumer(TaskBase):
    MS1_ISOTOPIC_PATTERN_WIDTH: float = 0.95
    MS1_IGNORE_BELOW: float = 0.05
    MSN_ISOTOPIC_PATTERN_WIDTH: float = 0.80
    MSN_IGNORE_BELOW: float = 0.05

    MS1_SCORE_THRESHOLD: float = 20.0
    MSN_SCORE_THRESHOLD: float = 10.0

    storage_type: Type[ScanStorageHandlerBase]
    scan_generator: ScanGenerator

    ms1_processing_args: Dict[str, Dict[str, Any]]
    msn_processing_args: Dict[str, Dict[str, Any]]

    ms1_averaging: int
    extract_only_tandem_envelopes: bool
    ignore_tandem_scans: bool
    deconvolute: bool

    n_processes: int

    start_scan_id: str
    end_scan_id: str
    start_scan_time: float
    end_scan_time: float

    verbose: bool

    stream_cls: Optional[Callable[[str], IO]]

    scan_generator_cls: ClassVar[Type[ScanGenerator]] = ScanGenerator

    def __init__(self, ms_file: str,
                 ms1_peak_picking_args: Optional[dict] = None,
                 msn_peak_picking_args: Optional[dict] = None,
                 ms1_deconvolution_args: Optional[dict] = None,
                 msn_deconvolution_args: Optional[dict] = None,
                 start_scan_id: Optional[str] = None,
                 end_scan_id: Optional[str] = None,
                 storage_path: Optional[str] = None,
                 sample_name: Optional[str]=None,
                 storage_type: Optional[Type[ScanStorageHandlerBase]]=None,
                 n_processes: int=5,
                 extract_only_tandem_envelopes: bool = False,
                 ignore_tandem_scans: bool = False,
                 ms1_averaging: int=0,
                 default_precursor_ion_selection_window: Optional[float] = 1.5,
                 deconvolute: bool=True,
                 verbose: bool=False,
                 start_scan_time: Optional[float] = None,
                 end_scan_time: Optional[float] = None,
                 stream_cls: Optional[Callable[[str], IO]] = None):

        if storage_type is None:
            storage_type = ThreadedMzMLScanStorageHandler

        self.ms_file = ms_file
        self.storage_path = storage_path
        self.sample_name = sample_name

        self.n_processes = n_processes
        self.storage_type = storage_type
        self.extract_only_tandem_envelopes = extract_only_tandem_envelopes
        self.ignore_tandem_scans = ignore_tandem_scans
        self.ms1_averaging = ms1_averaging
        self.default_precursor_ion_selection_window = default_precursor_ion_selection_window

        # for display purposes only
        self.ms1_processing_args = {
            "peak_picking": ms1_peak_picking_args,
        }
        self.msn_processing_args = {
            "peak_picking": msn_peak_picking_args,
        }

        self.deconvolute = deconvolute

        if deconvolute:
            self.ms1_processing_args["deconvolution"] = ms1_deconvolution_args
            self.msn_processing_args["deconvolution"] = msn_deconvolution_args

        n_helpers = max(self.n_processes - 1, 0)
        self.scan_generator = self.scan_generator_cls(
            ms_file,
            number_of_helpers=n_helpers,
            ms1_peak_picking_args=ms1_peak_picking_args,
            msn_peak_picking_args=msn_peak_picking_args,
            ms1_deconvolution_args=ms1_deconvolution_args,
            msn_deconvolution_args=msn_deconvolution_args,
            extract_only_tandem_envelopes=extract_only_tandem_envelopes,
            ignore_tandem_scans=ignore_tandem_scans,
            ms1_averaging=ms1_averaging,
            default_precursor_ion_selection_window=default_precursor_ion_selection_window,
            deconvolute=deconvolute,
            verbose=verbose)

        self.start_scan_id = start_scan_id
        self.end_scan_id = end_scan_id
        self.start_scan_time = start_scan_time
        self.end_scan_time = end_scan_time
        self.stream_cls = stream_cls
        self.sample_run = None

    @classmethod
    def default_processing_configuration(cls, averagine=ms_deisotope.peptide, msn_averagine=None):
        if msn_averagine is None:
            msn_averagine = averagine

        ms1_peak_picking_args = {
            "transforms": [
            ]
        }

        ms1_deconvolution_args = {
            "scorer": ms_deisotope.scoring.PenalizedMSDeconVFitter(20, 2.),
            "max_missed_peaks": 3,
            "averagine": averagine,
            "truncate_after": cls.MS1_ISOTOPIC_PATTERN_WIDTH,
            "ignore_below": cls.MS1_IGNORE_BELOW,
            "deconvoluter_type": ms_deisotope.AveraginePeakDependenceGraphDeconvoluter,
            "use_quick_charge": True,
        }

        msn_peak_picking_args = {}

        msn_deconvolution_args = {
            "scorer": ms_deisotope.scoring.MSDeconVFitter(10),
            "averagine": msn_averagine,
            "max_missed_peaks": 1,
            "truncate_after": cls.MSN_ISOTOPIC_PATTERN_WIDTH,
            "ignore_below": cls.MSN_IGNORE_BELOW,
            "deconvoluter_type": ms_deisotope.AveraginePeakDependenceGraphDeconvoluter,
            "use_quick_charge": True,
        }

        return (ms1_peak_picking_args, msn_peak_picking_args,
                ms1_deconvolution_args, msn_deconvolution_args)

    def run(self):
        self.log("Initializing Generator")
        sink = ScanSink(self.scan_generator, self.storage_type)
        self.log("Initializing Storage")
        sink.configure_storage(self.storage_path, self.sample_name,
                               self.scan_generator, stream_cls=self.stream_cls)
        self.scan_generator.configure_iteration(self.start_scan_id, self.end_scan_id)
        self.log("Setting Sink")

        self.log("Begin Processing")
        last_scan_time = 0
        i = 0
        for scan in sink:
            i += 1
            if (scan.scan_time - last_scan_time > 1.0) or (i % 1000 == 0):
                percent_complete = None
                if self.end_scan_time is not None:
                    percent_complete = (
                        scan.scan_time - self.start_scan_time) / (self.end_scan_time - self.start_scan_time)
                if percent_complete is not None:
                    self.log("Processed %s (time: %0.3f %0.2f%% Done)" % (
                        scan.id, scan.scan_time, percent_complete * 100))
                else:
                    self.log("Processed %s (time: %0.3f)" % (
                        scan.id, scan.scan_time,))
                last_scan_time = scan.scan_time

        self.log("Finished Recieving Scans")
        sink.complete()
        self.log("Completed Sample %s" % (self.sample_name,))
        sink.commit()
