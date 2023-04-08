import os
import logging
import sys
import faulthandler
import multiprocessing

from typing import Callable, List, Optional, Tuple, Union

import click

import ms_deisotope
from ms_deisotope.data_source.infer_type import MSFileLoader

from ms_deisotope.data_source.scan.loader import RandomAccessScanSource, ScanIterator
from ms_deisotope.data_source.scan.scan_iterator import MSEIterator
from ms_deisotope.data_source.scan.mobility_frame import (
    IonMobilityFrame, Generic3DIonMobilityFrameSource, IonMobilitySource, IonMobilitySourceRandomAccessFrameSource)

from ms_deisotope.peak_set import IonMobilityProfileDeconvolutedPeakSolution

from ms_deisotope.output.mzml import MzMLSerializer, IonMobilityAware3DMzMLSerializer, ProcessedGeneric3DIonMobilityFrameSource

from ms_deisotope.feature_map.mobility_frame_processor import IonMobilityFrameProcessor
from ms_deisotope.feature_map.feature_map import IonMobilityProfileDeconvolutedLCMSFeatureForest
from ms_deisotope.feature_map.feature_graph import IonMobilityProfileDeconvolutedFeatureGraph
from ms_deisotope.feature_map.precursor_product_correlation import PrecursorProductCorrelationGraph, NaiveIonMobilityOverlapIterator

from ms_deisotope.tools.deisotoper.process import ScanIDYieldingProcess, ScanBunchLoader, DeconvolutingScanTransformingProcess
from ms_deisotope.tools.deisotoper.scan_generator import ScanGenerator
from ms_deisotope.tools.deisotoper.workflow import SampleConsumer
from ms_deisotope.tools.deisotoper.output import ThreadedMzMLScanStorageHandler

from ms_deisotope.task import TaskBase
from ms_deisotope.tools.utils import processes_option, register_debug_hook, progress
from ms_deisotope.tools.maintenance import register_waters_masslynx


faulthandler.enable()


logger = logging.getLogger("mse_deconvolute")
logger.addHandler(logging.NullHandler())


def make_iterator(reader, start_index, stop_index=float('inf'), low_energy_function=1, lock_mass_function=3):
    iterator = MSEIterator(
        reader.start_from_frame(
            index=start_index, require_ms1=False, grouped=False),
        lambda x: x, low_energy_config=low_energy_function, lock_mass_config=lock_mass_function)
    for bunch in iterator:
        if bunch.precursor:
            i = bunch.precursor.index
        elif bunch.products:
            i = bunch.products[0].index
        else:
            i = 0
        if i >= stop_index:
            break
        yield bunch


def open_mse_file(path, **kwargs):
    reader = MSFileLoader(path, **kwargs)
    if not isinstance(reader, IonMobilitySource):
        reader = Generic3DIonMobilityFrameSource(reader)
    return reader


class MSEFrameIDYieldingProcess(ScanIDYieldingProcess):
    loader: IonMobilitySourceRandomAccessFrameSource

    _iterator: MSEIterator

    low_energy_function: int = 1
    lock_mass_function: int = 3

    def __init__(self, ms_file_path: os.PathLike, scan_id_queue: multiprocessing.JoinableQueue, start_scan: str = None,
                 max_scans: Optional[int] = None, end_scan: str = None, no_more_event: Optional[multiprocessing.Event] = None,
                 ignore_tandem_scans: bool = False, batch_size: int = 1, log_handler: Callable = None,
                 output_queue: Optional[multiprocessing.JoinableQueue] = None, low_energy_function=1, lock_mass_function=3):
        super().__init__(
            ms_file_path=ms_file_path, scan_id_queue=scan_id_queue, start_scan=start_scan, max_scans=max_scans,
            end_scan=end_scan, no_more_event=no_more_event, ignore_tandem_scans=ignore_tandem_scans, batch_size=batch_size,
            log_handler=log_handler, output_queue=output_queue
        )
        self.low_energy_function = low_energy_function
        self.lock_mass_function = lock_mass_function

    def _open_ms_file(self) -> Union[ScanIterator, RandomAccessScanSource]:
        path = self.ms_file_path
        reader = open_mse_file(path)
        self.loader = reader
        return reader

    def _initialize_iterator(self):
        if self.start_scan is not None:
            self.loader.start_from_frame(
                self.start_scan, require_ms1=False, grouped=False)
        else:
            self.loader.reset()
            self.loader.make_frame_iterator(grouped=False)
        self._iterator = MSEIterator(
            self.loader,
            lambda x: x,
            low_energy_config=self.low_energy_function,
            lock_mass_config=self.lock_mass_function,
            on_lock_mass_scan=self.skip_scan)

    def _prepare_end_scan_marker(self) -> Optional[str]:
        end_scan = self.end_scan
        if end_scan is None:
            try:
                self.end_scan_index = len(self.loader)
            except AttributeError:
                self.end_scan_index = sys.maxint
            self.log_handler(
                f"End scan not specified, defaulting to index {self.end_scan_index}")
        else:
            self.end_scan_index = self.loader.get_frame_by_id(
                self.end_scan).index
        return end_scan


class FrameBunchLoader(ScanBunchLoader):
    loader: IonMobilitySourceRandomAccessFrameSource

    def get(self) -> Tuple[IonMobilityFrame, List[IonMobilityFrame]]:
        scan_id, product_scan_ids = self.queue.popleft()
        if scan_id is not None:
            precursor = self.loader.get_frame_by_id(scan_id)
        else:
            precursor = None
        products = [self.loader.get_frame_by_id(pid)
                    for pid in product_scan_ids if pid is not None]
        return (precursor, products)


class MSEDeconvolutingFrameTransformingProcess(DeconvolutingScanTransformingProcess):
    loader: IonMobilitySourceRandomAccessFrameSource
    transformer: IonMobilityFrameProcessor

    _loggers_to_silence = ["ms_deisotope.frame_processor"]

    def __init__(self, ms_file_path, input_queue, output_queue, no_more_event=None, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None, ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 ms1_averaging=0, log_handler=None, deconvolute=True, verbose=False, reader_options=None):
        reader_options = reader_options or {}
        self.reader_options = reader_options
        super().__init__(
            ms_file_path, input_queue, output_queue, no_more_event,
            ms1_peak_picking_args, msn_peak_picking_args,
            ms1_deconvolution_args, msn_deconvolution_args,
            None,
            ms1_averaging=ms1_averaging,
            log_handler=log_handler,
            deconvolute=deconvolute,
            too_many_peaks_threshold=0,
            default_precursor_ion_selection_window=0)

    def make_scan_transformer(self, loader: IonMobilitySourceRandomAccessFrameSource = None) -> IonMobilityFrameProcessor:
        self.transformer = IonMobilityFrameProcessor(
            loader,
            ms1_peak_picking_args=self.ms1_peak_picking_args,
            msn_peak_picking_args=self.msn_peak_picking_args,
            ms1_deconvolution_args=self.ms1_deconvolution_args,
            msn_deconvolution_args=self.msn_deconvolution_args,
            loader_type=lambda x: x,
            ms1_averaging=self.ms1_averaging)
        return self.transformer

    def _process_ms1(self, scan, product_scans) -> Tuple[IonMobilityFrame, List, List[IonMobilityFrame]]:
        scan, priorities, product_scans = self.transformer.process_frame_group(
            scan, product_scans)
        return scan, priorities, product_scans

    def _deconvolute_ms1(self, scan: IonMobilityFrame, priorities: List, product_scans: List[IonMobilityFrame]):
        self.transformer.deconvolute_precursor_features(scan)
        scan.features = None

    def _process_msn(self, product_scan: IonMobilityFrame):
        self.transformer.extract_product_features(product_scan)

    def _deconvolute_msn(self, product_scan: IonMobilityFrame):
        self.transformer.deconvolute_product_features(product_scan)
        product_scan.features = None

    def _open_ms_file(self) -> IonMobilitySourceRandomAccessFrameSource:
        self.loader = open_mse_file(self.ms_file_path, **self.reader_options)
        return self.loader

    def _make_batch_loader(self, loader: IonMobilitySourceRandomAccessFrameSource) -> FrameBunchLoader:
        return FrameBunchLoader(loader)

    def send_scan(self, scan: IonMobilityFrame):
        if scan.deconvoluted_features is not None:
            scan.features = None
        return super().send_scan(scan)


class MSEFrameGenerator(ScanGenerator):
    def __init__(self, ms_file, number_of_helpers=4,
                 ms1_peak_picking_args=None, msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 ms1_averaging=0, deconvolute=True, verbose=False,
                 reader_options=None, low_energy_function=1, lock_mass_function=3):
        if reader_options is None:
            reader_options = dict()
        self.ms_file = ms_file

        self.scan_ids_exhausted_event = multiprocessing.Event()
        self.reader_options = reader_options
        self._iterator = None

        self._scan_yielder_process = None
        self._deconv_process = None

        self._input_queue = None
        self._output_queue = None
        self._deconv_helpers = None
        self._order_manager = None

        self.number_of_helpers = number_of_helpers

        self.ms1_peak_picking_args = ms1_peak_picking_args
        self.msn_peak_picking_args = msn_peak_picking_args
        self.ms1_averaging = ms1_averaging

        self.deconvoluting = deconvolute
        self.ms1_deconvolution_args = ms1_deconvolution_args
        self.msn_deconvolution_args = msn_deconvolution_args

        self.low_energy_function = low_energy_function
        self.lock_mass_function = lock_mass_function

        self.extract_only_tandem_envelopes = False
        self.default_precursor_ion_selection_window = 0
        self.ignore_tandem_scans = False

        self._scan_interval_tree = None
        self.verbose = verbose
        self.log_controller = self.ipc_logger()

    def _open_ms_file(self, **kwargs):
        return open_mse_file(self.ms_file, **kwargs)

    def _make_scan_id_yielder(self, start_scan: str, end_scan: str, max_scans: int) -> MSEFrameIDYieldingProcess:
        return MSEFrameIDYieldingProcess(
            self.ms_file, self._input_queue, start_scan=start_scan, end_scan=end_scan,
            max_scans=max_scans, no_more_event=self.scan_ids_exhausted_event,
            ignore_tandem_scans=self.ignore_tandem_scans, batch_size=1,
            output_queue=self._output_queue,
            low_energy_function=self.low_energy_function,
            lock_mass_function=self.lock_mass_function)

    def _make_transforming_process(self) -> MSEDeconvolutingFrameTransformingProcess:
        return MSEDeconvolutingFrameTransformingProcess(
            self.ms_file,
            self._input_queue,
            self._output_queue,
            self.scan_ids_exhausted_event,
            ms1_peak_picking_args=self.ms1_peak_picking_args,
            msn_peak_picking_args=self.msn_peak_picking_args,
            ms1_deconvolution_args=self.ms1_deconvolution_args,
            msn_deconvolution_args=self.msn_deconvolution_args,
            log_handler=self.log_controller.sender(),
            ms1_averaging=self.ms1_averaging,
            deconvolute=self.deconvoluting,
            verbose=self.verbose,
            reader_options=self.reader_options)


class MSESampleConsumer(SampleConsumer):
    def __init__(self, ms_file,
                 ms1_peak_picking_args=None, msn_peak_picking_args=None, ms1_deconvolution_args=None,
                 msn_deconvolution_args=None, start_scan_id=None, end_scan_id=None, storage_path=None,
                 sample_name=None, storage_type=None, n_processes=5,
                 ms1_averaging=0,
                 deconvolute=True,
                 verbose=False,
                 start_scan_time=None,
                 end_scan_time=None,
                 reader_options=None,
                 low_energy_function=1,
                 lock_mass_function=3):

        if storage_type is None:
            storage_type = IonMobilityAware3DThreadedMzMLScanStorageHandler

        self.ms_file = ms_file
        self.storage_path = storage_path
        self.sample_name = sample_name

        self.n_processes = n_processes
        self.storage_type = storage_type
        self.ms1_averaging = ms1_averaging
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
        self.scan_generator = MSEFrameGenerator(
            ms_file,
            number_of_helpers=n_helpers,
            ms1_peak_picking_args=ms1_peak_picking_args,
            msn_peak_picking_args=msn_peak_picking_args,
            ms1_deconvolution_args=ms1_deconvolution_args,
            msn_deconvolution_args=msn_deconvolution_args,
            ms1_averaging=ms1_averaging,
            deconvolute=deconvolute,
            verbose=verbose,
            reader_options=reader_options,
            low_energy_function=low_energy_function,
            lock_mass_function=lock_mass_function)

        self.start_scan_id = start_scan_id
        self.end_scan_id = end_scan_id
        self.start_scan_time = start_scan_time
        self.end_scan_time = end_scan_time

        self.sample_run = None


class IonMobilityAwareMzMLSerializer(MzMLSerializer):
    def _prepare_extra_arrays(self, scan, **kwargs):
        extra_arrays = super(IonMobilityAwareMzMLSerializer,
                             self)._prepare_extra_arrays(scan, **kwargs)
        if scan.deconvoluted_peak_set is not None:
            # This is sensitive to units used? Shouldn't there be a unit key?
            # Waters uses milliseconds
            extra_arrays.append(("mean drift time array", [
                p.drift_time for p in scan.deconvoluted_peak_set
            ]))
        return extra_arrays


class IonMobilityAware3DThreadedMzMLScanStorageHandler(ThreadedMzMLScanStorageHandler):
    def _make_writer(self, n_spectra: int, sample_name: str, deconvoluted: bool, stream_cls):
        self.handle = stream_cls(self.path, 'wb')
        serializer = IonMobilityAware3DMzMLSerializer(
            self.handle,
            n_spectra=n_spectra,
            sample_name=sample_name,
            deconvoluted=True)
        return serializer


averagine_map = {
    "glycopeptide": ms_deisotope.glycopeptide,
    "heparin": ms_deisotope.heparin,
    "peptide": ms_deisotope.peptide,
    "glycan": ms_deisotope.glycan,
    "heparan_sulfate": ms_deisotope.heparan_sulfate,
    "permethylated_glycan": ms_deisotope.permethylated_glycan,
}


def _default_log_handler():
    formatter = logging.Formatter(
        '%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(formatter)
    return handler



@click.group("waters-cyclic-deconvolute")
def cli():
    pass


@cli.command("feature-deconvolution", short_help="Deisotope and charge state deconvolute LC-IM-MSe run")
@click.argument("input_path", type=click.Path())
@click.argument("output_path", type=click.Path(writable=True))
@click.option("-m", "--lockmass-config", type=float, help="The lock mass used")
@click.option("-s", "--start-time", type=float, help="The time to start processing cycles from", default=0)
@click.option("-e", "--end-time", type=float, help="The time to stop processing cycles at", default=None)
@click.option("-a", "--averagine", type=click.Choice(list(averagine_map)), default='glycopeptide',
              help='The isotopic model to use. Defaults to the glycopeptide averagine.')
@click.option("-i", "--minimum-intensity", type=float, default=10.0, help="The minimum intensity to accept a peak")
@click.option("-w", "--isolation-window-width", type=float, default=0.0,
              help="The isolation window size on either side of the set mass.")
@processes_option
@click.option("-k", "--lock-mass-function", type=int, default=3, help="The number of the lock mass function. For normal low-high MSE this is 3.")
@click.option("-d", "--denoise", type=float, default=1.0, help="Aggressiveness of background noise removal. Defaults to 1.0")
@click.option("-g", "--signal-averaging", type=int, default=2, help="The number of adjacent IM-MS spectra to average within a single cycle.")
def feature_deconvolution(input_path, output_path, lockmass_config=None, start_time=0, end_time=None, averagine='glycopeptide',
                          minimum_intensity=10.0, lock_mass_function=3, processes: int = 4,
                          isolation_window_width=0.0, denoise=1.0, signal_averaging=2):
    """Extract features from each IM-MS cycle followed by deisotoping and charge state deconvolution.
    """
    sample_name = os.path.basename(input_path).rsplit(".", 1)[0]
    logging.basicConfig(
        level="INFO", format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
        filemode='w',
        filename="cyclic_deconvolute_%s_%s.log" % (sample_name, start_time))
    logging.getLogger().addHandler(_default_log_handler())
    input_path = str(input_path)

    print(f"Running on PID {os.getpid()}")

    open_options = {}
    if lockmass_config is not None:
        open_options['lockmass_config'] = lockmass_config

    reader = open_mse_file(input_path, **open_options)

    if start_time is not None:
        start_id = reader.get_frame_by_time(start_time).id
    else:
        start_id = None
    if end_time is not None:
        end_id = reader.get_frame_by_time(end_time).id
    else:
        end_id = None

    averagine = averagine_map[averagine]

    task = MSESampleConsumer(
        input_path, storage_path=output_path,
        ms1_peak_picking_args={"error_tolerance": 4e-5,
                               "minimum_intensity": minimum_intensity / 2, "denoise": denoise},
        msn_peak_picking_args={
            "average_within": signal_averaging, "error_tolerance": 4e-5, "denoise": denoise},
        ms1_deconvolution_args={
            "averagine": averagine,
            "truncate_after": 0.95,
            "scorer": ms_deisotope.PenalizedMSDeconVFitter(5, 1),
            "minimum_intensity": minimum_intensity,
            "copy": False
        },
        msn_deconvolution_args={
            "averagine": averagine,
            "truncate_after": 0.8,
            "scorer": ms_deisotope.MSDeconVFitter(1),
            "minimum_intensity": minimum_intensity / 2,
            "copy": False
        }, ms1_averaging=signal_averaging,
        reader_options={"lockmass_config": lockmass_config,
                        "default_isolation_width": isolation_window_width},
        sample_name=sample_name,
        deconvolute=True,
        n_processes=processes,
        start_scan_id=start_id,
        end_scan_id=end_id,
        start_scan_time=start_time,
        end_scan_time=end_time,
        lock_mass_function=lock_mass_function)

    task.start()


@cli.command("deconvolve-products", short_help="Deconvolve correlated precursor-product ion relationships")
@click.argument("input_path", type=click.Path())
@click.argument("output_path", type=click.Path(writable=True))
@click.option("-e", "--edges-per-feature", type=int, default=1000,
              help='The maximum number of precursor features a product feature can connect to.')
@click.option("-w", "--weight-scaling", type=float, default=1.0,
              help="A weight scaling factor. >1 to incrase the overall abundance of deconvolved peaks")
def precursor_product_deconvolution(input_path, output_path, edges_per_feature: int=1000, weight_scaling: float=1.0):
    """Takes a deconvolved LC-IM-MSe run and generate pseudospectra for
    precursor ions using correlated product ion features enclosed in the IM
    and RT dimensions.
    """
    logging.basicConfig(
        level="INFO", format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
        filemode='w',
        filename="precursor_product_deconvolution_%s.log" % (
            os.path.basename(input_path).rsplit(".", 1)[0]))

    print(f"Running on PID {os.getpid()}")

    logging.getLogger().addHandler(_default_log_handler())
    input_path = str(input_path)

    scan_reader = ms_deisotope.MSFileLoader(input_path)
    frame_reader = ProcessedGeneric3DIonMobilityFrameSource(scan_reader)


    precursor_forest = IonMobilityProfileDeconvolutedLCMSFeatureForest()
    product_forest = IonMobilityProfileDeconvolutedLCMSFeatureForest()

    frame_reader.make_frame_iterator(grouped='mse')
    logger.info("Constructing LC-IM-MSe features")
    for i, bunch in enumerate(frame_reader):
        if bunch.precursor is None:
            continue
        if i % 100 == 0:
            logger.info("... Processing %r %0.3f", bunch.precursor.id, bunch.precursor.time)
        if bunch.precursor:
            for f in bunch.precursor.deconvoluted_features:
                p = IonMobilityProfileDeconvolutedPeakSolution.from_feature(f)
                precursor_forest.handle_peak(p, bunch.precursor.time)

        for product in bunch.products:
            for f in product.deconvoluted_features:
                p = IonMobilityProfileDeconvolutedPeakSolution.from_feature(f)
                product_forest.handle_peak(p, product.time)

    logger.info("Smoothing precursors")
    precursor_forest.smooth_overlaps().split_sparse(0.5)
    logger.info("Smoothing products")
    product_forest.smooth_overlaps().split_sparse(0.5)

    logger.info("Building precursor-product graph")
    precursor_graph = IonMobilityProfileDeconvolutedFeatureGraph(precursor_forest)
    product_graph = IonMobilityProfileDeconvolutedFeatureGraph(product_forest)
    corr_graph = PrecursorProductCorrelationGraph(
        precursor_graph, product_graph, max_edge_count_per_node=edges_per_feature)
    corr_graph.build()

    logger.info("Generating pseudo-spectra")
    fh = open(output_path, 'wb')
    with MzMLSerializer(fh, len(scan_reader), sample_name=os.path.basename(output_path)) as writer:
        writer.copy_metadata_from(scan_reader)
        proc_method = writer.build_processing_method()
        writer.add_data_processing(proc_method)
        for bunch in corr_graph.iterspectra(weight_scale_factor=weight_scaling):
            writer.save(bunch)


@cli.command("ion-mobility-overlap-pseudospectra", short_help="Generate naive pseudo-spectra from ion mobility overlaps")
@click.argument("input_path", type=click.Path())
@click.argument("output_path", type=click.Path(writable=True))
def naive_ion_mobility_overlap_pseudospectra(input_path, output_path):
    """Generate pseudo-spectra from deconvolved cycles where the precursor ion spans
    some or all of the ion dimension of the product ion's mobility dimension and has
    a larger neutral mass. Makes no use of retention time.
    """
    logging.basicConfig(
        level="INFO", format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
        filemode='w',
        filename="naive_ion_mobility_overlap_pseudospectra_%s.log" % (
            os.path.basename(input_path).rsplit(".", 1)[0]))

    print(f"Running on PID {os.getpid()}")

    logging.getLogger().addHandler(_default_log_handler())
    input_path = str(input_path)

    scan_reader = ms_deisotope.MSFileLoader(input_path)
    frame_reader = ProcessedGeneric3DIonMobilityFrameSource(scan_reader)
    frame_reader.make_frame_iterator(grouped='mse')

    n_frames = len(scan_reader)

    iterator = NaiveIonMobilityOverlapIterator(frame_reader)

    logger.info("Generating pseudo-spectra")
    fh = open(output_path, 'wb')
    with click.progressbar(length=n_frames, item_show_func=lambda x: x if x else "-") as prog:
        with MzMLSerializer(fh, len(scan_reader) * 50, sample_name=os.path.basename(output_path)) as writer:
            writer.copy_metadata_from(scan_reader)
            proc_method = writer.build_processing_method()
            writer.add_data_processing(proc_method)
            for bunch in iterator:
                prog.update(n_steps=2, current_item=bunch.precursor.id)
                writer.save(bunch)


cli.add_command(register_waters_masslynx)


if __name__ == "__main__":
    import multiprocessing
    register_debug_hook()
    multiprocessing.freeze_support()
    logging.captureWarnings(True)
    TaskBase.log_with_logger(logger)
    cli.main()
