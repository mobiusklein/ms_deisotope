import os
import logging
import sys
import faulthandler
import multiprocessing
import pickle

from typing import List, Optional, Tuple, Union


from itertools import chain

import psutil
import numpy as np

import ms_deisotope
from ms_deisotope.data_source.infer_type import MSFileLoader
from ms_deisotope.data_source.scan.loader import RandomAccessScanSource, ScanIterator
from ms_deisotope.peak_dependency_network.intervals import Interval, IntervalTreeNode
from ms_deisotope.peak_set import IonMobilityDeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.data_source.scan import ProcessedScan, PrecursorInformation, ScanBunch
from ms_deisotope.data_source.scan.scan_iterator import MSEIterator
from ms_deisotope.output.mzml import MzMLSerializer, IonMobilityAware3DMzMLSerializer
from ms_deisotope.data_source.scan.mobility_frame import (
    IonMobilityFrame, Generic3DIonMobilityFrameSource, IonMobilitySource, IonMobilitySourceRandomAccessFrameSource)
from ms_deisotope.feature_map.mobility_frame_processor import IonMobilityFrameProcessor

from ms_deisotope.tools.deisotoper.process import ScanIDYieldingProcess, ScanBunchLoader, DeconvolutingScanTransformingProcess
from ms_deisotope.tools.deisotoper.scan_generator import ScanGenerator
from ms_deisotope.tools.deisotoper.workflow import ScanSink, SampleConsumer
from ms_deisotope.tools.deisotoper.output import ThreadedMzMLScanStorageHandler

from ms_deisotope.tools.utils import register_debug_hook


import click


faulthandler.enable()


logger = logging.getLogger("mse_deconvolute")
logger.addHandler(logging.NullHandler())


def make_iterator(reader, start_index, stop_index=float('inf'), low_energy_function=1, lock_mass_function=3):
    iterator = MSEIterator(
        reader.start_from_frame(
            index=start_index, require_ms1=False, grouped=False),
        lambda x: x, low_energy_function, lock_mass_function)
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

    _loggers_to_silence = ["deconvolution_frame_processor"]

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


class MSEFrameGenerator(ScanGenerator):
    def __init__(self, ms_file, number_of_helpers=4,
                 ms1_peak_picking_args=None, msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 ms1_averaging=0, deconvolute=True, verbose=False, reader_options=None):
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
            output_queue=self._output_queue)

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
                 deconvolute=True, verbose=False, start_scan_time=None, end_scan_time=None, reader_options=None):

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
            reader_options=reader_options)

        self.start_scan_id = start_scan_id
        self.end_scan_id = end_scan_id
        self.start_scan_time = start_scan_time
        self.end_scan_time = end_scan_time

        self.sample_run = None


def connect_products_to_precursors(parent, child):
    cft = IntervalTreeNode.build(
        [Interval(f.start_time, f.end_time, [f]) for f in child.deconvoluted_features])

    prec_prod_rels = []
    for f in parent.deconvoluted_features:
        if len(f) < 3:
            continue
        prods = []
        for prod in (chain.from_iterable(cft.overlaps(f.start_time, f.end_time))):
            if prod.charge > f.charge:
                continue
            if prod.neutral_mass > f.neutral_mass:
                continue
            chunk = prod._copy_chunk(
                prod[prod.find_time(f.start_time)[1]:
                     prod.find_time(f.end_time)[1] + 1])
            prods.append(chunk)
        prods.sort(key=lambda x: x.neutral_mass)
        prec_prod_rels.append((f, prods))
    prec_prod_rels.sort(key=lambda x: (
        x[0].start_time, x[0].neutral_mass, len(x[1])))
    return prec_prod_rels


def weighted_centroid(feature):
    total = 0
    normalizer = 0
    for node in feature:
        weight = node.total_intensity()
        total += node.time * weight
        normalizer += weight
    return total / normalizer


def merge_envelopes(envelopes):
    base = envelopes[0].clone()
    for env in envelopes[1:]:
        for i, p in enumerate(env):
            base[i].intensity += p.intensity

    return base


def feature_to_peak(feature):
    peak_cluster = feature.peaks
    peak_cluster = [pi for p in peak_cluster for pi in p]
    total_intensity = sum(p.intensity for p in peak_cluster)
    mz = sum(p.mz * p.intensity for p in peak_cluster) / total_intensity
    neutral_mass = sum(
        p.neutral_mass * p.intensity for p in peak_cluster) / total_intensity
    most_abundant_mass = sum(
        p.most_abundant_mass * p.intensity for p in peak_cluster) / total_intensity
    a_to_a2_ratio = sum(
        p.a_to_a2_ratio * p.intensity for p in peak_cluster) / total_intensity
    average_mass = sum(
        p.average_mass * p.intensity for p in peak_cluster) / total_intensity
    signal_to_noise = sum(p.signal_to_noise *
                          p.intensity for p in peak_cluster) / total_intensity
    fwhm = sum(p.full_width_at_half_max *
               p.intensity for p in peak_cluster) / total_intensity
    area = sum(p.area * p.intensity for p in peak_cluster) / total_intensity
    score = sum(p.score * p.intensity for p in peak_cluster) / total_intensity
    charge = peak_cluster[0].charge
    envelope = merge_envelopes([p.envelope for p in peak_cluster])
    return IonMobilityDeconvolutedPeak(
        neutral_mass=neutral_mass, intensity=total_intensity, charge=charge, signal_to_noise=signal_to_noise,
        full_width_at_half_max=fwhm, index=-1, a_to_a2_ratio=a_to_a2_ratio, most_abundant_mass=most_abundant_mass,
        average_mass=average_mass, score=score, envelope=envelope, mz=mz, fit=None, chosen_for_msms=False,
        area=area, drift_time=weighted_centroid(feature))


def frame_pair_to_scan_bunch(parent, child, prec_prod_rels, scan_index):
    precursor = ProcessedScan(
        parent.id, parent.id, None, 1, parent.time, scan_index, None,
        DeconvolutedPeakSet([feature_to_peak(f)
                             for f in parent.deconvoluted_features]).reindex(),
        # Hard-coded positive mode for the moment while IonMobilityFrame-Source-like does not have a polarity
        1,
    )

    scan_index += 1

    pseudo_ms2s = []

    for prec_feat, prod_feats in prec_prod_rels:
        if abs(prec_feat.charge) == 1:
            continue
        pinfo = PrecursorInformation(
            prec_feat.mz, prec_feat.total_signal, prec_feat.charge, precursor.id, None, annotations={
                "ion mobility drift time": weighted_centroid(prec_feat),
                "source frame id": child.id,
            })
        assert pinfo.has_ion_mobility()
        # Eventually, create an artificial ScanAcquisitionInformation instance here to let us set the
        # drift time of the product scan to match the precursor drift time too?
        prod = ProcessedScan(
            'merged=%d' % scan_index, child.id + '.%d' % scan_index, pinfo, 2, child.time,
            scan_index, None,
            DeconvolutedPeakSet([feature_to_peak(f)
                                 for f in prod_feats]).reindex(),
            # Hard-coded positive mode for the moment while IonMobilityFrame-Source-like does not have a polarity
            1, activation=child.activation)
        prod.annotations['spectrum title'] = child.id + '.%d' % scan_index
        pseudo_ms2s.append(prod)
        scan_index += 1

    return ScanBunch(precursor, pseudo_ms2s), scan_index


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


@click.command("cyclic_deconvolute")
@click.argument("input_path", type=click.Path())
@click.argument("output_path", type=click.Path(writable=True))
@click.option("-m", "--lockmass-config", type=float, help="The lock mass used", default=785.8421)
@click.option("-s", "--start-time", type=float, help="The time to start processing cycles from", default=0)
@click.option("-e", "--end-time", type=float, help="The time to stop processing cycles at", default=None)
@click.option("-a", "--averagine", type=click.Choice(list(averagine_map)), default='glycopeptide',
              help='The isotopic model to use. Defaults to the glycopeptide averagine.')
@click.option("-i", "--minimum-intensity", type=float, default=10.0, help="The minimum intensity to accept a peak")
@click.option("-n", "--no-product-splitting", is_flag=True, default=False,
              help=(
                  'If true, do not split high energy frames into pseudo-MS2 spectra,'
                  ' creating continuous ion mobility features'))
@click.option("-k", "--lockmass-function", type=int, default=3, help="The number of the lock mass function. For normal low-high MSE this is 3.")
def main(input_path, output_path, lockmass_config, start_time=0, end_time=None, averagine='glycopeptide',
         minimum_intensity=10.0, no_product_splitting=False, lockmass_function=3):
    logging.basicConfig(
        level="INFO", format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
        filemode='w',
        filename="cyclic_deconvolute_%s_%s.log" % (os.path.basename(input_path).rsplit(".", 1)[0], start_time))
    logging.getLogger().addHandler(logging.StreamHandler(sys.stderr))
    input_path = str(input_path)

    print(os.getpid())

    reader = open_mse_file(
        input_path, lockmass_config=lockmass_config)

    if start_time is not None:
        start_id = reader.get_frame_by_time(start_time).id
    else:
        start_id = None
    if end_time is not None:
        end_id = reader.get_frame_by_time(end_time).id
    else:
        end_id = None

    task = MSESampleConsumer(
        input_path, storage_path=output_path, ms1_peak_picking_args={}, msn_peak_picking_args={},
        ms1_deconvolution_args={
            "averagine": ms_deisotope.glycopeptide,
            "truncate_after": 0.95,
            "scorer": ms_deisotope.PenalizedMSDeconVFitter(5, 2),
            "minimum_intensity": minimum_intensity,
            "copy": False
        },
        msn_deconvolution_args={
            "averagine": ms_deisotope.glycopeptide,
            "truncate_after": 0.95,
            "scorer": ms_deisotope.MSDeconVFitter(1),
            "minimum_intensity": minimum_intensity,
            "copy": False
        }, ms1_averaging=1, reader_options={"lockmass_config": lockmass_config},
        deconvolute=True,
        n_processes=8,
        start_scan_id=start_id,
        end_scan_id=end_id,
        start_scan_time=start_time,
        end_scan_time=end_time)

    task.start()
