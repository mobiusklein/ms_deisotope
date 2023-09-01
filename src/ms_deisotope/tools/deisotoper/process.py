import os
import signal
import logging
import sys
import multiprocessing
import traceback
import pickle

from typing import Any, Callable, Deque, Dict, List, Optional, Tuple, Union, TYPE_CHECKING
from collections import deque
from multiprocessing import Process
from queue import Empty as QueueEmpty

from ms_deisotope.data_source import Scan, RandomAccessScanSource, ScanIterator


import ms_deisotope
from ms_deisotope.data_source.scan.scan import ProcessedScan

from ms_deisotope.processor import (
    ScanProcessor, MSFileLoader,
    NoIsotopicClustersError, EmptyScanError)

from ms_deisotope.task import show_message
from ms_deisotope.output.mzml import _PeakPacker

try:
    from pyzstd import decompress, compress
except ImportError:
    from gzip import decompress, compress

if TYPE_CHECKING:
    import multiprocessing.synchronize

class CompressedPickleMessage(object):
    def __init__(self, obj=None):
        self.obj = obj

    def __getstate__(self):
        state = {
            "payload": compress(pickle.dumps(self.obj, -1))
        }
        return state

    def __setstate__(self, state):
        self.obj = pickle.loads(decompress(state['payload']))

    def __reduce__(self):
        return self.__class__, (None, ), self.__getstate__()

    def __iter__(self):
        yield self.obj


DONE = b"--NO-MORE--"
SCAN_STATUS_GOOD = b"good"
SCAN_STATUS_SKIP = b"skip"


class ScanTransmissionMixin(object):
    output_queue: multiprocessing.JoinableQueue
    scan_packer: Optional[_PeakPacker] = None

    def skip_entry(self, index: int, ms_level: int):
        self.output_queue.put((SCAN_STATUS_SKIP, index, ms_level))

    def skip_scan(self, scan: Scan):
        self.output_queue.put((SCAN_STATUS_SKIP, scan.index, scan.ms_level))

    def send_scan(self, scan: Scan):
        scan = scan.pack()
        # this attribute is not needed, and for MS1 scans is dangerous
        # to pickle.
        # It can pull other scans which may not yet have been packed
        # into the message sent back to the main process which in
        # turn can form a reference cycle and eat a lot of memory
        scan.product_scans = []
        if isinstance(scan, ProcessedScan):
            if scan.deconvoluted_peak_set is not None:
                scan.peak_set = None
                for peak in scan.deconvoluted_peak_set:
                    if peak.fit is not None:
                        peak.fit = None
        if self.scan_packer is not None:
            scan = self.scan_packer(scan)
        self.output_queue.put(
            (CompressedPickleMessage(scan), scan.index, scan.ms_level))


class _ProcessHelper:
    _error_occurred: "multiprocessing.synchronize.Event"

    def try_set_process_name(self, name=None):
        """
        This helper method may be used to try to change a process's name
        in order to make discriminating which role a particular process is
        fulfilling. This uses a third-party utility library that may not behave
        the same way on all platforms, and therefore this is done for convenience
        only.

        Parameters
        ----------
        name : str, optional
            A name to set. If not provided, will check the attribute ``process_name``
            for a non-null value, or else have no effect.
        """
        if name is None:
            name = getattr(self, 'process_name', None)
        if name is None:
            return
        try:
            import setproctitle
            setproctitle.setproctitle(name)
        except (ImportError, AttributeError):
            pass

    def error_occurred(self) -> bool:
        return self._error_occurred.is_set()


class ScanIDYieldingProcess(Process, ScanTransmissionMixin, _ProcessHelper):
    ms_file_path: os.PathLike
    scan_id_queue: multiprocessing.JoinableQueue
    loader: Union[ScanIterator, RandomAccessScanSource]

    start_scan: str
    end_scan: str
    max_scans: Optional[int]
    end_scan_index: Optional[int]

    ignore_tandem_scans: bool
    batch_size: int

    no_more_event: Optional[multiprocessing.Event]
    log_handler: Callable

    def __init__(self, ms_file_path: os.PathLike, scan_id_queue: multiprocessing.JoinableQueue, start_scan: str = None,
                 max_scans: Optional[int]=None, end_scan: str=None, no_more_event: Optional[multiprocessing.Event]=None,
                 ignore_tandem_scans: bool=False, batch_size: int=1, log_handler: Callable=None,
                 output_queue: Optional[multiprocessing.JoinableQueue]=None):
        if log_handler is None:
            log_handler = show_message
        Process.__init__(self)
        self.daemon = True
        self.ms_file_path = ms_file_path
        self.scan_id_queue = scan_id_queue
        self.loader = None

        self.start_scan = start_scan
        self.max_scans = max_scans
        self.end_scan = end_scan
        self.end_scan_index = None
        self.passed_first_batch = False
        self.ignore_tandem_scans = ignore_tandem_scans
        self.batch_size = batch_size

        self.log_handler = log_handler

        self.no_more_event = no_more_event
        self.output_queue = output_queue
        self._error_occurred = multiprocessing.Event()

    def _make_scan_batch(self) -> Tuple[
            List[Tuple[str, List[str]]],
            List[str]
        ]:
        batch = []
        scan_ids = []
        for _ in range(self.batch_size):
            try:
                bunch = next(self._iterator)
                scan, products = bunch
                if scan is not None:
                    scan_id = scan.id
                    if scan.index > self.end_scan_index:
                        break
                else:
                    scan_id = None
                products = [prod for prod in products if prod.index <= self.end_scan_index]
                product_scan_ids = [p.id for p in products]
            except StopIteration:
                break
            except Exception as e:
                self.log_handler("An error occurred in _make_scan_batch", e)
                break
            if not self.ignore_tandem_scans:
                batch.append((scan_id, product_scan_ids, True))
            else:
                batch.append((scan_id, product_scan_ids, False))
            scan_ids.append(scan_id)
        return batch, scan_ids

    def _initialize_iterator(self):
        if self.start_scan is not None:
            try:
                self.loader.start_from_scan(
                    self.start_scan, require_ms1=self.loader.has_ms1_scans(), grouped=True)
            except IndexError as e:
                self.log_handler(
                    "An error occurred while locating start scan", e)
                self.loader.reset()
                self.loader.make_iterator(grouped=True)
            except AttributeError as e:
                self.log_handler(
                    "The reader does not support random access, start time will be ignored", e)
                self.loader.reset()
                self.loader.make_iterator(grouped=True)
        else:
            self.loader.make_iterator(grouped=True)
        self._iterator = self.loader

    def _prepare_end_scan_marker(self) -> Optional[str]:
        end_scan = self.end_scan
        if end_scan is None:
            try:
                self.end_scan_index = len(self.loader)
            except AttributeError:
                self.end_scan_index = sys.maxint
        else:
            self.end_scan_index = self.loader.get_scan_by_id(end_scan).index
        return end_scan

    def _open_ms_file(self) -> Union[ScanIterator, RandomAccessScanSource]:
        self.loader = MSFileLoader(self.ms_file_path, decode_binary=False)
        return self.loader

    def run(self):
        signal.signal(signal.SIGINT, signal.Handlers.SIG_IGN)
        self.try_set_process_name("ms-deisotope-sched")
        self._open_ms_file()
        self._initialize_iterator()

        count: int = 0
        last: int = 0
        if self.max_scans is None:
            max_scans = float('inf')
        else:
            max_scans = self.max_scans

        end_scan = self._prepare_end_scan_marker()
        while count < max_scans:
            try:
                batch, ids = self._make_scan_batch()
                if len(batch) > 0:
                    self.scan_id_queue.put(batch)
                count += len(ids)
                if (count - last) > 1000:
                    last = count
                    self.scan_id_queue.join()
                if (end_scan in ids and end_scan is not None) or len(ids) == 0:
                    self.log_handler("End Scan Found")
                    break
            except StopIteration:
                break
            except Exception as e:
                self.log_handler("An error occurred while fetching scans", e)
                break

        if self.no_more_event is not None:
            self.no_more_event.set()
            self.log_handler("All Scan IDs have been dealt. %d scan bunches." % (count,))
        else:
            self.scan_id_queue.put(DONE)


class ScanTransformMixin(object):
    _batch_store: Deque[Tuple[str, List[str], bool]]

    input_queue: multiprocessing.JoinableQueue

    def log_error(self, error: Exception, scan_id: str, scan: Scan, product_scan_ids: List[str]):
        tb = traceback.format_exc()
        self.log_handler(
            "An %r occurred for %s (index %r) in Process %r\n%s" % (
                error, scan_id, scan.index, multiprocessing.current_process(),
                tb))

    def _init_batch_store(self):
        self._batch_store = deque()

    def get_work(self, block: bool=True, timeout: float=30) -> Tuple[str, List[str], bool]:
        if self._batch_store:
            return self._batch_store.popleft()
        else:
            batch = self.input_queue.get(block, timeout)
            self.input_queue.task_done()
            self._batch_store.extend(batch)
            result = self._batch_store.popleft()
            return result

    def log_message(self, message):
        self.log_handler(message + ", %r" %
                         (multiprocessing.current_process().name))

    def all_work_done(self) -> bool:
        return self._work_complete.is_set()

    def make_scan_transformer(self, loader=None):
        raise NotImplementedError()

    _loggers_to_silence = ["ms_deisotope.scan_processor"]

    def _silence_loggers(self):
        nologs = self._loggers_to_silence
        if not self.deconvolute:
            nologs.append("deconvolution")

        debug_mode = os.getenv("MS_DEISOTOPE_DEBUG")
        if debug_mode:
            handler = logging.FileHandler(
                "ms-deisotope-deconvolution-debug-%s.log" % (os.getpid()), 'w')
            fmt = logging.Formatter(
                "%(asctime)s - %(name)s:%(filename)s:%(lineno)-4d - %(levelname)s - %(message)s",
                "%H:%M:%S")
            handler.setFormatter(fmt)
        for logname in nologs:
            logger_to_silence = logging.getLogger(logname)
            if debug_mode:
                logger_to_silence.setLevel("DEBUG")
                logger_to_silence.addHandler(handler)
            else:
                logger_to_silence.propagate = False
                logger_to_silence.setLevel("CRITICAL")
                logger_to_silence.addHandler(logging.NullHandler())


class ScanBunchLoader(object):
    queue: Deque[Tuple[str, List[str]]]
    loader: Union[RandomAccessScanSource, ScanIterator]

    def __init__(self, mzml_loader):
        self.loader = mzml_loader
        self.queue = deque()

    def put(self, scan_id: str, product_scan_ids: List[str]):
        self.queue.append((scan_id, product_scan_ids))

    def get(self) -> Tuple[Scan, List[Scan]]:
        scan_id, product_scan_ids = self.queue.popleft()
        if scan_id is not None:
            precursor = self.loader.get_scan_by_id(scan_id)
        else:
            precursor = None
        products = [self.loader.get_scan_by_id(
            pid) for pid in product_scan_ids if pid is not None]
        if precursor:
            precursor.product_scans = products
        return (precursor, products)


class DeconvolutingScanTransformingProcess(Process, ScanTransformMixin, ScanTransmissionMixin, _ProcessHelper):
    """
    DeconvolutingScanTransformingProcess describes a child process that consumes scan id bunches
    from a shared input queue, retrieves the relevant scans, and preprocesses them using an
    instance of :class:`ms_deisotope.processor.ScanProcessor`, sending the reduced result
    to a shared output queue.

    Attributes
    ----------
    input_queue : multiprocessing.JoinableQueue
        A shared input queue which contains payloads of bunches of
        scan ids
    ms1_deconvolution_args : dict
        Parameters passed to :class:`ms_deisotope.processor.ScanProcessor`
    ms1_peak_picking_args : dict
        Parameters passed to :class:`ms_deisotope.processor.ScanProcessor`
    msn_deconvolution_args : dict
        Parameters passed to :class:`ms_deisotope.processor.ScanProcessor`
    msn_peak_picking_args : dict
        Parameters passed to :class:`ms_deisotope.processor.ScanProcessor`
    ms_file_path : str
        Path to the spectral data file on disk
    no_more_event : multiprocessing.Event
        An event which will be set when the process feeding the input
        queue has run out of items to add, indicating that any QueueEmptyException
        should be treated as a signal to finish rather than to wait for
        new input
    output_queue : multiprocessing.JoinableQueue
        A shared output queue which this object will put
        :class:`ms_deisotope.data_source.common.ProcessedScan` bunches onto.
    """

    ms_file_path: os.PathLike
    input_queue: multiprocessing.Queue
    output_queue: multiprocessing.Queue

    loader: Union[ScanIterator, RandomAccessScanSource]

    ms1_averaging: int
    envelope_selector: Callable
    transformer: Optional[ScanProcessor]

    ms1_peak_picking_args: Dict[str, Any]
    msn_peak_picking_args: Dict[str, Any]

    ms1_deconvolution_args: Dict[str, Any]
    msn_deconvolution_args: Dict[str, Any]

    too_many_peaks_threshold: int
    default_precursor_ion_selection_window: float
    deconvolute: bool

    _work_complete: multiprocessing.Event
    no_more_event: Optional[multiprocessing.Event]

    log_handler: Callable

    def __init__(self, ms_file_path, input_queue, output_queue,
                 no_more_event=None, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 envelope_selector=None, ms1_averaging=0, log_handler=None,
                 deconvolute=True, verbose=False, too_many_peaks_threshold=7000,
                 default_precursor_ion_selection_window=1.5,
                 scan_packer=None):
        if log_handler is None:
            log_handler = show_message

        if ms1_peak_picking_args is None:
            ms1_peak_picking_args = {
                "transforms": [],
                "start_mz": 250
            }
        if msn_peak_picking_args is None:
            msn_peak_picking_args = {
                "transforms": []
            }
        if ms1_deconvolution_args is None:
            ms1_deconvolution_args = {
                "scorer": ms_deisotope.scoring.PenalizedMSDeconVFitter(20., 2),
                "charge_range": (1, 8),
                "averagine": ms_deisotope.peptide,
                "use_quick_charge": True,
            }
        if msn_deconvolution_args is None:
            msn_deconvolution_args = {
                "scorer": ms_deisotope.scoring.MSDeconVFitter(10.),
                "charge_range": (1, 8),
                "averagine": ms_deisotope.peptide,
                "use_quick_charge": True,
            }

        Process.__init__(self)
        self.daemon = True
        self.verbose = verbose
        self._init_batch_store()
        self.ms_file_path = ms_file_path
        self.input_queue = input_queue
        self.output_queue = output_queue

        self.ms1_peak_picking_args = ms1_peak_picking_args
        self.msn_peak_picking_args = msn_peak_picking_args
        self.ms1_deconvolution_args = ms1_deconvolution_args
        self.msn_deconvolution_args = msn_deconvolution_args
        self.envelope_selector = envelope_selector
        self.ms1_averaging = ms1_averaging
        self.deconvolute = deconvolute
        self.scan_packer = scan_packer

        self.transformer = None

        self.no_more_event = no_more_event
        self._work_complete = multiprocessing.Event()
        self._error_occurred = multiprocessing.Event()
        self.log_handler = log_handler
        self.too_many_peaks_threshold = too_many_peaks_threshold
        self.default_precursor_ion_selection_window = default_precursor_ion_selection_window

    def make_scan_transformer(self, loader: Union[ScanIterator, RandomAccessScanSource] = None) -> ScanProcessor:
        self.transformer = ScanProcessor(
            loader,
            ms1_peak_picking_args=self.ms1_peak_picking_args,
            msn_peak_picking_args=self.msn_peak_picking_args,
            ms1_deconvolution_args=self.ms1_deconvolution_args,
            msn_deconvolution_args=self.msn_deconvolution_args,
            loader_type=lambda x: x,
            envelope_selector=self.envelope_selector,
            ms1_averaging=self.ms1_averaging,
            default_precursor_ion_selection_window=self.default_precursor_ion_selection_window)
        return self.transformer

    def _process_ms1(self, scan, product_scans) -> Tuple[Scan, List, List[Scan]]:
        scan, priorities, product_scans = self.transformer.process_scan_group(
            scan, product_scans)
        return scan, priorities, product_scans

    def _deconvolute_ms1(self, scan: Scan, priorities: List, product_scans: List[Scan]):
        self.transformer.deconvolute_precursor_scan(scan, priorities, product_scans)

    def _handle_ms1_scan(self, scan: Scan, product_scans: List[Scan], scan_id: str, product_scan_ids: List[str]) -> Tuple[Scan, List[Scan]]:
        try:
            # Check if the m/z array is empty, if so skip the scan. This may trigger
            # arbitrary data loading and decompression behavior, so it can fail.
            if len(scan.arrays[0]) == 0:
                self.skip_scan(scan)
            else:
                try:
                    scan, priorities, product_scans = self._process_ms1(
                        scan, product_scans)
                    if scan is None:
                        # no way to report skip
                        pass
                    else:
                        if self.verbose:
                            self.log_message("Handling Precursor Scan %r with %d peaks" % (scan.id, len(scan.peak_set)))
                        if self.deconvolute:
                            self._deconvolute_ms1(scan, priorities, product_scans)
                        self.send_scan(scan)
                except (KeyboardInterrupt, SystemExit) as e:
                    raise
                except NoIsotopicClustersError as e:
                    self.log_message("No isotopic clusters were extracted from scan %s (%r peaks)" % (
                        e.scan_id, len(scan.peak_set)))
                    self.skip_scan(scan)
                except EmptyScanError as e:
                    self.skip_scan(scan)
                except Exception as e:
                    self.skip_scan(scan)
                    self.log_error(e, scan_id, scan, (product_scan_ids))
        except (KeyboardInterrupt, SystemExit) as e:
            raise
        except Exception as err:
            self.skip_scan(scan)
            self.log_error(err, scan_id, scan, product_scan_ids)
        return scan, product_scans

    def _process_msn(self, product_scan: Scan):
        self.transformer.pick_product_scan_peaks(product_scan)

    def _deconvolute_msn(self, product_scan: Scan):
        self.transformer.deconvolute_product_scan(product_scan)

    def _handle_msn(self, product_scan: Scan, precursor_scan: Scan):
        try:
            self._process_msn(product_scan)
            if self.verbose:
                self.log_message("Handling Product Scan %r with %d peaks (%0.3f/%0.3f, %r)" % (
                    product_scan.id, len(product_scan.peak_set), product_scan.precursor_information.mz,
                    product_scan.precursor_information.extracted_mz,
                    product_scan.precursor_information.defaulted))
            if self.deconvolute:
                self._deconvolute_msn(product_scan)
                if precursor_scan is None and product_scan.precursor_information:
                    product_scan.precursor_information.default(orphan=True)
            self.send_scan(product_scan)
        except (KeyboardInterrupt, SystemExit) as e:
            raise
        except NoIsotopicClustersError as e:
            self.log_message("No isotopic clusters were extracted from scan %s (%r)" % (
                e.scan_id, len(product_scan.peak_set)))
            self.skip_scan(product_scan)
        except EmptyScanError as e:
            self.skip_scan(product_scan)
        except Exception as e:
            self.skip_scan(product_scan)
            self.log_error(e, product_scan.id, product_scan, ())

    def handle_scan_bunch(self, scan: Scan, product_scans: List[Scan], scan_id: str, product_scan_ids: List[str], process_msn: bool=True):
        transformer = self.transformer
        # handle the MS1 scan if it is present
        if scan is not None:
            scan, product_scans = self._handle_ms1_scan(
                scan, product_scans, scan_id, product_scan_ids)
        for product_scan in product_scans:
            # no way to report skip
            try:
                if product_scan is None:
                    continue
                if len(product_scan.arrays[0]) == 0 or (not process_msn):
                    self.skip_scan(product_scan)
                    continue
                self._handle_msn(product_scan, scan)
            except (KeyboardInterrupt, SystemExit) as e:
                raise
            except Exception as err:
                self.skip_scan(product_scan)
                self.log_error(err, product_scan.id,
                               product_scan, [])

    def _open_ms_file(self) -> Union[RandomAccessScanSource, ScanIterator]:
        self.loader = MSFileLoader(self.ms_file_path, decode_binary=False)
        return self.loader

    def _make_batch_loader(self, loader: Union[ScanIterator, RandomAccessScanSource]) -> ScanBunchLoader:
        return ScanBunchLoader(loader)

    def run(self):
        signal.signal(signal.SIGINT, signal.Handlers.SIG_IGN)
        self.try_set_process_name("ms-deisotope-deconv")
        self._silence_loggers()
        loader = self._open_ms_file()
        queued_loader = self._make_batch_loader(loader)
        transformer = self.make_scan_transformer(loader)

        has_input: bool = True
        i: int = 0
        last: int = 0
        scan = None
        product_scans = []
        while has_input:
            try:
                scan_id, product_scan_ids, process_msn = self.get_work(True, 10)
            except QueueEmpty:
                if self.no_more_event is not None and self.no_more_event.is_set():
                    has_input = False
                continue

            i += 1 + len(product_scan_ids)
            if scan_id == DONE:
                has_input = False
                break

            try:
                queued_loader.put(scan_id, product_scan_ids)
                scan, product_scans = queued_loader.get()
            except (KeyboardInterrupt, SystemExit) as e:
                self.log_message("Interrupt received.")
                self._error_occurred.set()
                break
            except Exception as e:
                self.log_message("Something went wrong when loading bunch (%s): %r.\nRecovery is not possible." % (
                    (scan_id, product_scan_ids), e))
                self._error_occurred.set()
                break

            self.handle_scan_bunch(
                scan, product_scans,
                scan_id, product_scan_ids,
                process_msn)

            if (i - last) > 1000:
                last = i
                self.output_queue.join()

        self.log_message("Done (%d scans)" % i)

        if self.no_more_event is None:
            self.output_queue.put((DONE, DONE, DONE))

        if os.getenv("MS_DEISOTOPE_DEBUG"):
            self._dump_averagine_caches(transformer)

        self._work_complete.set()

    def _dump_averagine_caches(self, transformer: ScanProcessor):
        pid = os.getpid()
        fname = "ms_deisotope_averagine_cache_ms$_%s.pkl" % pid
        with open(fname.replace("$", '1'), 'wb') as fh:
            pickle.dump(transformer.ms1_deconvolution_args.get('averagine'), fh, 2)
        with open(fname.replace("$", 'n'), 'wb') as fh:
            pickle.dump(transformer.msn_deconvolution_args.get('averagine'), fh, 2)
