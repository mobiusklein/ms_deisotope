import os
import logging
import sys
import multiprocessing
import traceback

from collections import deque
from multiprocessing import Process
try:
    from Queue import Empty as QueueEmpty
except ImportError:
    from queue import Empty as QueueEmpty
try:
    import cPickle as pickle
except ImportError:
    import pickle


import ms_deisotope

from ms_deisotope.processor import (
    ScanProcessor, MSFileLoader,
    NoIsotopicClustersError, EmptyScanError)

from ms_deisotope.task import show_message


DONE = b"--NO-MORE--"
SCAN_STATUS_GOOD = b"good"
SCAN_STATUS_SKIP = b"skip"


class ScanIDYieldingProcess(Process):

    def __init__(self, ms_file_path, queue, start_scan=None, max_scans=None, end_scan=None,
                 no_more_event=None, ignore_tandem_scans=False, batch_size=1, log_handler=None):
        if log_handler is None:
            log_handler = show_message
        Process.__init__(self)
        self.daemon = True
        self.ms_file_path = ms_file_path
        self.queue = queue
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

    def _make_scan_batch(self):
        batch = []
        scan_ids = []
        for _ in range(self.batch_size):
            try:
                bunch = next(self.loader)
                scan, products = bunch
                products = [prod for prod in products if prod.index <= self.end_scan_index]
                if scan is not None:
                    scan_id = scan.id
                else:
                    scan_id = None
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

    def run(self):
        self.loader = MSFileLoader(self.ms_file_path, decode_binary=False)

        if self.start_scan is not None:
            try:
                self.loader.start_from_scan(
                    self.start_scan, require_ms1=self.loader.has_ms1_scans(), grouped=True)
            except IndexError as e:
                self.log_handler("An error occurred while locating start scan", e)
                self.loader.reset()
                self.loader.make_iterator(grouped=True)
            except AttributeError:
                self.log_handler("The reader does not support random access, start time will be ignored", e)
                self.loader.reset()
                self.loader.make_iterator(grouped=True)
        else:
            self.loader.make_iterator(grouped=True)

        count = 0
        last = 0
        if self.max_scans is None:
            max_scans = float('inf')
        else:
            max_scans = self.max_scans

        end_scan = self.end_scan
        if end_scan is None:
            try:
                self.end_scan_index = len(self.loader)
            except AttributeError:
                self.end_scan_index = sys.maxint
        else:
            self.end_scan_index = self.loader.get_scan_by_id(end_scan).index
        while count < max_scans:
            try:
                batch, ids = self._make_scan_batch()
                if len(batch) > 0:
                    self.queue.put(batch)
                count += len(ids)
                if (count - last) > 1000:
                    last = count
                    self.queue.join()
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
            self.queue.put(DONE)


class ScanTransformMixin(object):
    def log_error(self, error, scan_id, scan, product_scan_ids):
        tb = traceback.format_exc()
        self.log_handler(
            "An %r occurred for %s (index %r) in Process %r\n%s" % (
                error, scan_id, scan.index, multiprocessing.current_process(),
                tb))

    def _init_batch_store(self):
        self._batch_store = deque()

    def get_work(self, block=True, timeout=30):
        if self._batch_store:
            return self._batch_store.popleft()
        else:
            batch = self.input_queue.get(block, timeout)
            self._batch_store.extend(batch)
            result = self._batch_store.popleft()
            return result

    def log_message(self, message):
        self.log_handler(message + ", %r" %
                         (multiprocessing.current_process().name))

    def skip_entry(self, index, ms_level):
        self.output_queue.put((SCAN_STATUS_SKIP, index, ms_level))

    def skip_scan(self, scan):
        self.output_queue.put((SCAN_STATUS_SKIP, scan.index, scan.ms_level))

    def send_scan(self, scan):
        scan = scan.pack()
        # this attribute is not needed, and for MS1 scans is dangerous
        # to pickle.
        # It can pull other scans which may not yet have been packed
        # into the message sent back to the main process which in
        # turn can form a reference cycle and eat a lot of memory
        scan.product_scans = []
        self.output_queue.put((scan, scan.index, scan.ms_level))

    def all_work_done(self):
        return self._work_complete.is_set()

    def make_scan_transformer(self, loader=None):
        raise NotImplementedError()


class ScanBunchLoader(object):

    def __init__(self, mzml_loader):
        self.loader = mzml_loader
        self.queue = deque()

    def put(self, scan_id, product_scan_ids):
        self.queue.append((scan_id, product_scan_ids))

    def get(self):
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


class DeconvolutingScanTransformingProcess(Process, ScanTransformMixin):
    """DeconvolutingScanTransformingProcess describes a child process that consumes scan id bunches
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

    def __init__(self, ms_file_path, input_queue, output_queue,
                 no_more_event=None, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 envelope_selector=None, ms1_averaging=0, log_handler=None,
                 deconvolute=True, verbose=False, too_many_peaks_threshold=7000,
                 default_precursor_ion_selection_window=1.5):
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

        self.transformer = None

        self.no_more_event = no_more_event
        self._work_complete = multiprocessing.Event()
        self.log_handler = log_handler
        self.too_many_peaks_threshold = too_many_peaks_threshold
        self.default_precursor_ion_selection_window = default_precursor_ion_selection_window

    def make_scan_transformer(self, loader=None):
        transformer = ScanProcessor(
            loader,
            ms1_peak_picking_args=self.ms1_peak_picking_args,
            msn_peak_picking_args=self.msn_peak_picking_args,
            ms1_deconvolution_args=self.ms1_deconvolution_args,
            msn_deconvolution_args=self.msn_deconvolution_args,
            loader_type=lambda x: x,
            envelope_selector=self.envelope_selector,
            ms1_averaging=self.ms1_averaging,
            default_precursor_ion_selection_window=self.default_precursor_ion_selection_window)
        return transformer

    def handle_scan_bunch(self, scan, product_scans, scan_id, product_scan_ids, process_msn=True):
        transformer = self.transformer
        # handle the MS1 scan if it is present
        if scan is not None:
            try:
                if len(scan.arrays[0]) == 0:
                    self.skip_scan(scan)
                else:
                    try:
                        scan, priorities, product_scans = transformer.process_scan_group(
                            scan, product_scans)
                        if scan is None:
                            # no way to report skip
                            pass
                        else:
                            if self.verbose:
                                self.log_message("Handling Precursor Scan %r with %d peaks" % (scan.id, len(scan.peak_set)))
                            if self.deconvolute:
                                transformer.deconvolute_precursor_scan(
                                    scan, priorities, product_scans)
                            self.send_scan(scan)
                    except (KeyboardInterrupt, SystemExit) as e:
                        raise
                    except NoIsotopicClustersError as e:
                        self.log_message("No isotopic clusters were extracted from scan %s (%r)" % (
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
        for product_scan in product_scans:
            # no way to report skip
            try:
                if product_scan is None:
                    continue
                if len(product_scan.arrays[0]) == 0 or (not process_msn):
                    self.skip_scan(product_scan)
                    continue
                try:
                    transformer.pick_product_scan_peaks(product_scan)
                    if self.verbose:
                        self.log_message("Handling Product Scan %r with %d peaks (%0.3f/%0.3f, %r)" % (
                            product_scan.id, len(product_scan.peak_set), product_scan.precursor_information.mz,
                            product_scan.precursor_information.extracted_mz,
                            product_scan.precursor_information.defaulted))
                    if self.deconvolute:
                        transformer.deconvolute_product_scan(product_scan)
                        if scan is None:
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
                    self.log_error(e, product_scan.id,
                                   product_scan, (product_scan_ids))
            except (KeyboardInterrupt, SystemExit) as e:
                raise
            except Exception as err:
                self.skip_scan(product_scan)
                self.log_error(err, product_scan.id,
                               product_scan, [])

    def _silence_loggers(self):
        nologs = ["deconvolution_scan_processor"]
        if not self.deconvolute:
            nologs.append("deconvolution")

        debug_mode = os.getenv("MS_DEISOTOPE_DEBUG")
        if debug_mode:
            handler = logging.FileHandler("ms-deisotope-deconvolution-debug-%s.log" % (os.getpid()), 'w')
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

    def run(self):
        loader = MSFileLoader(self.ms_file_path, decode_binary=False)
        queued_loader = ScanBunchLoader(loader)

        has_input = True
        transformer = self.make_scan_transformer(loader)
        self.transformer = transformer
        self._silence_loggers()
        i = 0
        last = 0
        while has_input:
            try:
                scan_id, product_scan_ids, process_msn = self.get_work(True, 10)
                self.input_queue.task_done()
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
                break
            except Exception as e:
                self.log_message("Something went wrong when loading bunch (%s): %r.\nRecovery is not possible." % (
                    (scan_id, product_scan_ids), e))

            self.handle_scan_bunch(scan, product_scans, scan_id, product_scan_ids, process_msn)
            if (i - last) > 1000:
                last = i
                self.output_queue.join()

        self.log_message("Done (%d scans)" % i)

        if self.no_more_event is None:
            self.output_queue.put((DONE, DONE, DONE))

        if os.getenv("MS_DEISOTOPE_DEBUG"):
            self._dump_averagine_caches(transformer)

        self._work_complete.set()

    def _dump_averagine_caches(self, transformer):
        pid = os.getpid()
        fname = "ms_deisotope_averagine_cache_ms$_%s.pkl" % pid
        with open(fname.replace("$", '1'), 'wb') as fh:
            pickle.dump(transformer.ms1_deconvolution_args.get('averagine'), fh, 2)
        with open(fname.replace("$", 'n'), 'wb') as fh:
            pickle.dump(transformer.msn_deconvolution_args.get('averagine'), fh, 2)
