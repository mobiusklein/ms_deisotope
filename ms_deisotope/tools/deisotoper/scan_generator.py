'''Defines the base class for organizing the multiprocessing deconvolution algorithm
encapsulating all behaviors from start to finish.
'''
import multiprocessing

from multiprocessing import JoinableQueue

from ms_deisotope.processor import MSFileLoader

from ms_deisotope.feature_map.quick_index import index as build_scan_index
from ms_deisotope.task import TaskBase

from .collator import ScanCollator
from .process import ScanIDYieldingProcess, DeconvolutingScanTransformingProcess


class ScanGeneratorBase(object):
    """The partial base class for the type that generates scan processing tasks.

    Defines getter and setters for properties expected to be used.

    Attributes
    ----------
    scan_source: :class:`~ms_deisotope.data_source.RandomAccessScanSource`
        Where to read scan data from.
    deconvoluting: bool
        Whether the transformation process will deisotope and charge state deconvolute.
    ms1_averaging: int
        The number of scans on either side of each MS1 scan to average with.
    extract_only_tandem_envelopes: bool
        Whether or not to process intervals around MSn precursors in MS1 scans or not.
    ignore_tandem_scans: bool
        Whether or not to ignore MSn scans.
    """
    def configure_iteration(self, start_scan=None, end_scan=None, max_scans=None):
        """Set :attr:`_iterator` to the result of :meth:`make_iterator`

        Parameters
        ----------
        start_scan : int, optional
            The starting scan index
        end_scan : int, optional
            The ending scan index
        max_scans : int, optional
            The maximum number of scans to process.
        """
        raise NotImplementedError()

    def make_iterator(self, start_scan=None, end_scan=None, max_scans=None):
        """Create an :class:`Iterator` object over the scan product stream

        Parameters
        ----------
        start_scan : int, optional
            The starting scan index
        end_scan : int, optional
            The ending scan index
        max_scans : int, optional
            The maximum number of scans to process.

        Yields
        ------
        :class:`~.ProcessedScan`
        """
        raise NotImplementedError()

    _iterator = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._iterator is None:
            self._iterator = self.make_iterator()
        return next(self._iterator)

    def next(self):
        return self.__next__()

    def close(self):
        pass

    @property
    def scan_source(self):
        return None

    _deconvoluting = False

    @property
    def deconvoluting(self):
        return self._deconvoluting

    @deconvoluting.setter
    def deconvoluting(self, value):
        self._deconvoluting = value

    _ms1_averaging = 0

    @property
    def ms1_averaging(self):
        return self._ms1_averaging

    @ms1_averaging.setter
    def ms1_averaging(self, value):
        self._ms1_averaging = value

    _ignore_tandem_scans = False

    @property
    def ignore_tandem_scans(self):
        return self._ignore_tandem_scans

    @ignore_tandem_scans.setter
    def ignore_tandem_scans(self, value):
        self._ignore_tandem_scans = value

    _extract_only_tandem_envelopes = False

    @property
    def extract_only_tandem_envelopes(self):
        return self._extract_only_tandem_envelopes

    @extract_only_tandem_envelopes.setter
    def extract_only_tandem_envelopes(self, value):
        self._extract_only_tandem_envelopes = value


class ScanGenerator(TaskBase, ScanGeneratorBase):
    def __init__(self, ms_file, number_of_helpers=4,
                 ms1_peak_picking_args=None, msn_peak_picking_args=None,
                 ms1_deconvolution_args=None, msn_deconvolution_args=None,
                 extract_only_tandem_envelopes=False, ignore_tandem_scans=False,
                 ms1_averaging=0, default_precursor_ion_selection_window=1.5,
                 deconvolute=True, verbose=False):
        self.ms_file = ms_file
        self.ignore_tandem_scans = ignore_tandem_scans

        self.scan_ids_exhausted_event = multiprocessing.Event()

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
        self.default_precursor_ion_selection_window = default_precursor_ion_selection_window

        self.deconvoluting = deconvolute
        self.ms1_deconvolution_args = ms1_deconvolution_args
        self.msn_deconvolution_args = msn_deconvolution_args
        self.extract_only_tandem_envelopes = extract_only_tandem_envelopes
        self._scan_interval_tree = None
        self.verbose = verbose
        self.log_controller = self.ipc_logger()

    @property
    def scan_source(self):
        return self.ms_file

    def join(self):
        if self._scan_yielder_process is not None:
            self._scan_yielder_process.join()
        if self._deconv_process is not None:
            self._deconv_process.join()
        if self._deconv_helpers is not None:
            for helper in self._deconv_helpers:
                helper.join()

    def _terminate(self):
        if self._scan_yielder_process is not None:
            self._scan_yielder_process.terminate()
        if self._deconv_process is not None:
            self._deconv_process.terminate()
        if self._deconv_helpers is not None:
            for helper in self._deconv_helpers:
                helper.terminate()

    def _preindex_file(self):
        reader = MSFileLoader(self.ms_file, use_index=False)
        try:
            reader.prebuild_byte_offset_file(self.ms_file)
        except AttributeError:
            # the type does not support this type of indexing
            pass
        except IOError:
            # the file could not be written
            pass
        except Exception as e:
            # something else went wrong
            self.error("An error occurred while pre-indexing.", e)

    def _make_interval_tree(self, start_scan, end_scan):
        reader = MSFileLoader(self.ms_file, decode_binary=False)
        if start_scan is not None:
            start_ix = reader.get_scan_by_id(start_scan).index
        else:
            start_ix = 0
        if end_scan is not None:
            end_ix = reader.get_scan_by_id(end_scan).index
        else:
            end_ix = len(reader)
        reader.reset()
        _index, interval_tree = build_scan_index(
            reader, self.number_of_helpers + 1, (start_ix, end_ix))
        self._scan_interval_tree = interval_tree

    def _make_transforming_process(self):
        return DeconvolutingScanTransformingProcess(
            self.ms_file,
            self._input_queue,
            self._output_queue,
            self.scan_ids_exhausted_event,
            ms1_peak_picking_args=self.ms1_peak_picking_args,
            msn_peak_picking_args=self.msn_peak_picking_args,
            ms1_deconvolution_args=self.ms1_deconvolution_args,
            msn_deconvolution_args=self.msn_deconvolution_args,
            envelope_selector=self._scan_interval_tree,
            log_handler=self.log_controller.sender(),
            ms1_averaging=self.ms1_averaging,
            default_precursor_ion_selection_window=self.default_precursor_ion_selection_window,
            deconvolute=self.deconvoluting,
            verbose=self.verbose)

    def _make_collator(self):
        return ScanCollator(
            self._output_queue, self.scan_ids_exhausted_event, self._deconv_helpers,
            self._deconv_process, input_queue=self._input_queue,
            include_fitted=not self.deconvoluting)

    def _initialize_workers(self, start_scan=None, end_scan=None, max_scans=None):
        try:
            self._input_queue = JoinableQueue(int(1e6))
            self._output_queue = JoinableQueue(int(1e6))
        except OSError:
            # Not all platforms permit limiting the size of queues
            self._input_queue = JoinableQueue()
            self._output_queue = JoinableQueue()

        self._preindex_file()

        if self.extract_only_tandem_envelopes:
            self.log("Constructing Scan Interval Tree")
            self._make_interval_tree(start_scan, end_scan)

        self._terminate()
        self._scan_yielder_process = ScanIDYieldingProcess(
            self.ms_file, self._input_queue, start_scan=start_scan, end_scan=end_scan,
            max_scans=max_scans, no_more_event=self.scan_ids_exhausted_event,
            ignore_tandem_scans=self.ignore_tandem_scans, batch_size=1)
        self._scan_yielder_process.start()

        self._deconv_process = self._make_transforming_process()

        self._deconv_helpers = []

        for _ in range(self.number_of_helpers):
            self._deconv_helpers.append(self._make_transforming_process())
        self._deconv_process.start()

        self._order_manager = self._make_collator()

    def make_iterator(self, start_scan=None, end_scan=None, max_scans=None):
        self._initialize_workers(start_scan, end_scan, max_scans)

        for scan in self._order_manager:
            yield scan
        self.log_controller.stop()
        self.join()
        self._terminate()

    def configure_iteration(self, start_scan=None, end_scan=None, max_scans=None):
        self._iterator = self.make_iterator(start_scan, end_scan, max_scans)

    def close(self):
        self._terminate()
