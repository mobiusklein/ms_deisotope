"""Manages keeping scans delivered out-of-order in-order for
writing to disk.
"""

try:
    from Queue import Empty as QueueEmpty
except ImportError:
    from queue import Empty as QueueEmpty

from ms_deisotope.data_source.common import ProcessedScan
from ms_deisotope.task import TaskBase, CallInterval

from .process import (
    SCAN_STATUS_SKIP, DONE)


class ScanCollator(TaskBase):
    """Collates incoming scan bunches from multiple
    ScanTransformingProcesses, passing them along in
    the correct order.

    Attributes
    ----------
    count_jobs_done : int
        The number of scan bunches taken from :attr:`queue`
    count_since_last : int
        The number of work-cycles since the last scan bunch
        has been yielded
    done_event : multiprocessing.Event
        An IPC Event to indicate that all scan ids have been
        sent to the worker processes
    helper_producers : list
        A list of ScanTransformingProcesses
    include_fitted : bool
        Whether or not to save the raw fitted peaks for each
        scan produced. When this is `False`, they will be
        discarded and memory will be saved
    last_index : int
        The index of the last scan yielded through the iterator
        loop. This controls the next scan to be yielded and any
        waiting conditions
    primary_worker : ScanTransformingProcess
        The first worker to start consuming scans which will dictate
        the first handled index. Is required to run in isolation
        from other worker processes to insure that the first scan
        arrives in order
    queue : multiprocessing.Queue
        The IPC queue that all workers place their results on
        to be consumed and yielded in order
    started_helpers : bool
        Whether or not the additional workers in :attr:`helper_producers`
        have been started or not
    waiting : dict
        A mapping from scan index to `Scan` object. Used to serve
        scans through the iterator when their index is called for
    """
    _log_received_scans = False

    def __init__(self, queue, done_event, helper_producers=None, primary_worker=None,
                 include_fitted=False, input_queue=None):
        if helper_producers is None:
            helper_producers = []
        self.queue = queue
        self.last_index = None
        self.count_jobs_done = 0
        self.count_since_last = 0
        self.waiting = {}
        self.done_event = done_event
        self.helper_producers = helper_producers
        self.started_helpers = False
        self.primary_worker = primary_worker
        self.include_fitted = include_fitted
        self.input_queue = input_queue

    def all_workers_done(self):
        '''
        Check if all of the worker processes have set their "work done"
        flag.

        Returns
        -------
        bool
        '''
        if self.done_event.is_set():
            if self.primary_worker.all_work_done():
                for helper in self.helper_producers:
                    if not helper.all_work_done():
                        return False
                return True
            else:
                return False
        return False

    def store_item(self, item, index):
        """Stores an incoming work-item for easy
        access by its `index` value. If configuration
        requires it, this will also reduce the number
        of peaks in `item`.

        Parameters
        ----------
        item : str or ProcessedScan
            Either a stub indicating why this work item
            is not
        index : int
            Scan index to store
        """
        if self._log_received_scans:
            self.log("-- received %d: %s" % (index, item))
        self.waiting[index] = item
        if not self.include_fitted and isinstance(item, ProcessedScan):
            item.peak_set = []

    def consume(self, timeout=10):
        """Fetches the next work item from the input
        queue :attr:`queue`, blocking for at most `timeout` seconds.

        Parameters
        ----------
        timeout : int, optional
            The duration to allow the process to block
            for while awaiting new work items.

        Returns
        -------
        bool
            Whether or not a new work item was found waiting
            on the :attr:`queue`
        """
        blocking = timeout != 0
        try:
            item, index, _ = self.queue.get(blocking, timeout)
            self.queue.task_done()
            # DONE message may be sent many times.
            while item == DONE:
                item, index, _ = self.queue.get(blocking, timeout)
                self.queue.task_done()
            self.store_item(item, index)
            return True
        except QueueEmpty:
            return False

    def start_helper_producers(self):
        """Starts the additional :class:`ScanTransformingProcess` workers
        in :attr:`helper_producers` if they have not been started already.

        Should only be invoked once
        """
        if self.started_helpers:
            return
        self.started_helpers = True
        for helper in self.helper_producers:
            if helper.is_alive():
                continue
            helper.start()

    def produce(self, scan):
        """Performs any final quality controls on the outgoing
        :class:`ProcessedScan` object and takes care of any internal
        details.

        Resets :attr:`count_since_last` to `0`.

        Parameters
        ----------
        scan : ProcessedScan
            The scan object being finalized for hand-off
            to client code

        Returns
        -------
        ProcessedScan
            The version of `scan` ready to be used by other
            parts of the program
        """
        self.count_since_last = 0
        return scan

    def count_pending_items(self):
        """Count the number of scans that are waiting to be added to the
        write queue.

        Returns
        -------
        int
        """
        return len(self.waiting)

    def drain_queue(self):
        """Try to read a lot of scans from the incoming result queue.

        If there are items pending to be sent to the write queue immediately,
        don't read as many.

        Returns
        -------
        int:
            The number of items drained.
        """
        i = 0
        has_next = self.last_index + 1 not in self.waiting
        while (self.count_pending_items() < (1000 if has_next else 10)
               and self.consume(.1)):
            self.count_jobs_done += 1
            has_next = self.last_index + 1 not in self.waiting
            i += 1
        if i > 15:
            self.log("Drained Output Queue of %d Items" % (i, ))
        return i

    def print_state(self):
        """Log the state of the collator, reporting pending items
        in all output stages.

        If a worker process is done, try to join it.
        """
        try:
            if self.queue.qsize() > 0:
                self.log("%d since last work item" % (self.count_since_last,))
                keys = sorted(self.waiting.keys())
                if len(keys) > 5:
                    self.log("Waiting Keys: %r..." % (keys[:5],))
                else:
                    self.log("Waiting Keys: %r" % (keys,))
                self.log("%d Keys Total" % (len(self.waiting),))
                self.log("The last index handled: %r" % (self.last_index,))
                self.log("Number of items waiting in the queue: %d" %
                         (self.queue.qsize(),))
        except NotImplementedError:
            # Some platforms do not support qsize
            pass
        for worker in ([self.primary_worker] + list(self.helper_producers)):
            code = worker.exitcode
            if code is not None and code != 0:
                self.log("%r has exit code %r" % (worker, code))
                worker.join(5)

    def __iter__(self):
        has_more = True
        # Log the state of the collator every 3 minutes
        status_monitor = CallInterval(60 * 3, self.print_state)
        status_monitor.start()
        while has_more:
            if self.consume(1):
                self.count_jobs_done += 1
                try:
                    if self.queue.qsize() > 500:
                        self.drain_queue()
                except NotImplementedError:
                    # Some platforms do not support qsize
                    self.drain_queue()
            if self.last_index is None:
                keys = sorted(self.waiting)
                if keys:
                    i = 0
                    n = len(keys)
                    found_content = False
                    while i < n:
                        scan = self.waiting.pop(keys[i])
                        if scan == SCAN_STATUS_SKIP:
                            self.last_index = keys[i]
                            i += 1
                            continue
                        else:
                            found_content = True
                            break
                    if found_content:
                        self.last_index = scan.index
                        yield self.produce(scan)
                    if self.last_index is not None:
                        self.start_helper_producers()
            elif self.last_index + 1 in self.waiting:
                while self.last_index + 1 in self.waiting:
                    scan = self.waiting.pop(self.last_index + 1)
                    if scan == SCAN_STATUS_SKIP:
                        self.last_index += 1
                        continue
                    else:
                        self.last_index = scan.index
                        yield self.produce(scan)
            elif len(self.waiting) == 0:
                if self.all_workers_done():
                    self.log("All Workers Claim Done.")
                    has_something = self.consume()
                    self.log("Checked Queue For Work: %r" % has_something)
                    if not has_something and len(self.waiting) == 0 and self.queue.empty():
                        has_more = False
            else:
                self.count_since_last += 1
                if self.count_since_last % 1000 == 0:
                    self.print_state()
        status_monitor.stop()
