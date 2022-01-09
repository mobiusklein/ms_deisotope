from __future__ import print_function
import logging
import multiprocessing
import threading
import traceback
import six

from datetime import datetime

try:
    from queue import Empty
except ImportError:
    from Queue import Empty


logger = logging.getLogger("ms_deisotope.task")


def ensure_text(obj):
    if six.PY2:
        return unicode(obj)
    else:
        return str(obj)


def fmt_msg(*message):
    return u"%s %s" % (ensure_text(datetime.now().isoformat(' ')), u', '.join(map(ensure_text, message)))


def printer(obj, *message):
    print(fmt_msg(*message))


def show_message(*message):
    print(fmt_msg(*message))


def debug_printer(obj, *message):
    if obj.in_debug_mode():
        print(u"DEBUG:" + fmt_msg(*message))


class CallInterval(object):
    """Call a function every `interval` seconds from
    a separate thread.

    Attributes
    ----------
    stopped: threading.Event
        A semaphore lock that controls when to run `call_target`
    call_target: callable
        The thing to call every `interval` seconds
    args: iterable
        Arguments for `call_target`
    interval: number
        Time between calls to `call_target`
    """

    def __init__(self, interval, call_target, *args):
        self.stopped = threading.Event()
        self.interval = interval
        self.call_target = call_target
        self.args = args
        self.thread = threading.Thread(target=self.mainloop)
        self.thread.daemon = True

    def mainloop(self):
        while not self.stopped.wait(self.interval):
            try:
                self.call_target(*self.args)
            except (KeyboardInterrupt, SystemExit):
                self.stop()
            except Exception as e:
                logger.exception("An error occurred in %r", self, exc_info=e)

    def start(self):
        self.thread.start()

    def stop(self):
        self.stopped.set()


class MessageSpooler(object):
    """An IPC-based logging helper

    Attributes
    ----------
    halting : bool
        Whether the object is attempting to
        stop, so that the internal thread can
        tell when it should stop and tell other
        objects using it it is trying to stop
    handler : Callable
        A Callable object which can be used to do
        the actual logging
    message_queue : multiprocessing.Queue
        The Inter-Process Communication queue
    thread : threading.Thread
        The internal listener thread that will consume
        message_queue work items
    """
    def __init__(self, handler):
        self.handler = handler
        self.message_queue = multiprocessing.Queue()
        self.halting = False
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True
        self.thread.start()

    def run(self):
        while not self.halting:
            try:
                message = self.message_queue.get(True, 2)
                self.handler(*message)
            except Empty:
                continue

    def stop(self):
        self.halting = True
        self.thread.join()

    def sender(self):
        return MessageSender(self.message_queue)


class MessageSender(object):
    """A simple callable for pushing objects into an IPC
    queue.

    Attributes
    ----------
    queue : multiprocessing.Queue
        The Inter-Process Communication queue
    """
    def __init__(self, queue):
        self.queue = queue

    def __call__(self, *message):
        self.send(*message)

    def send(self, *message):
        self.queue.put(message)


class LogUtilsMixin(object):

    logger_state = None
    print_fn = printer
    debug_print_fn = debug_printer
    error_print_fn = printer
    warn_print_fn = printer

    @classmethod
    def log_with_logger(cls, logger):
        cls.logger_state = logger
        cls.print_fn = logger.info
        cls.debug_print_fn = logger.debug
        cls.error_print_fn = logger.error
        cls.warn_print_fn = logger.warn

    def pipe_from(self, task):
        task.print_fn = self.log
        task.debug_print_fn = self.debug
        task.error_print_fn = self.error
        task.warn_print_fn = self.warn

    def in_debug_mode(self):
        if self._debug_enabled is None:
            logger_state = self.logger_state
            if logger_state is not None:
                self._debug_enabled = logger_state.isEnabledFor("DEBUG")
        return bool(self._debug_enabled)

    def log(self, *message):
        self.print_fn(u', '.join(map(ensure_text, message)))

    def warn(self, *message):
        self.warn_print_fn(u', '.join(map(ensure_text, message)))

    def debug(self, *message):
        self.debug_print_fn(u', '.join(map(ensure_text, message)))

    def error(self, *message, **kwargs):
        exception = kwargs.get("exception")
        self.error_print_fn(u', '.join(map(ensure_text, message)))
        if exception is not None:
            self.error_print_fn(traceback.format_exc(exception))

    def ipc_logger(self, handler=None):
        if handler is None:
            def default_handler(message):
                self.log(message)
            handler = default_handler
        return MessageSpooler(handler)


class ProgressUpdater(LogUtilsMixin):
    def __init__(self, *args, **kwargs):
        kwargs['item_show_func'] = self._prepare_message
        try:
            import click
            self.progress_bar = click.progressbar(*args, **kwargs)
        except ImportError:
            self.progress_bar = None

    def _prepare_message(self, *message):
        if len(message) == 1:
            if message[0] is None:
                return ""
        return ', '.join(map(str, message))

    def log(self, *message):
        if self.progress_bar is not None:
            if not self.progress_bar.is_hidden:
                self.progress_bar.current_item = self._prepare_message(
                    *message)
                return
        super(ProgressUpdater, self).log(*message)

    def warn(self, *message):
        if self.progress_bar is not None:
            if not self.progress_bar.is_hidden:
                self.progress_bar.current_item = "WARN: " + self._prepare_message(
                    *message)
                return
        super(ProgressUpdater, self).warn(*message)

    def debug(self, *message):
        if self.progress_bar is not None:
            if not self.progress_bar.is_hidden:
                self.progress_bar.current_item = "DEBUG: " + self._prepare_message(
                    *message)
                return
        super(ProgressUpdater, self).debug(*message)

    def begin(self):
        if self.progress_bar is not None:
            self.progress_bar.__enter__()

    def end(self, exc_type=None, exc_value=None, tb=None):
        if self.progress_bar is not None:
            self.progress_bar.__exit__(exc_type, exc_value, tb)

    def __enter__(self):
        return self.begin()

    def __exit__(self, exc_type, exc_value, tb):
        return self.end(exc_type, exc_value, tb)

    def update(self, i, *message):
        if self.progress_bar is not None:
            self.progress_bar.update(i)
        if message:
            self.log(*message)
