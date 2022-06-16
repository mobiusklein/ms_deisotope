from __future__ import print_function

import logging
import logging.handlers
import multiprocessing
import sys
import threading
import traceback
import re

from typing import Dict

from datetime import datetime

from queue import Empty


logger = logging.getLogger("ms_deisotope.task")


def ensure_text(obj):
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

    _debug_enabled = None

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
        kwargs.setdefault("color", True)
        kwargs.setdefault("fill_char", click.style('-', 'blue'))
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


class ProcessAwareFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        d = record.__dict__
        try:
            if d['processName'] == "MainProcess":
                d['maybeproc'] = ''
            else:
                d['maybeproc'] = ":%s:" % d['processName'].replace(
                    "Process", '')
        except KeyError:
            d['maybeproc'] = ''
        return super(ProcessAwareFormatter, self).format(record)


class LevelAwareColoredLogFormatter(ProcessAwareFormatter):
    try:
        from colorama import Fore, Style
        # GREY = Fore.WHITE
        GREY = ''
        BLUE = Fore.BLUE
        GREEN = Fore.GREEN
        YELLOW = Fore.YELLOW
        RED = Fore.RED
        BRIGHT = Style.BRIGHT
        DIM = Style.DIM
        BOLD_RED = Fore.RED + Style.BRIGHT
        RESET = Style.RESET_ALL
    except ImportError:
        GREY = ''
        BLUE = ''
        GREEN = ''
        YELLOW = ''
        RED = ''
        BRIGHT = ''
        DIM = ''
        BOLD_RED = ''
        RESET = ''

    def _colorize_field(self, fmt: str, field: str, color: str) -> str:
        return re.sub("(" + field + ")", color + r"\1" + self.RESET, fmt)


    def _patch_fmt(self, fmt: str, level_color: str) -> str:
        fmt = self._colorize_field(fmt, r"%\(asctime\)s", self.GREEN)
        fmt = self._colorize_field(fmt, r"%\(name\).*?s", self.BLUE)
        fmt = self._colorize_field(fmt, r"%\(message\).*?s", self.GREY)
        if level_color:
            fmt = self._colorize_field(fmt, r"%\(levelname\).*?s", level_color)
        return fmt

    def __init__(self, fmt, level_color=None, **kwargs):
        fmt = self._patch_fmt(fmt, level_color=level_color)
        super().__init__(fmt, **kwargs)


class ColoringFormatter(logging.Formatter):
    level_to_color = {
        logging.INFO: LevelAwareColoredLogFormatter.GREEN,
        logging.DEBUG: LevelAwareColoredLogFormatter.GREY + LevelAwareColoredLogFormatter.DIM,
        logging.WARN: LevelAwareColoredLogFormatter.YELLOW + LevelAwareColoredLogFormatter.BRIGHT,
        logging.ERROR: LevelAwareColoredLogFormatter.BOLD_RED,
        logging.CRITICAL: LevelAwareColoredLogFormatter.BOLD_RED,
        logging.FATAL: LevelAwareColoredLogFormatter.RED + LevelAwareColoredLogFormatter.DIM,
    }

    _formatters: Dict[int, LevelAwareColoredLogFormatter]

    def __init__(self, fmt: str, **kwargs):
        self._formatters = {}
        for level, style in self.level_to_color.items():
            self._formatters[level] = LevelAwareColoredLogFormatter(fmt, level_color=style, **kwargs)

    def format(self, record: logging.LogRecord) -> str:
        fmtr = self._formatters[record.levelno]
        return fmtr.format(record)



def init_logging(filename=None, queue=None):
    logging.basicConfig(
        level="INFO", format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt="%H:%M:%S", handlers=[]
    )
    logging.captureWarnings(True)

    logger = logging.getLogger('ms_deisotope')
    format_string = '[%(asctime)s] %(levelname).1s | %(name)s | %(message)s'
    formatter = ProcessAwareFormatter(format_string, datefmt="%H:%M:%S")
    colorized_formatter = ColoringFormatter(format_string, datefmt="%H:%M:%S")

    # If there was a queue, don't add any other handlers, route all logging through
    # the queue
    if queue:
        queue_handler = logging.handlers.QueueHandler(queue)
        queue_handler.setFormatter(formatter)
        logger.addHandler(queue_handler)
    else:
        # Otherwise, configure handlers for `ms_deisotope`
        stderr_handler = logging.StreamHandler(sys.stderr)
        if sys.stderr.isatty():
            stderr_handler.setFormatter(colorized_formatter)
        else:
            stderr_handler.setFormatter(formatter)
        logger.addHandler(stderr_handler)
        if filename:
            file_handler = logging.FileHandler(filename=filename, mode='w', encoding='utf8')
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
    LogUtilsMixin.log_with_logger(logger)
