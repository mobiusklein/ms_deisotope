"""
Helper functions and types for CLI tools.

Dependent upon :mod:`click`.
"""
import ast
import itertools
import json
import multiprocessing
import os
import platform
import re
import sys
import warnings

from typing import Dict, Iterable, Iterator, Generic, TYPE_CHECKING, Union, TypeVar

import click
from click._termui_impl import ProgressBar

from brainpy import periodic_table

from ms_deisotope.averagine import (
    Averagine, glycan as n_glycan_averagine, permethylated_glycan,
    peptide, glycopeptide, heparin, heparan_sulfate)

from ms_deisotope.task.log_utils import fmt_msg
from ms_deisotope.data_source import ScanBase, ScanBunch


T = TypeVar('T')


def processes_option(f) -> click.Option:
    """
    A re-usable process number option.

    Returns
    -------
    click.Option
    """
    opt = click.option(
        "-p", "--processes", 'processes', type=click.IntRange(1, multiprocessing.cpu_count()),
        default=min(multiprocessing.cpu_count(), 4), help=('Number of worker processes to use. Defaults to 4 '
                                                           'or the number of CPUs, whichever is lower'))
    return opt(f)


def validate_element(element: str) -> bool:
    """Validate that an element string is known"""
    valid = element in periodic_table
    if not valid:
        raise click.Abort("%r is not a valid element" % element)
    return valid


def parse_averagine_formula(formula: str) -> Averagine:
    """Parse an formula into an :class:`Averagine` instance"""
    if isinstance(formula, Averagine):
        return formula
    return Averagine({k: float(v) for k, v in re.findall(r"([A-Z][a-z]*)([0-9\.]*)", formula)
                      if float(v or 0) > 0 and validate_element(k)})


averagines: Dict[str, Averagine] = {
    'glycan': n_glycan_averagine,
    'permethylated-glycan': permethylated_glycan,
    'peptide': peptide,
    'glycopeptide': glycopeptide,
    'heparin': heparin,
    "heparan-sulfate": heparan_sulfate
}


def validate_averagine(averagine_string: Union[str, Averagine]) -> Averagine:
    """
    Validate that the input is an averagine or can be converted into one.

    Returns
    -------
    Averagine
    """
    if isinstance(averagine_string, Averagine):
        return averagine_string
    if averagine_string in averagines:
        return averagines[averagine_string]
    else:
        return parse_averagine_formula(averagine_string)


class AveragineParamType(click.types.StringParamType):
    """A parameter type that converts to an :class:`Averagine` object."""

    name = "MODEL"

    models = averagines

    def convert(self, value, param, ctx):
        return validate_averagine(value)

    def get_metavar(self, param):
        return '[%s]' % '|'.join(sorted(averagines.keys()))

    def get_missing_message(self, param):
        return 'Choose from %s, or provide a formula.' % ', '.join(self.choices)


def register_debug_hook():
    """A hook to automatically drop into an interactive debugger on an error."""
    import traceback

    def info(type, value, tb):
        if hasattr(sys, 'ps1') or not sys.stderr.isatty():
            sys.__excepthook__(type, value, tb)
        else:
            try:
                import ipdb as pdb_api
            except ImportError:
                import pdb as pdb_api
            traceback.print_exception(type, value, tb)
            pdb_api.post_mortem(tb)

    sys.excepthook = info
    # logging.basicConfig(level="DEBUG")


def envar_to_bool(env_val: str) -> bool:
    """Coerce a string easily put into an env var into a boolean"""
    env_val = env_val.lower()
    if not env_val:
        return False
    if env_val in ('0', 'no', 'false', 'off'):
        return False
    elif env_val in ('1', 'yes', 'true', 'on'):
        return True
    return None


def is_debug_mode() -> bool:
    """Detect if the program is being told to launch in debug mode by the environment"""
    env_val = os.environ.get('MS_DEISOTOPE_DEBUG', '').lower()
    val = envar_to_bool(env_val)
    if val is not None:
        return val
    else:
        warnings.warn("MS_DEISOTOPE_DEBUG value %r was not recognized. Enabling debug mode" % (env_val, ))
        return True


class Spinner(object):
    """Super-simple synchronous CLI spinner implementation."""

    _default_symbols = ['-', '/', '|', '\\']

    def __init__(self, stream=None, symbols=None, title=None):
        if symbols is None:
            symbols = list(self._default_symbols)
        if stream is None:
            stream = click.get_text_stream('stdout')
        self.symbol_cycle = itertools.cycle(symbols)
        self.stream = stream
        self.title = title

    def _next_symbol(self):
        symbol = next(self.symbol_cycle)
        return symbol

    def update(self, *args, **kwargs):
        if self.stream.isatty():
            return self._update(*args, **kwargs)

    def _update(self, *args, **kwargs):
        self.stream.write("\b")
        self.stream.flush()
        symbol = self._next_symbol()
        self.stream.write(symbol)
        self.stream.flush()

    def __enter__(self):
        if self.title is not None:
            self.stream.write(self.title)
            self.stream.write("  ")
            self.stream.flush()
        return self

    def __exit__(self, *args):
        self.stream.write("\n")
        self.stream.flush()

    def __repr__(self):
        template = "{self.__class__.__name__}({self.stream}, {self.symbol_cycle})"
        return template.format(self=self)

spinner = Spinner


class ProgressLogger(Iterable[T]):
    """
    A simple text logger that wraps an iterable and logs update messages as chunks are requested.

    This class tries to emulate :func:`click.progressbar` for use when the attached STDERR stream is not
    a terminal.
    """

    iterable: Iterator[T]

    def __init__(self, iterable=None, length=None, label=None, item_show_func=None, interval=None,
                 file=None, writer=None, **kwargs):
        if iterable is not None:
            try:
                length = len(iterable)
            except TypeError:
                pass
        if interval is None:
            if length is None:
                interval = 1000
            else:
                interval = int(length * 0.05)
        if interval == 0:
            interval = 1
        if label is None:
            label = ''
        if file is None:
            file = click.get_text_stream('stderr')
        self.iterable = iterable
        self.length = length
        self.item_show_func = item_show_func
        self.current_item = None
        self.count = 0
        self.interval = interval
        self.last_update = 0
        self.label = label
        self.file = file
        self.writer = writer

    def update(self, n_steps: int, current_item=None):
        self.count += n_steps
        self.current_item = current_item
        if self.count > self.last_update:
            self._log()
            self.last_update += self.interval * (int(n_steps // self.interval) + 1)

    def render_progress(self):
        self._log()

    def _log(self):
        item = self.current_item
        if self.item_show_func is not None:
            show = self.item_show_func(item)
        else:
            show = str(item)
        i = self.count
        label = self.label
        if self.length is not None and self.length != 0:
            prog_label = "%s %d/%d (%0.2f%%)" % (label, i, self.length, i * 100.0 / self.length)
        else:
            prog_label = "%s %d" % (label, i)
        message = fmt_msg("%s: %s" % (prog_label, show))
        if self.writer is None:
            click.echo(message, file=self.file, err=True)
        else:
            self.writer(message)

    def __iter__(self):
        item: T
        for item in self.iterable:
            self.update(1, current_item=item)
            yield item

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


def progress(iterable: Iterator[T]=None, *args, **kwargs) -> Union[ProgressBar[T], ProgressLogger[T]]:
    """
    A wrapper that will dispatch to :func:`click.progressbar` when `sys.stdout` is a
    TTY and :class:`ProgressLogger` otherwise.
    """
    if sys.stdout.isatty():
        return click.progressbar(iterable=iterable, *args, **kwargs)
    else:
        return ProgressLogger(iterable=iterable, *args, **kwargs)


def progress_iter(iterable: Iterator[T]=None, *args, **kwargs) -> Iterator[T]:
    prog = progress(iterable, *args, **kwargs)
    with prog:
        yield from prog


def type_cast(value: str) -> Union[str, int, float, list, dict]:
    try:
        return json.loads(value)
    except (TypeError, ValueError):
        try:
            return ast.literal_eval(value)
        except (TypeError, ValueError):
            return str(value)


_MP_CONFIG = False


def configure_mp_context():
    global _MP_CONFIG
    if _MP_CONFIG:
        return
    multiprocessing.freeze_support()
    current_method = multiprocessing.get_start_method()
    if platform.system() not in ("Windows", "Darwin"):
        if current_method != 'forkserver':
            multiprocessing.set_start_method('forkserver')
    else:
        if current_method != 'spawn':
            multiprocessing.set_start_method('spawn')
    _MP_CONFIG = True
