import itertools
import multiprocessing
import os
import re
import sys
import warnings
import logging

import click

from brainpy import periodic_table

from ms_deisotope.averagine import (
    Averagine, glycan as n_glycan_averagine, permethylated_glycan,
    peptide, glycopeptide, heparin, heparan_sulfate)

from ms_deisotope.task.log_utils import fmt_msg


def processes_option(f):
    opt = click.option(
        "-p", "--processes", 'processes', type=click.IntRange(1, multiprocessing.cpu_count()),
        default=min(multiprocessing.cpu_count(), 4), help=('Number of worker processes to use. Defaults to 4 '
                                                           'or the number of CPUs, whichever is lower'))
    return opt(f)


def validate_element(element):
    valid = element in periodic_table
    if not valid:
        raise click.Abort("%r is not a valid element" % element)
    return valid


def parse_averagine_formula(formula):
    if isinstance(formula, Averagine):
        return formula
    return Averagine({k: float(v) for k, v in re.findall(r"([A-Z][a-z]*)([0-9\.]*)", formula)
                      if float(v or 0) > 0 and validate_element(k)})


averagines = {
    'glycan': n_glycan_averagine,
    'permethylated-glycan': permethylated_glycan,
    'peptide': peptide,
    'glycopeptide': glycopeptide,
    'heparin': heparin,
    "heparan-sulfate": heparan_sulfate
}


def validate_averagine(averagine_string):
    if isinstance(averagine_string, Averagine):
        return averagine_string
    if averagine_string in averagines:
        return averagines[averagine_string]
    else:
        return parse_averagine_formula(averagine_string)


class AveragineParamType(click.types.StringParamType):
    name = "MODEL"

    models = averagines

    def convert(self, value, param, ctx):
        return validate_averagine(value)

    def get_metavar(self, param):
        return '[%s]' % '|'.join(sorted(averagines.keys()))

    def get_missing_message(self, param):
        return 'Choose from %s, or provide a formula.' % ', '.join(self.choices)


def register_debug_hook():
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
    logging.basicConfig(level="DEBUG")


def is_debug_mode():
    env_val = os.environ.get('MS_DEISOTOPE_DEBUG', '').lower()
    if not env_val:
        return False
    if env_val in ('0', 'no', 'false', 'off'):
        return False
    elif env_val in ('1', 'yes', 'true', 'on'):
        return True
    else:
        warnings.warn("MS_DEISOTOPE_DEBUG value %r was not recognized. Enabling debug mode" % (env_val, ))
        return True


class Spinner(object):
    """Super-simple synchronous CLI spinner implementation.
    """
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


class ProgressLogger(object):
    """A simple text logger that wraps an iterable and logs update messages as chunks are requested.

    This class tries to emulate :func:`click.progressbar` for use when the attached STDERR stream is not
    a terminal.
    """
    def __init__(self, iterable=None, length=None, label=None, item_show_func=None, interval=None, file=None, **kwargs):
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

    def update(self, n):
        self.count += n
        if self.count > self.last_update:
            self._log()
            self.last_update += self.interval * (int(n // self.interval) + 1)

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
        click.echo(message, file=self.file, err=True)

    def __iter__(self):
        for item in self.iterable:
            self.current_item = item
            self.update(1)
            yield item

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass


def progress(*args, **kwargs):
    """A wrapper that will dispatch to :func:`click.progressbar` when `sys.stdout` is a TTY and :class:`ProgressLogger`
    otherwise.
    """
    if sys.stdout.isatty():
        return click.progressbar(*args, **kwargs)
    else:
        return ProgressLogger(*args, **kwargs)
