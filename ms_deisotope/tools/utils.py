import multiprocessing
import os
import re
import sys
import warnings

import click

from brainpy import periodic_table

from ms_deisotope.averagine import (
    Averagine, glycan as n_glycan_averagine, permethylated_glycan,
    peptide, glycopeptide, heparin, heparan_sulfate)


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
            except ImportError
                import pdb as pdb_api
            traceback.print_exception(type, value, tb)
            pdb_api.post_mortem(tb)

    sys.excepthook = info


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
