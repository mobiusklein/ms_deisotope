import unittest
import os
import io
import tempfile

from click.testing import CliRunner

from ms_deisotope.data_source.mzml import MzMLLoader
from ms_deisotope.output import ProcessedMzMLDeserializer
from ms_deisotope.data_source import _compression
from ms_deisotope.tools import indexing, conversion, deisotoper
from ms_deisotope.test.common import datafile


def run_ms_deisotope():
    runner = CliRunner(mix_stderr=False)
    path = datafile("20150710_3um_AGP_001_29_30.mzML.gz")
    reference = datafile("20150710_3um_AGP_001_29_30.preprocessed.mzML.gz")
    result = runner.invoke(deisotoper.deisotope, [
        "-b", 0, "-t", 20, "-tn", 10, "-m", 3, "-mn", 1, path, reference
    ])
    print(result.stdout)


if __name__ == "__main__":
    run_ms_deisotope()