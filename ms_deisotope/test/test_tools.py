import unittest
import io

from click.testing import CliRunner

from ms_deisotope.data_source.mzml import MzMLLoader
from ms_deisotope.tools import indexing, conversion
from ms_deisotope.test.common import datafile


def test_describe():
    runner = CliRunner()

    path = datafile("small.mzML")
    result = runner.invoke(indexing.describe, [path])
    lines = result.output.splitlines()
    assert "small.mzML" in lines[0]
    assert lines[1] == "File Format: mzML format"
    assert lines[2] == "ID Format: Thermo nativeID format"
    assert lines[3] == "Format Supports Random Access: True"


def test_mgf():
    runner = CliRunner()

    path = datafile("small.mzML")
    result = runner.invoke(conversion.mgf, [path, '-'])
    lines = result.output.splitlines()
    count = 0
    for line in lines:
        if "BEGIN" in line:
            count += 1
    assert count == 34


def test_mzml():
    runner = CliRunner()

    path = datafile("small.mzML")
    result = runner.invoke(conversion.mzml, ['-p', '-c', path, '-'])
    buff = io.BytesIO(result.output.encode("utf-8"))
    reader = MzMLLoader(buff)
    n = len(reader)
    assert n == 48
