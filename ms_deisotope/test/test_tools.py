import unittest

from click.testing import CliRunner

from ms_deisotope.tools import indexing

from ms_deisotope.test.common import datafile


def test_task():
    runner = CliRunner()

    path = datafile("small.mzML")
    result = runner.invoke(indexing.describe, [path])
    print(result.output)
    lines = result.output.splitlines()
    assert "small.mzML" in lines[0]
    assert lines[1] == "File Format: mzML format"
    assert lines[2] == "ID Format: Thermo nativeID format"
    assert lines[3] == "Format Supports Random Access: True"
