import unittest
import io

from click.testing import CliRunner

from ms_deisotope.data_source.mzml import MzMLLoader
from ms_deisotope.data_source import _compression
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

    result = runner.invoke(conversion.mgf, [path, '-z', '-'])
    assert _compression.starts_with_gz_magic(result.stdout_bytes)
    buff = io.BytesIO(result.stdout_bytes)
    reader = _compression.GzipFile(fileobj=buff, mode='rb')
    count = 0
    for line in reader:
        if b"BEGIN" in line:
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

    result = runner.invoke(
        conversion.mzml, ['-p', '-z', '-c', path, '-'], catch_exceptions=False)
    buff = io.BytesIO(result.stdout_bytes)
    reader = MzMLLoader(_compression.get_opener(buff))
    n = len(reader)
    assert n == 48



def test_idzip():
    runner = CliRunner(mix_stderr=False)
    path = datafile("20150710_3um_AGP_001_29_30.mzML.gz")
    stdin_data = io.BytesIO(open(path, 'rb').read())
    result = runner.invoke(
        indexing.idzip_compression,
        ['-'],
        input=stdin_data,
        mix_stderr=False)
    assert b"Detected gzip input file" in result.stderr_bytes
    outbuff = io.BytesIO(result.stdout_bytes)
    outstream = _compression.GzipFile(fileobj=outbuff, mode='rb')
    instream = _compression.GzipFile(path, mode='rb')
    in_data = instream.read()
    out_data = outstream.read()
    assert in_data == out_data

    path = datafile("small.mzML")
    stdin_data = io.BytesIO(open(path, 'rb').read())
    result = runner.invoke(
        indexing.idzip_compression,
        ['-'],
        input=stdin_data,
        mix_stderr=False)
    assert b"Detected gzip input file" not in result.stderr_bytes
    outbuff = io.BytesIO(result.stdout_bytes)
    outstream = _compression.GzipFile(fileobj=outbuff, mode='rb')
    instream = io.open(path, mode='rb')
    in_data = instream.read()
    out_data = outstream.read()
    assert in_data == out_data
