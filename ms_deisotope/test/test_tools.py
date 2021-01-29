import unittest
import os
import io
import tempfile

import pytest
from click.testing import CliRunner

from ms_deisotope.data_source.mzml import MzMLLoader
from ms_deisotope.output import ProcessedMzMLDeserializer
from ms_deisotope.data_source import _compression
from ms_deisotope.tools import indexing, conversion, deisotoper
from ms_deisotope.test.common import datafile


def test_describe():
    runner = CliRunner(mix_stderr=False)

    path = datafile("small.mzML")
    result = runner.invoke(indexing.describe, [path])
    lines = result.output.splitlines()
    assert "small.mzML" in lines[0]
    assert lines[1] == "File Format: mzML format"
    assert lines[2] == "ID Format: Thermo nativeID format"
    assert lines[3] == "Format Supports Random Access: True"


def test_mgf():
    runner = CliRunner(mix_stderr=False)
    if os.path.exists("-idx.json"):
        raise IOError("Orphan index file exists before running test")
    path = datafile("small.mzML")
    result = runner.invoke(conversion.mgf, [path, '-'], catch_exceptions=False)
    lines = result.output.splitlines()
    count = 0
    for line in lines:
        if "BEGIN" in line:
            count += 1
    assert count == 34
    if os.path.exists("-idx.json"):
        raise IOError("Orphan index file exists after running uncompressed test")
    result = runner.invoke(conversion.mgf, [path, '-z', '-'])
    assert _compression.starts_with_gz_magic(result.stdout_bytes)
    buff = io.BytesIO(result.stdout_bytes)
    reader = _compression.GzipFile(fileobj=buff, mode='rb')
    count = 0
    for line in reader:
        if b"BEGIN" in line:
            count += 1
    assert count == 34
    if os.path.exists("-idx.json"):
        raise IOError("Orphan index file exists after running compressed test")


def test_mzml():
    runner = CliRunner(mix_stderr=False)
    if os.path.exists("-idx.json"):
        raise IOError("Orphan index file exists before running test")
    path = datafile("small.mzML")
    result = runner.invoke(conversion.mzml, ['-p', '-c', path, '-'])
    buff = io.BytesIO(result.output.encode("utf-8"))
    reader = MzMLLoader(buff)
    n = len(reader)
    assert n == 48
    if os.path.exists("-idx.json"):
        raise IOError(
            "Orphan index file exists after running uncompressed test")

    result = runner.invoke(
        conversion.mzml, ['-p', '-z', '-c', path, '-'], catch_exceptions=False)
    buff = io.BytesIO(result.stdout_bytes)
    reader = MzMLLoader(_compression.get_opener(buff))
    n = len(reader)
    assert n == 48
    if os.path.exists("-idx.json"):
        raise IOError("Orphan index file exists after running compressed test")



def compare_peaks(peaks_a, peaks_b):
    missing = []
    for peak in peaks_a:
        peaks = peaks_b.all_peaks_for(
            peak.neutral_mass, 1e-6)
        for ref in peaks:
            if peak == ref:
                break
        else:
            missing.append(peak)
    return missing


def diff_deconvoluted_peak_set(peaks_a, peaks_b):
    a_missing = compare_peaks(peaks_a, peaks_b)
    b_missing = compare_peaks(peaks_b, peaks_a)
    return a_missing, b_missing


@pytest.mark.slow
def test_ms_deisotope():
    runner = CliRunner(mix_stderr=False)
    path = datafile("20150710_3um_AGP_001_29_30.mzML.gz")
    reference = datafile("20150710_3um_AGP_001_29_30.preprocessed.mzML.gz")
    outpath = tempfile.mktemp()
    result = runner.invoke(deisotoper.deisotope, [
        "-b", 0, "-t", 20, "-tn", 10, "-m", 3, "-mn", 1, path, outpath
    ])
    result_reader = ProcessedMzMLDeserializer(outpath)
    reference_reader = ProcessedMzMLDeserializer(_compression.get_opener(reference))
    assert len(result_reader) == len(reference_reader)
    for a_bunch, b_bunch in zip(result_reader, reference_reader):
        assert len(a_bunch.products) == len(b_bunch.products)
        aprec = a_bunch.precursor
        bprec = b_bunch.precursor
        assert aprec.id == bprec.id
        diffa, diffb = diff_deconvoluted_peak_set(
            aprec.deconvoluted_peak_set, bprec.deconvoluted_peak_set)
        assert len(aprec.deconvoluted_peak_set) == len(
            bprec.deconvoluted_peak_set), "Peak Counts Diff On %r, (%r, %r)" % (aprec.id, diffa, diffb)
        assert aprec.deconvoluted_peak_set == bprec.deconvoluted_peak_set, "Peaks Diff On %r, (%r, %r)" % (
            aprec.id, diffa, diffb)

        for aprod, bprod in zip(a_bunch.products, b_bunch.products):
            assert aprod.id == bprod.id
            diffa, diffb = diff_deconvoluted_peak_set(aprod.deconvoluted_peak_set, bprod.deconvoluted_peak_set)
            assert len(aprod.deconvoluted_peak_set) == len(
                bprod.deconvoluted_peak_set), "Peak Counts Diff On %r, (%r, %r)" % (aprod.id, diffa, diffb)
            assert aprod.deconvoluted_peak_set == bprod.deconvoluted_peak_set, "Peaks Diff On %r" % (
                aprod.id, diffa, diffb)

    result_reader.close()
    reference_reader.close()
    os.remove(outpath)


def test_idzip():
    runner = CliRunner(mix_stderr=False)
    path = datafile("20150710_3um_AGP_001_29_30.mzML.gz")
    stdin_data = io.BytesIO(open(path, 'rb').read())
    result = runner.invoke(
        indexing.idzip_compression,
        ['-'],
        input=stdin_data)
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
        input=stdin_data)
    assert b"Detected gzip input file" not in result.stderr_bytes
    outbuff = io.BytesIO(result.stdout_bytes)
    outstream = _compression.GzipFile(fileobj=outbuff, mode='rb')
    instream = io.open(path, mode='rb')
    in_data = instream.read()
    out_data = outstream.read()
    assert in_data == out_data
