import os
import gzip
import pickle
import sys

try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass

data_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "test_data"))


def datafile(name):
    return os.path.join(data_path, name)


def gzload(path):
    with gzip.open(path, 'rb') as fh:
        if sys.version_info.major > 2:
            return pickle.load(fh, encoding='latin1')
        else:
            return pickle.load(fh)


def example_scan_bunch():
    import ms_deisotope
    reader = ms_deisotope.MSFileLoader(
        datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
    return reader.next()
