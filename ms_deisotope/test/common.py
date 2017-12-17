import os
import gzip
import pickle
import sys


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
