import os
import gzip
import pickle
import sys
import logging

try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass

from ms_deisotope.config import get_config_dir
from ms_deisotope.task.log_utils import ColoringFormatter
from urllib.request import urlopen

_data_path = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "test_data"))

_http_data_uri_prefix = "https://raw.githubusercontent.com/mobiusklein/ms_deisotope/master//test_data/"


def _make_logger():
    format_string = '[%(asctime)s] %(levelname).1s | %(name)s | %(message)s'
    colorized_formatter = ColoringFormatter(format_string, datefmt="%H:%M:%S")
    logger = logging.getLogger(__name__)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(colorized_formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    return logger


logger = _make_logger()


def datafile(name: str, force_download: bool = False) -> os.PathLike:
    path = os.path.join(_data_path, name)
    if not os.path.exists(path) or force_download:
        config_dir = get_config_dir()
        test_data_dir = os.path.join(config_dir, "test_data")
        if not os.path.exists(test_data_dir):
            logger.info(
                "Creating user-level test data directory %s", test_data_dir)
            os.makedirs(test_data_dir)
        user_path = os.path.join(test_data_dir, name)
        if not os.path.exists(user_path):
            logger.info(
                "Downloading %s to %s from the remote repository", name, user_path)
            size = 0
            with open(user_path, 'wb') as fh:
                data = urlopen(_http_data_uri_prefix + name).read()
                size = len(data)
                fh.write(data)
            logger.info("Wrote %0.2f MB to %s", size / 1e6, user_path)
        path = user_path
    return path


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
