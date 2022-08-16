import logging

try:
    from ms_deisotope.data_source._vendor.AgilentD import (
        AgilentDLoader, register_dll_dir, AgilentDScanPtr,
        AgilentDDataInterface,
        log as _api_logger)
    from comtypes import COMError

    comtypes_logger = logging.getLogger("comtypes")
    comtypes_logger.setLevel("INFO")
    _api_logger.setLevel("INFO")

    def is_agilent_d_dir(path):
        try:
            AgilentDLoader(path)
            return True
        except (WindowsError, IOError, ImportError, COMError):
            return False

    def infer_reader(path):
        if is_agilent_d_dir(path):
            return AgilentDLoader
        raise ValueError("Not Agilent .d Directory")

    def determine_if_available():
        try:
            AgilentDLoader.create_com_object()
            return True
        except (WindowsError, COMError):
            return False
except ImportError as e:
    message = str(e)

    def is_agilent_d_dir(path):
        return False

    def infer_reader(path):
        raise ValueError(message)

    def register_dll_dir(paths):
        print("no-op: %s" % (message,))
        return False

    def determine_if_available():
        print("no-op: %s" % (message,))
        return False

    class AgilentDLoader(object):
        def __init__(self, *args, **kwargs):
            raise ImportError(message)
