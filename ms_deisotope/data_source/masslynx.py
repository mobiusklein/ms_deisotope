
try:
    from ms_deisotope.data_source._vendor.masslynx.libload import (
        proxy, register_dll as register_waters_masslynx_dll, _register_dll)

    from ms_deisotope.data_source._vendor.masslynx.loader import (
        MassLynxRawLoader, is_waters_raw_dir,
        determine_if_available, infer_reader,
        IndexEntry, Cycle
    )
except ImportError as e:
    message = str(e)

    def is_waters_raw_dir(path):
        return False

    def infer_reader(path):
        raise ValueError(message)

    def register_waters_masslynx_dll(paths):
        print("no-op: %s" % (message,))
        return False

    def determine_if_available():
        print("no-op: %s" % (message,))
        return False

    class MassLynxRawLoader(object):
        def __init__(self, *args, **kwargs):
            raise ImportError(message)
