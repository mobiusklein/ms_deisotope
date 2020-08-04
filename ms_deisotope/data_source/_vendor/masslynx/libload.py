'''
The Windows DLL needed for the MassLynx reader to work can be obtained from:
https://interface.waters.com/masslynx/developers-area/sdks/, and by linking
with them we agree to follow the associated licensing.
'''

import os
from ctypes import WinDLL
from ctypes.util import find_library

dll = None


def register_dll(search_paths=None, override=True):
    if search_paths is None:
        search_paths = []
    elif isinstance(search_paths, str):
        search_paths = [search_paths]
    from ms_deisotope.config import get_config
    search_paths.extend(get_config().get('vendor_readers', {}).get('waters-masslynx', []))
    search_paths.append(find_library("MassLynx.dll"))
    global dll
    if dll is None or override:
        for lib_path in search_paths:
            try:
                dll = _load_library(lib_path)
            except (Exception) as err:
                continue
            if dll is not None:
                break
    return dll


def _load_library(lib_path):
    dll = WinDLL(os.path.realpath(lib_path))
    return dll


def determine_if_available():
    return register_dll() is not None


class DLLProxy(object):

    def __getattribute__(self, attribute):
        if dll is None:
            register_dll()
        return getattr(dll, attribute)

    def __getattr__(self, attribute):
        if dll is None:
            register_dll()
        return getattr(dll, attribute)

proxy = DLLProxy()
