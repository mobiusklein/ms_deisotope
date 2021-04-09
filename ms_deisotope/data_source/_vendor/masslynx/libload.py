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
    result = _register_dll(search_paths, override)
    if not result:
        raise ImportError("Could not load MassLynx SDK from search paths!")



def _register_dll(search_paths=None, override=True):
    if search_paths is None:
        search_paths = []
    elif isinstance(search_paths, str):
        search_paths = [search_paths]
    elif not isinstance(search_paths, list):
        raise TypeError("Expected a list of paths or a single path string")
    from ms_deisotope.config import get_config
    search_paths = list(search_paths)
    search_paths.extend(get_config().get('vendor_readers', {}).get('waters-masslynx', []))
    found = find_library("MassLynx.dll")
    if found:
        search_paths.append(found)
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
    '''Checks whether or not the Waters
    RAW directory reading feature is available.

    Returns
    -------
    :class:`bool`:
        Whether or not the feature is enabled.
    '''
    try:
        return _register_dll(override=False)
    except (OSError, ImportError):
        return False


class DLLProxy(object):
    def __getattr__(self, attribute):
        if dll is None:
            register_dll()
        return getattr(dll, attribute)

    def _is_loaded(self):
        return dll is not None

proxy = DLLProxy()
