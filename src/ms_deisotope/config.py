'''Simple library-wide configuration management to handle
the tracking of external 3rd-party libraries and machine state.

'''
import os
import sys
import copy
import json
import warnings

from six import PY2


CONFIG_FILE_NAME = 'config.json'


def _get_home_dir():
    """Find user's home directory if possible.
    Otherwise, returns None.

    :see:
        http://mail.python.org/pipermail/python-list/2005-February/325395.html
    """
    if PY2 and sys.platform == 'win32':
        path = os.path.expanduser(b"~").decode(sys.getfilesystemencoding())
    else:
        path = os.path.expanduser("~")
    if os.path.isdir(path):
        return path
    for evar in ('HOME', 'USERPROFILE', 'TMP'):
        path = os.environ.get(evar)
        if path is not None and os.path.isdir(path):
            return path
    return None


def get_config_dir():
    """Get the configuration directory path.

    Tries the following routes:
        1. The environment variable "MS_DEISOTOPE_CONFIGDIR"
        2. If on an XDG-compliant platform (Linux, Free BSD), the environment
           variable "XDG_CONFIG_HOME"/ms_deisotope
        3. The user's home directory/.ms_deisotope

    If the configuration directory does not exist, it will
    be created.

    Returns
    -------
    str
    """
    FALLBACK_DIR = '.'
    configdir = os.environ.get('MS_DEISOTOPE_CONFIGDIR')
    if configdir is not None:
        configdir = os.path.abspath(configdir)
        if not os.path.exists(configdir):
            try:
                os.makedirs(configdir)
            except OSError:
                return FALLBACK_DIR
        return configdir

    home_dir = _get_home_dir()
    if home_dir is not None:
        configdir = os.path.join(home_dir, '.ms_deisotope')
    else:
        configdir = None
    if sys.platform.startswith(('linux', 'freebsd')):
        configdir = None
        xdg_base = os.environ.get('XDG_CONFIG_HOME')
        if xdg_base is None:
            xdg_base = home_dir
            if xdg_base is not None:
                xdg_base = os.path.join(xdg_base, '.config')
        if xdg_base is not None:
            configdir = os.path.join(xdg_base, 'ms_deisotope')
    if configdir is not None:
        if os.path.exists(configdir):
            return configdir
        try:
            os.makedirs(configdir)
            return configdir
        except OSError:
            return FALLBACK_DIR
    return FALLBACK_DIR


_DEFAULT_CONFIG = {
    'vendor_readers': {
        "thermo-com": [],
        "thermo-net": [],
        "agilent-com": [],
        "waters-masslynx": [],
    },
    'schema_version': '1.0.0'
}


def get_config():
    """Load the config.json file from the configuration
    directory given by :func:`get_config_dir`.

    If the config file does not exist, a default one
    will be created.

    Returns
    -------
    dict
    """
    confdir = get_config_dir()
    path = os.path.join(confdir, CONFIG_FILE_NAME)
    if not os.path.exists(path):
        try:
            save_config(_DEFAULT_CONFIG)
        except OSError as err:
            warnings.warn("Encountered the error %r when trying to write the default configuration" % (err, ))
            return copy.deepcopy(_DEFAULT_CONFIG)
    try:
        with open(path, 'rt') as fh:
            config = json.load(fh)
    except OSError as err:
        warnings.warn(
            "Encountered the error %r when trying to read configuration" % (err, ))
        config = copy.deepcopy(_DEFAULT_CONFIG)
    if "schema_version" not in config:
        config['schema_version'] = '1.0.0'
    return config


def save_config(config=None):
    """Save the configuration in `config` to disk in the
    configuration directory given by :func:`get_config_dir`.

    Parameters
    ----------
    config : dict, optional
        The configuration dictionary (the default is None, which will use the default, empty config)

    """
    if config is None:
        config = copy.deepcopy(_DEFAULT_CONFIG)
    confdir = get_config_dir()
    path = os.path.join(confdir, CONFIG_FILE_NAME)
    with open(path, 'wt') as fh:
        json.dump(config, fh)
