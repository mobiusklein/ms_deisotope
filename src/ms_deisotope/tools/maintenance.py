import os
import click

from ms_deisotope.config import get_config_dir, get_config, save_config
from ms_deisotope.tools.utils import register_debug_hook


@click.group()
def maintenance():
    '''Tools for checking features and managing optional dependencies
    of ms_deisotope.
    '''
    pass


@maintenance.command('vendor-readers-available', short_help="Check if vendor readers are available")
def has_vendor_readers():
    """Logs whether libraries for using mass spectrometry vendor libraries to read
    directly from their native file formats are installed.
    """
    from ms_deisotope.data_source.agilent_d import determine_if_available as agilent_d_available
    from ms_deisotope.data_source.thermo_raw import determine_if_available as thermo_com_available
    from ms_deisotope.data_source.thermo_raw_net import determine_if_available as thermo_net_available
    from ms_deisotope.data_source.masslynx import determine_if_available as waters_masslynx_available

    try:
        import comtypes
        click.secho('comtypes %s detected' % (comtypes.__version__, ))
    except ImportError:
        click.secho("comtypes not installed", fg='yellow')
    try:
        import clr
        click.secho("pythonnet %s detected" % (clr.__version__, ))
    except ImportError:
        click.secho("pythonnet not installed", fg='yellow')

    click.echo("Thermo COM Reader Available: %s" % (thermo_com_available(), ))
    click.echo("Thermo RawFileReader Available: %s" % (thermo_net_available(), ))
    click.echo("Agilent COM Reader Available: %s" % (agilent_d_available(), ))
    click.echo("Waters MassLynx Available: %s" % (bool(waters_masslynx_available()), ))


@maintenance.command('register-thermo-com', short_help="Register Thermo COM Reader Library")
@click.option('--path', type=click.Path(file_okay=True, dir_okay=False), help="Specify a DLL to try to register", multiple=True)
def register_thermo_com(path=None):
    '''Register an installed MSFileReader or XRawFile DLL. Searches the standard installation
    paths, but the `--path` option can be used to specify additional search paths.
    '''
    if path is None:
        path = []
    try:
        import comtypes
    except ImportError:
        click.secho("comtypes not installed", fg='yellow')
        raise click.ClickException(
            "This feature requires the comtypes library. Please install it.")
    import logging
    from ms_deisotope.data_source._vendor.MSFileReader import _register_dll, log, DLL_IS_LOADED
    from ms_deisotope.data_source.thermo_raw import determine_if_available
    if DLL_IS_LOADED:
        click.secho("DLL Already Registered and Loaded")
        return
    if log is not None:
        log.addHandler(logging.StreamHandler())
        log.setLevel('DEBUG')
    result = _register_dll(path)
    if result:
        click.secho("DLL Registration Successful")
        if determine_if_available():
            click.secho("DLL Load Succesful", fg='cyan')
            if path is not None:
                config = get_config()
                config['vendor_readers']['thermo-com'].extend(path)
                save_config(config)
        else:
            click.secho("DLL Load Unsuccessful", fg='red')
    else:
        click.secho("DLL Registration Unsuccessful")
    if log is not None:
        log.handlers = log.handlers[:-1]


@maintenance.command('register-thermo-net', short_help="Register Thermo .NET Reader Library")
@click.option('--path', type=click.Path(file_okay=False, dir_okay=True),
              help="Specify a Directory of DLLs to try to register", multiple=True)
def register_thermo_net(path):
    '''Register a bundle of Thermo's RawFileReader library's DLLs. The `--path` option can
    be used to specify additional search paths.
    '''
    if path is None:
        path = []
    path = list(map(os.path.abspath, path))
    try:
        import clr
    except ImportError:
        click.secho("pythonnet not installed", fg='yellow')
        raise click.ClickException(
            "This feature requires the pythonnet library. Please install it.")
    from ms_deisotope.data_source.thermo_raw_net import _register_dll, _RawFileReader
    from ms_deisotope.data_source.thermo_raw_net import determine_if_available
    if _RawFileReader:
        click.secho("DLL Already Registered and Loaded")
        return
    result = _register_dll(path)
    if result:
        click.secho("DLL Registration Successful")
        if determine_if_available():
            click.secho("DLL Load Succesful", fg='cyan')
            if path is not None:
                config = get_config()
                rest = sorted(set(path) - set(config['vendor_readers']['thermo-net']))
                config['vendor_readers']['thermo-net'].extend(rest)
                save_config(config)
        else:
            click.secho("DLL Load Unsuccessful", fg='red')
    else:
        click.secho("DLL Registration Unsuccessful")


@maintenance.command('register-waters-masslynx', short_help="Register Waters MassLynx SDK")
@click.option('--path', type=click.Path(file_okay=True, dir_okay=False),
              help="Specify the MassLynx.dll to try to register", multiple=True)
def register_waters_masslynx(path):
    if path is None:
        path = []
    from ms_deisotope.data_source._vendor.masslynx import libload
    path = list(map(os.path.abspath, path))
    result = libload._register_dll(path)
    if result:
        click.secho("DLL Registration Successful")
        out = libload.determine_if_available()
        if out:
            click.secho("DLL Load Succesful", fg='cyan')
            if path is not None:
                config = get_config()
                config.setdefault('vendor_readers', {})
                config['vendor_readers'].setdefault("waters-masslynx", [])
                rest = sorted(
                    filter(bool, set(path) - set(config.get('vendor_readers', {}).get('waters-masslynx', []))))
                config['vendor_readers']['waters-masslynx'].extend(rest)
                save_config(config)
        else:
            click.secho("DLL Load Unsuccessful", fg='red')
    else:
        click.secho("DLL Registration Unsuccessful")


@maintenance.command('show-config', short_help="Display the config file's contents")
def show_config():
    '''Load the config file and write it to STDOUT
    '''
    import json
    click.echo(get_config_dir())
    config = get_config()
    click.echo(json.dumps(config, sort_keys=True, indent=2))


if __name__ == "__main__":
    register_debug_hook()
    maintenance.main()
