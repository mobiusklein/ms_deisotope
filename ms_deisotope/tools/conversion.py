'''A collection of utilities for using :mod:`ms_deisotope` to convert between mass spectrometry
data file formats.
'''

import click
from six import string_types as basestring

from ms_peak_picker.scan_filter import parse as parse_filter
from pyteomics.xml import unitfloat

import ms_deisotope
from ms_deisotope.data_source._compression import GzipFile
from ms_deisotope.data_source.metadata import activation as activation_module, data_transformation
from ms_deisotope.output import MzMLSerializer, MGFSerializer

from ms_deisotope.tools.utils import is_debug_mode, register_debug_hook


@click.group()
def ms_conversion():
    """A command line tool for converting between mass spectrometry data formats
    """
    pass


def to_mgf(reader, outstream, msn_filters=None):
    """Translate the spectra from `reader` into MGF format written to `outstream`.

    As MGF files do not usually contain MS1 spectra, these will be omitted. Additionally,
    MSn spectra will be centroided if they are not already.

    Parameters
    ----------
    reader : :class:`~.ScanIterator`
        The source of spectra to iterate over
    outstream : file-like
        The output stream to write to
    msn_filters : list, optional
        An optional list of strings or :class:`~.ScanFilterBase` instances which will be
        used to transform the m/z and intennsity arrays of MSn spectra before they are futher
        procssed. (the default is None, which results in no transformations)

    """
    if not msn_filters:
        msn_filters = []
    reader.make_iterator(grouped=False)
    writer = MGFSerializer(outstream, deconvoluted=False)
    try:
        n_spectra = len(reader)
    except TypeError:
        n_spectra = None
    progbar = click.progressbar(
        reader,
        label="Processed Spectra", length=n_spectra,
        item_show_func=lambda x: str(x.id) if x else '')
    with progbar:
        for scan in progbar:
            if scan.ms_level == 1:
                continue
            if msn_filters:
                scan = scan.transform(msn_filters)
            if scan.peak_set is None:
                scan.pick_peaks()
            writer.save_scan(scan)


@ms_conversion.command('mgf', short_help="Convert a mass spectrometry data file to MGF")
@click.argument("source")
@click.argument("output", type=click.File(mode='wb'))
@click.option("-rn", "--msn-filter", "msn_filters", multiple=True, type=parse_filter)
def mgf(source, output, msn_filters=None):
    """Convert a mass spectrometry data file to MGF. MGF can only represent centroid spectra
    and generally does not contain any MS1 information.
    """
    reader = ms_deisotope.MSFileLoader(source)
    to_mgf(reader, output, msn_filters=msn_filters)


def to_mzml(reader, outstream, pick_peaks=False, ms1_filters=None, msn_filters=None, default_activation=None,
            correct_precursor_mz=False):
    """Translate the spectra from `reader` into mzML format written to `outstream`.

    Wraps the process of iterating over `reader`, performing a set of simple data transformations if desired,
    and then converts each :class:`~.Scan` into mzML format. Any data transformations are described in the
    appropriate metadata section. All other metadata from `reader` is copied to into `outstream`.

    Parameters
    ----------
    reader : :class:`~.ScanIterator`
        The source of spectra to iterate over
    outstream : file-like
        The output stream to write mzML to.
    pick_peaks : bool, optional
        Whether to centroid profile spectra (the default is False)
    ms1_filters : list, optional
        An optional list of strings or :class:`~.ScanFilterBase` instances which will be
        used to transform the m/z and intensity arrays of MS1 spectra before they are further
        processed (the default is None, which results in no transformations)
    msn_filters : list, optional
        An optional list of strings or :class:`~.ScanFilterBase` instances which will be
        used to transform the m/z and intennsity arrays of MSn spectra before they are futher
        procssed. (the default is None, which results in no transformations)
    default_activation : :class:`str` or :class:`dict`, optional
        A default activation type to use when `reader` does not contain that information (the default is None)
    correct_precursor_mz : bool, optional
        Whether or not to assign the precursor m/z of each product scan to the nearest peak
        m/z in the precursor's peak list. (the default is False, which results in no correction)

    """
    if ms1_filters is None:
        ms1_filters = []
    if msn_filters is None:
        msn_filters = []
    reader.make_iterator(grouped=True)
    writer = MzMLSerializer(outstream, len(reader), deconvoluted=False)
    writer.copy_metadata_from(reader)
    method = data_transformation.ProcessingMethod(software_id='ms_deisotope_1')
    if pick_peaks:
        method.add('MS:1000035')
    if correct_precursor_mz:
        method.add('MS:1000780')
    method.add('MS:1000544')
    writer.add_data_processing(method)
    if default_activation is not None:
        if isinstance(default_activation, basestring):
            default_activation = activation_module.ActivationInformation(
                default_activation, unitfloat(0, 'electronvolt'))
        elif isinstance(default_activation, dict):
            default_activation = activation_module.ActivationInformation(**default_activation)
        else:
            click.secho("Could not convert %r into ActivationInformation" % (default_activation, ), err=1,
                        fg='yellow')
            default_activation = None
    if pick_peaks:
        try:
            writer.remove_file_contents("profile spectrum")
        except KeyError:
            pass
        writer.add_file_contents("centroid spectrum")
    n_spectra = len(reader)
    progbar = click.progressbar(
        label="Loading Spectra", length=n_spectra,
        item_show_func=lambda x: str(x.precursor.id if x.precursor else
                                     x.products[0].id) if x else '')
    with progbar:
        for bunch in reader:
            progbar.current_item = bunch
            progbar.update((bunch.precursor is not None) + len(bunch.products))
            discard_peaks = False
            if bunch.precursor is not None:
                if ms1_filters:
                    bunch = bunch._replace(precursor=bunch.precursor.transform(ms1_filters))
                if (pick_peaks or not bunch.precursor.is_profile):
                    bunch.precursor.pick_peaks()
                if correct_precursor_mz:
                    if not pick_peaks:
                        bunch.precursor.pick_peaks()
                        if bunch.precursor.is_profile:
                            discard_peaks = True

            for i, product in enumerate(bunch.products):
                if msn_filters:
                    product = bunch.products[i] = product.transform(msn_filters)
                if pick_peaks or not product.is_profile:
                    product.pick_peaks()
                if product.activation is None and default_activation is not None:
                    product.activation = default_activation
                if correct_precursor_mz:
                    product.precursor_information.correct_mz()
            if discard_peaks:
                bunch.precursor.peak_set = None
            writer.save_scan_bunch(bunch)
    writer.complete()
    writer.format()


@ms_conversion.command("mzml", short_help="Convert a mass spectrometry data file to mzML")
@click.argument("source")
@click.argument("output", type=click.Path(writable=True))
@click.option("-r", "--ms1-filter", "ms1_filters", multiple=True, type=parse_filter)
@click.option("-rn", "--msn-filter", "msn_filters", multiple=True, type=parse_filter)
@click.option("-p", "--pick-peaks", is_flag=True, help=("Enable peak picking, centroiding profile data"))
@click.option("-z", "--compress", is_flag=True, help=("Compress the output file using gzip"))
@click.option("-c", "--correct-precursor-mz", is_flag=True, help=(
    "Adjust the precursor m/z of each MSn scan to the nearest peak m/z in the precursor"))
def mzml(source, output, ms1_filters=None, msn_filters=None, pick_peaks=False, compress=False,
         correct_precursor_mz=False):
    """Convert `source` into mzML format written to `output`, applying a collection of optional data
    transformations along the way.
    """
    reader = ms_deisotope.MSFileLoader(source)
    if compress:
        if not output.endswith(".gz"):
            output += '.gz'
        stream = click.open_file(output, 'wb')
        stream = GzipFile(fileobj=stream, mode='wb')
    else:
        stream = click.open_file(output, 'wb')
    with stream:
        to_mzml(reader, stream, pick_peaks=pick_peaks, ms1_filters=ms1_filters,
                msn_filters=msn_filters, correct_precursor_mz=correct_precursor_mz)


if is_debug_mode():
    register_debug_hook()

if __name__ == '__main__':
    click.secho("Running Debug Mode", fg='yellow')
    register_debug_hook()
    ms_conversion.main()
