import sys

import click
from six import string_types as basestring

import ms_deisotope
from ms_deisotope.data_source.metadata import activation as activation_module, data_transformation
from ms_deisotope.output import MzMLSerializer, MGFSerializer
from ms_peak_picker.scan_filter import parse as parse_filter
from pyteomics.xml import unitfloat


@click.group()
def ms_conversion():
    pass


def to_mgf(reader, outstream, msn_filters=None):
    if not msn_filters:
        msn_filters = []
    reader.make_iterator(grouped=False)
    writer = MGFSerializer(outstream, deconvoluted=False)
    for scan in reader:
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


def to_mzml(reader, outstream, pick_peaks=False, ms1_filters=None, msn_filters=None, default_activation=None):
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
    for bunch in reader:
        if bunch.precursor is not None:
            if ms1_filters:
                bunch = bunch._replace(precursor=bunch.precursor.transform(ms1_filters))
            if (pick_peaks or not bunch.precursor.is_profile):
                bunch.precursor.pick_peaks()
        for i, product in enumerate(bunch.products):
            if msn_filters:
                product = bunch.products[i] = product.transform(msn_filters)
            if pick_peaks or not product.is_profile:
                product.pick_peaks()
            if product.activation is None and default_activation is not None:
                product.activation = default_activation
        writer.save_scan_bunch(bunch)
    writer.complete()
    writer.format()


@ms_conversion.command("mzml", short_help="Convert a mass spectrometry data file to mzML")
@click.argument("source")
@click.argument("output", type=click.File(mode='wb'))
@click.option("-r", "--ms1-filter", "ms1_filters", multiple=True, type=parse_filter)
@click.option("-rn", "--msn-filter", "msn_filters", multiple=True, type=parse_filter)
@click.option("-p", "--pick-peaks", is_flag=True, help=("Enable peak picking, centroiding profile data"))
def mzml(source, output, ms1_filters=None, msn_filters=None, pick_peaks=False):
    reader = ms_deisotope.MSFileLoader(source)
    to_mzml(reader, output, pick_peaks=pick_peaks, ms1_filters=ms1_filters, msn_filters=msn_filters)


if __name__ == '__main__':
    import traceback

    def info(type, value, tb):
        if hasattr(sys, 'ps1') or not sys.stderr.isatty():
            sys.__excepthook__(type, value, tb)
        else:
            import ipdb
            traceback.print_exception(type, value, tb)
            ipdb.post_mortem(tb)

    sys.excepthook = info

    ms_conversion.main()
