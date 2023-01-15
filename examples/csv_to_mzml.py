"""Convert a CSV file of peak centroids to a single spectrum mzML file.

Useful for when you need to convert some scattering of points from a viewer's
text export into something that a search engine can read.
"""
import sys
import os

import click

from ms_deisotope.data_source import text, PrecursorInformation
from ms_deisotope.data_source.metadata import file_information

from ms_deisotope.output.mzml import MzMLSerializer

@click.command("csv_to_mzml", short_help="Convert an m/z,intensity CSV file into a centroided mzML file")
@click.argument("csv_path", type=click.Path(readable=True))
@click.option("-m", "--precursor-mz", type=float, default=None, required=False, help="The selected ion m/z to record")
@click.option("-z", "--precursor-charge", type=int, default=None, required=False, help="The selected ion charge to record")
@click.option("-p", "--polarity", type=click.Choice(['1', '-1']), default=1, help="The scan polarity to record. Default is 1 (positive mode)")
def main(csv_path, precursor_mz=None, precursor_charge=None, polarity=1):
    """Convert an m/z,intensity CSV file into a centroided mzML file.

    Example:
        csv_to_mzml path/to/data.csv -m 571.63 -z -3 -p -1
    Writes to `path/to/data.mzML`
    """
    polarity = int(polarity)
    with open(csv_path, 'rt') as fh:
        scan = text.scan_from_csv(
            fh, is_profile=False, precursor_mz=precursor_mz, precursor_charge=precursor_charge,
            polarity=polarity)

    scan.pick_peaks()

    dest = os.path.splitext(csv_path)[0] + '.mzML'

    with MzMLSerializer(open(dest, 'wb'), 1, deconvoluted=False,
                        sample_name=os.path.basename(csv_path), build_extra_index=False) as writer:
        writer.add_file_contents(file_information.MS_MSn_Spectrum.name)
        writer.add_data_processing(writer.build_processing_method())
        writer.save(scan)

if __name__ == "__main__":
    main()
