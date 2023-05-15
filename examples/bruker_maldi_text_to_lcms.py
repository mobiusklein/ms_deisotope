"""
Convert a single spectrum TSV file exported from Bruker's MALDI software
into an mzML file with a false LC-MS peak shape.
"""
import os

import click

from ms_deisotope.data_source import text
from ms_deisotope.data_source.metadata import file_information

from ms_deisotope.output.mzml import MzMLSerializer


@click.command("csv_to_mzml", short_help="Convert an m/z->intensity MALDI TSV file into a profile LC-MS mzML file")
@click.argument("tsv_path", type=click.Path(readable=True))
@click.option("-n", "--number-of-points", type=int, default=5, required=False, help="The number of LC points")
@click.option("-p", "--polarity", type=click.Choice(['1', '-1']), default=1, help="The scan polarity to record. Default is 1 (positive mode)")
def main(tsv_path, number_of_points: int, polarity=1):
    """
    Convert a single spectrum TSV file exported from Bruker's MALDI software
    into an mzML file with a false LC-MS peak shape.

    Example:
        csv_to_mzml path/to/data.csv -n 10
    Writes to `path/to/data.mzML`
    """
    polarity = int(polarity)
    with open(tsv_path, 'rt') as fh:
        scan = text.scan_from_csv(
            fh,
            delimiter='\t',
            is_profile=True,
            ms_level=1,
            polarity=polarity)

    dest = os.path.splitext(tsv_path)[0] + '.mzML'

    with MzMLSerializer(open(dest, 'wb'), deconvoluted=False, n_spectra=number_of_points,
                        sample_name=os.path.basename(tsv_path), build_extra_index=False) as writer:
        writer.add_file_contents(file_information.MS_MS1_Spectrum.name)
        method = writer.build_processing_method(picked_peaks=False, baseline_reduction=True)
        writer.add_data_processing(method)
        n = number_of_points // 2
        for i in range(number_of_points):
            scan_i = scan.copy()
            scan_i.id = f"index={i + 1}"
            scan_i.scan_time = i * 0.1
            scan_i.arrays *= max(1 - (abs(n - i) / n + 0.1), 0.01)
            writer.save(scan_i)


if __name__ == "__main__":
    main()
