'''Convert a CSV file of peak centroids to a single spectrum mzML file.

Useful for when you need to convert some scattering of points from a viewer's
text export into something that a search engine can read.
'''
import sys
import os

from ms_deisotope.data_source import text
from ms_deisotope.data_source.metadata import file_information

from ms_deisotope.output.mzml import MzMLSerializer

csv_path = sys.argv[1]

scan = text.scan_from_csv(open(csv_path, 'rt'), is_profile=False)

dest = os.path.splitext(csv_path)[0] + '.mzML'

with MzMLSerializer(open(dest, 'wb'), 1, deconvoluted=False,
                    sample_name=os.path.basename(csv_path), build_extra_index=False) as writer:
    writer.add_file_contents(file_information.MS_MSn_Spectrum.name)
    writer.add_data_processing(writer.build_processing_method())
    writer.save(scan)
