import csv
from pyteomics import mgf

from StringIO import StringIO


def serialize_peaklist_csv(scan_source, out=None):
    if out is None:
        out = StringIO()
    writer = csv.writer(out)
    writer.writerow(["scan_id", "mz", "charge", "intensity", "fwhm", "score", "a_to_a2_ratio"])
    for scan in scan_source:
        for peak in scan.deconvoluted_peak_set:
            writer.writerow(map(str, [scan.id, peak.mz, peak.charge, peak.intensity, peak.full_width_at_half_max,
                                   peak.score, peak.a_to_a2_ratio]))
    return out


def serialize_mgf(scan_source, out=None):
    pass
