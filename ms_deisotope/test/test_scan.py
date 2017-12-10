import unittest

import numpy as np

from ms_deisotope.data_source import common, mzml
from ms_deisotope.averagine import peptide

from ms_peak_picker import FittedPeak
from ms_peak_picker.peak_statistics import gaussian_shape


points = [(576.5, 3, 8e4), (862.1, 2, 15e4)]
fwhm = 0.05


def make_profile(points, fwhm):
    peaks = []
    i = 0
    for point in points:
        tid = peptide.isotopic_cluster(point[0], point[1], truncate_after=0.99)
        for tp in tid:
            fp = FittedPeak(tp.mz, tp.intensity * point[2], 0, i, i,
                            fwhm, tp.intensity * point[2])
            peaks.append(fp)
    mz = np.array([0])
    intensity = np.array([0])

    for p in peaks:
        x, y = gaussian_shape(p)
        mz = np.concatenate([mz, [x[0] - 0.0001], x, [x[-1] + 0.0001]])
        intensity = np.concatenate([intensity, [0], y, [0]])
    return mz, intensity


class TestScanMachinery(unittest.TestCase):
    def make_scan(self):
        mz, intensity = make_profile(points, fwhm)
        scan = common.Scan(
            {
                "id": "spam",
                "index": 0,
                "ms level": 1,
                "m/z array": mz,
                "intensity array": intensity,
                "profile spectrum": "",
                "positive scan": "",
                "scanList": {
                    "scan": [
                        {"scan start time": 0}
                    ]
                }
            },
            mzml.MzMLDataInterface())
        return scan

    def test_pick_peaks(self):
        scan = self.make_scan()
        self.assertIsNone(scan.peak_set)
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        self.assertEqual(len(scan.peak_set), 10)
        self.assertIsNotNone(scan.has_peak(576.5))

    def test_equality(self):
        scan = self.make_scan()
        scan2 = self.make_scan()
        self.assertEqual(scan, scan2)


if __name__ == '__main__':
    unittest.main()
