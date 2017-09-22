import unittest

import numpy as np


def make_peak_set_test_suite(peak_set_cls, peak_cls):

    class TestDeconvolutedPeakSet(unittest.TestCase):
        def test_search(self):
            x = np.arange(1000, 1200, 0.001)
            y = np.ones_like(x)

            peaks = [peak_cls(x[i], y[i], 1, 1, None, 0) for i in range(len(x))]
            ps = peak_set_cls(peaks)
            ps.reindex()
            t = ps[5000].neutral_mass
            assert ps.has_peak(t) is ps[5000]

            for peak in ps.all_peaks_for(t):
                 assert abs((peak.neutral_mass - t) / t) < 1e-5
    return TestDeconvolutedPeakSet


try:
    from ms_deisotope._c.peak_set import (
        DeconvolutedPeakSet as CDeconvolutedPeakSet,
        DeconvolutedPeak as CDeconvolutedPeak,
        DeconvolutedPeakSetIndexed as CDeconvolutedPeakSetIndexed)

    TestCDeconvolutedPeakSet = make_peak_set_test_suite(CDeconvolutedPeakSet, CDeconvolutedPeak)
    TestCDeconvolutedPeakSetIndexed = make_peak_set_test_suite(CDeconvolutedPeakSetIndexed, CDeconvolutedPeak)

    from ms_deisotope.peak_set import (
        _DeconvolutedPeak, _DeconvolutedPeakSet)

    TestPythonDeconvolutedPeakSet = make_peak_set_test_suite(_DeconvolutedPeakSet, _DeconvolutedPeak)

except ImportError:
    from ms_deisotope.peak_set import (
        _DeconvolutedPeak, _DeconvolutedPeakSet)


    TestPythonDeconvolutedPeakSet = make_peak_set_test_suite(_DeconvolutedPeakSet, _DeconvolutedPeak)


if __name__ == '__main__':
    unittest.main()
