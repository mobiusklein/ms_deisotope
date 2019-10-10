import unittest
import pickle

import numpy as np


from ms_deisotope import peak_set as module

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

        def test_missing_interval(self):
            x = np.arange(1000, 1200, 0.5)
            y = np.ones_like(x)

            peaks = [peak_cls(x[i], y[i], 1, 1, None, 0) for i in range(len(x))]
            ps = peak_set_cls(peaks)
            ps.reindex()

            for xi in np.arange(1000, 1200, 0.01):
                for p in ps.all_peaks_for(xi):
                    assert abs((xi - p.neutral_mass) / p.neutral_mass) < 1e-5

        if module.DeconvolutedPeakSet == peak_set_cls:

            def test_pickle(self):
                x = np.arange(1000, 1200, 0.5)
                y = np.ones_like(x)

                peaks = [peak_cls(x[i], y[i], 1, 1, None, 0) for i in range(len(x))]
                ps = peak_set_cls(peaks)
                ps.reindex()

                reloaded = pickle.loads(pickle.dumps(ps, -1))
                assert ps == reloaded
                assert ps._mz_ordered == reloaded._mz_ordered


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
