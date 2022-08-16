import unittest

import numpy as np

from ms_peak_picker import FittedPeak
from brainpy._c.isotopic_distribution import TheoreticalPeak as Peak

from ms_deisotope.scoring import (
        PenalizedMSDeconVFitter, MSDeconVFitter, DotProductFitter, ScaledGTestFitter,
        GTestFitter, LeastSquaresFitter)

experimental = [
    FittedPeak(mz=739.920, intensity=8356.829, signal_to_noise=100.000,
               peak_count=30, index=30588, full_width_at_half_max=0.050, area=425.445),
    FittedPeak(mz=740.255, intensity=8006.456, signal_to_noise=100.000,
               peak_count=31, index=30622, full_width_at_half_max=0.051, area=416.640),
    FittedPeak(mz=740.589, intensity=4970.605, signal_to_noise=100.000,
               peak_count=32, index=30655, full_width_at_half_max=0.050, area=253.103),
    FittedPeak(mz=740.923, intensity=2215.961, signal_to_noise=100.000,
               peak_count=33, index=30688, full_width_at_half_max=0.051, area=118.020)
]

theoretical = [
    Peak(mz=739.920306, intensity=8310.933747, charge=-3),
    Peak(mz=740.254733, intensity=8061.025466, charge=-3),
    Peak(mz=740.588994, intensity=4926.998052, charge=-3),
    Peak(mz=740.923235, intensity=2250.893651, charge=-3)
]


class IsotopicFitScoringTests(unittest.TestCase):

    def test_penalized_msdeconv(self):
        scorer = PenalizedMSDeconVFitter(20, 2.0)
        # score = scorer.evaluate(None, experimental, theoretical)
        scores = [scorer.evaluate(None, experimental, theoretical) for i in range(10)]
        score = scores[0]
        assert all([np.isclose(s, score) for s in scores[1:]]), scores
        self.assertAlmostEqual(score, 293.47483621051316, 3)
        score = scorer(None, experimental, theoretical)
        self.assertAlmostEqual(score, 293.47483621051316, 3)

    def test_msdeconv(self):
        scorer = MSDeconVFitter()
        score = scorer.evaluate(None, experimental, theoretical)
        self.assertAlmostEqual(score, 293.5135479544602, 3)
        score = scorer(None, experimental, theoretical)
        self.assertAlmostEqual(score, 293.5135479544602, 3)

    def test_dotproduct(self):
        scorer = DotProductFitter()
        score = scorer.evaluate(None, experimental, theoretical)
        self.assertAlmostEqual(score, 163471351.56044182, 3)
        score = scorer(None, experimental, theoretical)
        self.assertAlmostEqual(score, 163471351.56044182, 3)

    def test_scaled_g_test(self):
        scorer = ScaledGTestFitter()
        score = scorer.evaluate(None, experimental, theoretical)
        self.assertAlmostEqual(score, 6.59454124295588e-05, 3)
        score = scorer(None, experimental, theoretical)
        self.assertAlmostEqual(score, 6.59454124295588e-05, 3)

    def test_g_test(self):
        scorer = GTestFitter()
        score = scorer.evaluate(None, experimental, theoretical)
        self.assertAlmostEqual(score, 1.5531726368555354, 3)
        score = scorer(None, experimental, theoretical)
        self.assertAlmostEqual(score, 1.5531726368555354, 3)

    def test_least_squares(self):
        scorer = LeastSquaresFitter()
        score = scorer.evaluate(None, experimental, theoretical)
        self.assertAlmostEqual(score, 7.463484119741042e-05, 3)
        score = scorer(None, experimental, theoretical)
        self.assertAlmostEqual(score, 7.463484119741042e-05, 3)


if __name__ == '__main__':
    unittest.main()
