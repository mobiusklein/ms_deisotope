import unittest

from ms_deisotope.averagine import (
    peptide, calculate_mass, average_compositions,
    _Averagine, Averagine, add_compositions,
    AveragineCache, _AveragineCache, TheoreticalIsotopicPattern,
    _TheoreticalIsotopicPattern, BasePeakToMonoisotopicOffsetEstimator)


tid1 = [
    (1000.00, 0.587117149327, 1),
    (1001.0029404, 0.314398842767, 1),
    (1002.005551, 0.0984840079067, 1)
]

tid2 = [
    (1000.0, 0.312953555329, 2),
    (1000.5014313042954, 0.339021408804, 2),
    (1001.002474283496, 0.213306639087, 2),
    (1001.503389167319, 0.0983516834107, 2),
    (1002.00427365722, 0.0363667133695, 2)
]

composition = peptide.base_composition


def make_averagine_suite(averagine_class):
    _peptide = averagine_class(composition)

    class TestAveragine(unittest.TestCase):
        def test_isotopic_cluster(self):
            tid = _peptide.isotopic_cluster(1000, 1)
            for i, peak in enumerate(tid):
                self.assertAlmostEqual(peak.mz, tid1[i][0], 3)
                self.assertAlmostEqual(peak.intensity, tid1[i][1], 3)
                self.assertEqual(peak.charge, 1)

            tid = _peptide.isotopic_cluster(1000, 2)
            for peak, match in zip(tid, tid2):
                self.assertAlmostEqual(peak.mz, match[0], 3)
                self.assertAlmostEqual(peak.intensity, match[1], 3)
                self.assertEqual(peak.charge, match[2])

        def test_truncate_after(self):
            tid = _peptide.isotopic_cluster(1000, 1, truncate_after=0.95)
            inst = _peptide.isotopic_cluster(1000, 1, truncate_after=1.0).truncate_after(0.95)

            for i, p in enumerate(inst):
                self.assertAlmostEqual(tid[i].mz, p.mz, 3)
                self.assertAlmostEqual(tid[i].intensity, p.intensity)

            self.assertAlmostEqual(tid.total(), 1.0, 3)
            self.assertAlmostEqual(inst.total(), 1.0, 3)

        def __repr__(self):
            r = super(TestAveragine, self).__repr__()
            return "%s(%s)" % (r, averagine_class)

    return TestAveragine


TestPurePythonAveragine = make_averagine_suite(_Averagine)
TestAveragine = make_averagine_suite(Averagine)
TestAveragineCache = make_averagine_suite(AveragineCache)
TestPurePythonAveragineCache = make_averagine_suite(_AveragineCache)


class TestSupportMethods(unittest.TestCase):
    def test_average_composition(self):
        avgd = average_compositions([composition, composition])
        for k, v in avgd.items():
            self.assertAlmostEqual(v, composition[k], 3)

    def test_add_composition(self):
        avgd = add_compositions({}, composition)
        for k, v in avgd.items():
            self.assertAlmostEqual(v, composition[k], 3)


class TestBasePeakToMonoisotopicOffsetEstimator(unittest.TestCase):
    def test_binned(self):
        est = BasePeakToMonoisotopicOffsetEstimator(peptide)
        assert est(1000) == 0
        assert est(2500) == 1
        assert est(3500) == 2

    def test_exact(self):
        est = BasePeakToMonoisotopicOffsetEstimator(peptide)
        assert est(3460.0, False) == 2


if __name__ == '__main__':
    unittest.main()
