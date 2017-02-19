import unittest

from ms_deisotope.averagine import (
    peptide, calculate_mass, average_compositions,
    _Averagine, Averagine, add_compositions,
    AveragineCache, _AveragineCache)


tid1 = [
    (1000.00, 0.571474, 1),
    (1001.0029404, 0.306022, 1),
    (1002.005551, 0.09586, 1)
]

tid2 = [
    (1000.0, 0.3082404601338808, 2),
    (1000.5014313042954, 0.3339157305154826, 2),
    (1001.002474283496, 0.21009423111625725, 2),
    (1001.503389167319, 0.09687050245407176, 2),
    (1002.00427365722, 0.0358190289635575, 2)
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


if __name__ == '__main__':
    unittest.main()
