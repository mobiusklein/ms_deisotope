import unittest

from ms_deisotope.peak_dependency_network import intervals
from ms_deisotope.peak_dependency_network.intervals import Interval, IntervalTreeNode


iv0 = Interval(0, 10)
iv1 = Interval(-10, -5)
iv2 = Interval(-10, 0)
iv3 = Interval(-10, 5)
iv4 = Interval(-10, 10)
iv5 = Interval(-10, 20)
iv6 = Interval(0, 20)
iv7 = Interval(5, 20)
iv8 = Interval(10, 20)
iv9 = Interval(15, 20)
iv10 = Interval(-5, 0)


class TestIntervals(unittest.TestCase):
    def test_interval_overlaps_interval(self):
        assert iv0.overlaps(iv0)
        assert not iv0.overlaps(iv1)
        assert iv0.overlaps(iv2)
        assert iv0.overlaps(iv3)
        assert iv0.overlaps(iv4)
        assert iv0.overlaps(iv5)
        assert iv0.overlaps(iv6)
        assert iv0.overlaps(iv7)
        assert iv0.overlaps(iv8)
        assert not iv0.overlaps(iv9)

    def test_contains_interval(self):
        assert iv0.contains_interval(iv0)
        assert not iv0.contains_interval(iv1)
        assert not iv0.contains_interval(iv2)
        assert not iv0.contains_interval(iv3)
        assert not iv0.contains_interval(iv4)
        assert not iv0.contains_interval(iv5)
        assert not iv0.contains_interval(iv6)
        assert not iv0.contains_interval(iv7)
        assert not iv0.contains_interval(iv8)
        assert not iv0.contains_interval(iv9)
        assert not iv0.contains_interval(iv10)

        assert not iv2.contains_interval(iv0)
        assert iv2.contains_interval(iv1)
        assert iv2.contains_interval(iv2)
        assert not iv2.contains_interval(iv3)
        assert not iv2.contains_interval(iv4)
        assert not iv2.contains_interval(iv5)
        assert not iv2.contains_interval(iv6)
        assert not iv2.contains_interval(iv7)
        assert not iv2.contains_interval(iv8)
        assert not iv2.contains_interval(iv9)
        assert iv2.contains_interval(iv10)


if __name__ == '__main__':
    unittest.main()
