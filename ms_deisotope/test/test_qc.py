import unittest

from ms_deisotope.qc import xrea


def test_xrea():
    x = [1, 22.0, 4, 11., 12, 25, 16]
    y = xrea(x)
    assert abs(0.3104 - y) < 1e-3
