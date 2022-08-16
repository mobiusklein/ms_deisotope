import unittest

from ms_deisotope.data_source.metadata import activation, cv

class TestTerm(unittest.TestCase):
    def test_cv_eq(self):
        assert activation.CID == activation.CID
        assert activation.CID == activation.ActivationInformation("CAD", 55.5).method
        assert activation.CID != activation.HCD
