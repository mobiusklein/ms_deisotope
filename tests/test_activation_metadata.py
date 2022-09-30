import unittest

from ms_deisotope.data_source.metadata import activation


class TestActivationMetadata(unittest.TestCase):
    def test_is_a_relationships(self):
        assert activation.HCD.is_a(activation.CID)
        assert not activation.HCD.is_a(activation.ETD)

    def test_activation_information(self):
        ai = activation.ActivationInformation(activation.HCD, 30)
        assert ai.has_dissociation_type(activation.HCD)
        assert not ai.has_dissociation_type(activation.ETD)

        mai = activation.MultipleActivationInformation([activation.ETD, activation.HCD], [26, 24])
        assert mai.has_dissociation_type(activation.HCD)
        assert mai.has_dissociation_type(activation.ETD)


if __name__ == '__main__':
    unittest.main()
