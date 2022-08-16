'''A simple data holder describing a sample which was used to generate
a mass spectrometry dataset.
'''
from ms_deisotope.utils import Base, dict_proxy


@dict_proxy("parameters")
class Sample(Base):
    """Describes a sample used to generate a dataset

    Attributes
    ----------
    id : :class:`object`
        A unique identifier across the samples with which to reference this sample
    name : :class:`str`
        A name for the sample
    parameters : :class:`dict`
        A collection of descriptive

    """
    def __init__(self, id, name=None, parameters=None, **kwargs):
        if parameters is None:
            parameters = dict(kwargs)
        else:
            parameters = dict(parameters)
            parameters.update(kwargs)
        self.id = id
        self.name = name
        self.parameters = parameters
