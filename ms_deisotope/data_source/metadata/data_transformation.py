"""Defines types for describing mass spectrometry data transformations,
data processing sequences, and a database of controlled vocabulary terms
for them.

.. note::
    The entities in this module are only for *naming* these transformations
    and how they were combined to process mass spectra. For implementations,
    see other features of :mod:`ms_deisotope` and :mod:`ms_peak_picker`, or
    other libraries.
"""

from collections import OrderedDict

from six import string_types as basestring

from ms_deisotope.utils import uid

from .cv import Term, TermSet
# from .software import Software, SoftwareName


class DataTransformation(Term):
    """Describes a named data transformation, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`DataTransformation` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


data_transformations = []


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('data transformation', term_cls_name="DataTransformation", writer=cog.out)
# ]]]
# CV Version: 4.1.55
data_transformations = TermSet([
    DataTransformation('file format conversion', 'MS:1000530',
                       ('Conversion of one file format to another.'),
                       'data transformation',
                       ['data transformation']),
    DataTransformation('data processing action', 'MS:1000543',
                       ('Data processing attribute used to describe the type of data '
                        'processing performed on the data file.'),
                       'data transformation',
                       ['data transformation']),
    DataTransformation('Conversion to mzML', 'MS:1000544',
                       ('Conversion of a file format to Proteomics Standards '
                        'Initiative mzML file format.'),
                       'data transformation',
                       ['file format conversion', 'data transformation']),
    DataTransformation('Conversion to mzXML', 'MS:1000545',
                       ('Conversion of a file format to Institute of Systems Biology '
                        'mzXML file format.'),
                       'data transformation',
                       ['file format conversion', 'data transformation']),
    DataTransformation('Conversion to mzData', 'MS:1000546',
                       ('Conversion of a file format to Proteomics Standards '
                        'Initiative mzData file format.'),
                       'data transformation',
                       ['file format conversion', 'data transformation']),
    DataTransformation('Conversion to dta', 'MS:1000741',
                       ('Conversion to dta format.'),
                       'data transformation',
                       ['file format conversion', 'data transformation']),
    DataTransformation('Conversion to mzMLb', 'MS:1002839',
                       ('Conversion of a file format to Proteomics Standards '
                        'Initiative mzMLb file format.'),
                       'data transformation',
                       ['file format conversion', 'data transformation']),
    DataTransformation('deisotoping', 'MS:1000033',
                       ('The removal of isotope peaks to represent the fragment ion '
                        'as one data point and is commonly done to reduce complexity. '
                        'It is done in conjunction with the charge state '
                        'deconvolution.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('charge deconvolution', 'MS:1000034',
                       ('The determination of the mass of an ion based on the mass '
                        'spectral peaks that represent multiple-charge ions.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('peak picking', 'MS:1000035',
                       ('Spectral peak processing conducted on the acquired data to '
                        'convert profile data to centroided data.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('smoothing', 'MS:1000592',
                       ('A process of reducing spikes of intensity in order to reduce '
                        'noise while preserving real peak signal. Many algorithms can '
                        'be applied for this process.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('baseline reduction', 'MS:1000593',
                       ('A process of removal of varying intensities generated due to '
                        'variable energy absorption before further processing can '
                        'take place. Baseline reduction facilitates meaningful '
                        'comparision between intensities of m/z values.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('retention time alignment', 'MS:1000745',
                       ('The correction of the spectrum scan times, as used e.g. in '
                        'label-free proteomics.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('charge state calculation', 'MS:1000778',
                       ("A process that infers the charge state of an MSn spectrum's "
                        'precursor(s) by the application of some algorithm.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('precursor recalculation', 'MS:1000780',
                       ('A process that recalculates existing precursor selected ions '
                        'with one or more algorithmically determined precursor '
                        'selected ions.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('intensity normalization', 'MS:1001484',
                       ('Normalization of data point intensities.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('m/z calibration', 'MS:1001485',
                       ('Calibration of data point m/z positions.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('data filtering', 'MS:1001486',
                       ('Filtering out part of the data.'),
                       'data transformation',
                       ['data processing action', 'data transformation']),
    DataTransformation('area peak picking', 'MS:1000801',
                       ('Spectral peak processing conducted on the acquired data to '
                        'convert profile data to centroided data. The area defined by '
                        'all raw data points that belong to the peak is reported.'),
                       'data transformation',
                       ['peak picking', 'data processing action', 'data transformation']),
    DataTransformation('height peak picking', 'MS:1000802',
                       ('Spectral peak processing conducted on the acquired data to '
                        'convert profile data to centroided data. The maximum '
                        'intensity of all raw data points that belong to the peak is '
                        'reported.'),
                       'data transformation',
                       ['peak picking', 'data processing action', 'data transformation']),
    DataTransformation('Savitzky-Golay smoothing', 'MS:1000782',
                       ('Reduces intensity spikes by applying local polynomial '
                        'regression (of degree k) on a distribution (of at least k+1 '
                        'equally spaced points) to determine the smoothed value for '
                        'each point. It tends to preserve features of the '
                        'distribution such as relative maxima, minima and width, '
                        "which are usually 'flattened' by other adjacent averaging "
                        'techniques.'),
                       'data transformation',
                       ['smoothing', 'data processing action', 'data transformation']),
    DataTransformation('LOWESS smoothing', 'MS:1000783',
                       ('Reduces intensity spikes by applying a modelling method '
                        'known as locally weighted polynomial regression. At each '
                        'point in the data set a low-degree polynomial is fit to a '
                        'subset of the data, with explanatory variable values near '
                        'the point whose response is being estimated. The polynomial '
                        'is fit using weighted least squares, giving more weight to '
                        'points near the point whose response is being estimated and '
                        'less weight to points further away. The value of the '
                        'regression function for the point is then obtained by '
                        'evaluating the local polynomial using the explanatory '
                        'variable values for that data point. The LOESS fit is '
                        'complete after regression function values have been computed '
                        'for each of the n data points. Many of the details of this '
                        'method, such as the degree of the polynomial model and the '
                        'weights, are flexible.'),
                       'data transformation',
                       ['smoothing', 'data processing action', 'data transformation']),
    DataTransformation('Gaussian smoothing', 'MS:1000784',
                       ('Reduces intensity spikes by convolving the data with a one- '
                        'dimensional Gaussian function.'),
                       'data transformation',
                       ['smoothing', 'data processing action', 'data transformation']),
    DataTransformation('moving average smoothing', 'MS:1000785',
                       ('Reduces intensity spikes by averaging each point with two or '
                        'more adjacent points. The more adjacent points that used, '
                        'the stronger the smoothing effect.'),
                       'data transformation',
                       ['smoothing', 'data processing action', 'data transformation']),
    DataTransformation('wavelet transformation smoothing', 'MS:1001997',
                       ('The random noise is removed by using the undecimated wavelet '
                        'transform." [DOI:10.1093/bioinformatics/btl355'),
                       'data transformation',
                       ['smoothing', 'data processing action', 'data transformation']),
    DataTransformation('top hat baseline reduction', 'MS:1001994',
                       ('Top-hat morphological filter based on the basic '
                        "morphological operations 'erosion' and 'dilatation'."),
                       'data transformation',
                       ['baseline reduction', 'data processing action', 'data transformation']),
    DataTransformation('convex hull baseline reduction', 'MS:1001995',
                       ('Constructs the baseline by fitting multiple parabolas to the '
                        'spectrum starting with the large scale structures.'),
                       'data transformation',
                       ['baseline reduction', 'data processing action', 'data transformation']),
    DataTransformation('median baseline reduction', 'MS:1001996',
                       ('The spectrum that will be baseline subtracted is divided '
                        'into a number of segments.'),
                       'data transformation',
                       ['baseline reduction', 'data processing action', 'data transformation']),
    DataTransformation('below precursor intensity dominance charge state calculation', 'MS:1000779',
                       ('Infers charge state as single or ambiguously multiple by '
                        'determining the fraction of intensity below the precursor '
                        'm/z.'),
                       'data transformation',
                       ['charge state calculation', 'data processing action', 'data transformation']),
    DataTransformation('msPrefix precursor recalculation', 'MS:1000781',
                       ('Recalculates one or more precursor selected ions by peak '
                        'detection in the isolation windows of high accuracy MS '
                        'precursor scans.'),
                       'data transformation',
                       ['precursor recalculation', 'data processing action', 'data transformation']),
    DataTransformation('area normalization', 'MS:1001999',
                       ('Normalization of areas below the curves.'),
                       'data transformation',
                       ['intensity normalization', 'data processing action', 'data transformation']),
    DataTransformation('low intensity data point removal', 'MS:1000594',
                       ('The removal of very low intensity data points that are '
                        'likely to be spurious noise rather than real signal.'),
                       'data transformation',
                       ['data filtering', 'data processing action', 'data transformation']),
    DataTransformation('high intensity data point removal', 'MS:1000746',
                       ('The removal of very high intensity data points.'),
                       'data transformation',
                       ['data filtering', 'data processing action', 'data transformation']),
    DataTransformation('sophisticated numerical annotation procedure', 'MS:1001998',
                       ('It searches for known patterns in the measured spectrum." '
                        '[DOI:10.1021/ac951158i'),
                       'data transformation',
                       ['area peak picking', 'peak picking', 'data processing action', 'data transformation']),
])
# [[[end]]]


def data_transformation(name):
    '''Translate a given name or identifier into a :class:`DataTransformation`
    instance.

    If no match is found in the database of known :class:`DataTransformation`
    types, a new dummy :class:`DataTransformation` is returned with all fields
    set to the value of ``name``

    Returns
    -------
    DataTransformation
    '''
    try:
        return data_transformations[name]
    except KeyError:
        return DataTransformation(name, name, name, name, [name])


class ProcessingMethod(object):
    '''Describes a single action applied to the associated spectra, composed of one
    or more data transformation operations.

    Attributes
    ----------
    operations: :class:`~.OrderedDict`
        The data transformations applied, mapping transformation term to any
        parameter value or the empty string.
    software_id: str
        A reference to a particular version of software used to carry out the
        transformations.
    order: int
        The order in which this action was performed, in the context of the
        parent :class:`DataProcessingInformation` object

    See Also
    --------
    :class:`DataTransformation`
    '''

    def __init__(self, operations=None, software_id=None, order=0):
        self.operations = operations or OrderedDict()
        self.software_id = software_id
        self.order = order

    def __getitem__(self, i):
        return list(self.operations.items())[i]

    def __len__(self):
        return len(self.operations)

    def __iter__(self):
        return iter(self.operations.items())

    def _get_names(self):
        return ["%r: %r" % (t.name, v) for t, v in self.operations.items()]

    def __repr__(self):
        template = "{self.__class__.__name__}({names}, order={self.order})"
        return template.format(self=self, names='[%s]' % ', '.join(self._get_names()))

    def add(self, operation, value=''):
        '''Add a new :class:`DataTransformation` operation to this method.

        Parameters
        ----------
        operation: :class:`DataTransformation` or :class:`str`
            The transformation applied. If passed as a :class:`str`, it will be translated
            into a :class:`DataTransformation`.
        value: :class:`object`
            A parameter value associated with the operation. Defaults to the empty string.
        '''
        if isinstance(operation, basestring):
            operation = data_transformation(operation)
        self.operations[operation] = value

    def update(self, operations):
        '''Add each operation in ``operations`` to this method.

        Parameters
        ----------
        operations: :class:`~.Mapping` or :class:`Iterable`
            If given a :class:`~.Mapping`, add each key-value
            pair using :meth:`add`. Otherwise extract keys or key-value pairs
            if :class:`Iterable` is supported, adding each in succession
        '''
        if isinstance(operations, dict):
            for op, val in operations.items():
                self.add(op, val)
        else:
            for op in operations:
                if isinstance(op, dict):
                    self.add(op['name'], op['value'])
                elif isinstance(op, (tuple, list)) and len(op) == 2:
                    self.add(*op)
                else:
                    self.add(op)


class DataProcessingInformation(object):
    '''A :class:`Sequence` of :class:`ProcessingMethod` instances, describing the sequence
    of data transformations applied to the spectra in the data associated with this instance.

    Attributes
    ----------
    methods: list
        A list of :class:`ProcessingMethod` instances describing the individual operations done
        in an ordered fashion.
    id: int
        A within-source unique identifier
    '''

    def __init__(self, methods=None, id=None):
        self.methods = methods or []
        self.id = id or uid()

    def __getitem__(self, i):
        return self.methods[i]

    def __len__(self):
        return len(self.methods)

    def __iter__(self):
        return iter(self.methods)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.methods}, id={self.id!r})"
        return template.format(self=self)


__all__ = [
    "DataTransformation", "data_transformations", 'data_transformation',
    "ProcessingMethod", "DataProcessingInformation"
]
