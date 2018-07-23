from collections import OrderedDict

from six import string_types as basestring

from .cv import Term, render_list
from .software import Software, SoftwareName


class DataTransformation(Term):
    """Describes a named data transformation, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`DataTransformation` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


def __generate_list_code():
    '''Prints the code to generate these static lists
    '''
    render_list('data transformation', term_cls_name="DataTransformation")


data_transformations = []


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('data transformation', term_cls_name="DataTransformation", writer=cog.out)
# ]]]
data_transformations = [
    DataTransformation(u'file format conversion', u'MS:1000530', u'Conversion of one file format to another.',
                       'data transformation', [u'data transformation']),
    DataTransformation(u'data processing action', u'MS:1000543', u'Data processing attribute used to describe the type of data processing performed on the data file.',
                       'data transformation', [u'data transformation']),
    DataTransformation(u'Conversion to mzML', u'MS:1000544', u'Conversion of a file format to Proteomics Standards Initiative mzML file format.',
                       'data transformation', [u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to mzMLb', u'MS:1002839', u'Conversion of a file format to Proteomics Standards Initiative mzMLb file format.',
                       'data transformation', [u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to dta', u'MS:1000741', u'Conversion to dta format.',
                       'data transformation', [u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to mzXML', u'MS:1000545', u'Conversion of a file format to Institute of Systems Biology mzXML file format.',
                       'data transformation', [u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to mzData', u'MS:1000546', u'Conversion of a file format to Proteomics Standards Initiative mzData file format.',
                       'data transformation', [u'file format conversion', u'data transformation']),
    DataTransformation(u'charge deconvolution', u'MS:1000034', u'The determination of the mass of an ion based on the mass spectral peaks that represent multiple-charge ions.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'peak picking', u'MS:1000035', u'Spectral peak processing conducted on the acquired data to convert profile data to centroided data.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'deisotoping', u'MS:1000033', u'The removal of isotope peaks to represent the fragment ion as one data point and is commonly done to reduce complexity. It is done in conjunction with the charge state deconvolution.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'baseline reduction', u'MS:1000593', u'A process of removal of varying intensities generated due to variable energy absorption before further processing can take place. Baseline reduction facilitates meaningful comparision between intensities of m/z values.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'smoothing', u'MS:1000592', u'A process of reducing spikes of intensity in order to reduce noise while preserving real peak signal. Many algorithms can be applied for this process.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'precursor recalculation', u'MS:1000780', u'A process that recalculates existing precursor selected ions with one or more algorithmically determined precursor selected ions.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'data filtering', u'MS:1001486', u'Filtering out part of the data.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'retention time alignment', u'MS:1000745', u'The correction of the spectrum scan times, as used e.g. in label-free proteomics.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'intensity normalization', u'MS:1001484', u'Normalization of data point intensities.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'm/z calibration', u'MS:1001485', u'Calibration of data point m/z positions.',
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'charge state calculation', u'MS:1000778', u"A process that infers the charge state of an MSn spectrum's precursor(s) by the application of some algorithm.",
                       'data transformation', [u'data processing action', u'data transformation']),
    DataTransformation(u'area peak picking', u'MS:1000801', u'Spectral peak processing conducted on the acquired data to convert profile data to centroided data. The area defined by all raw data points that belong to the peak is reported.',
                       'data transformation', [u'peak picking', u'data processing action', u'data transformation']),
    DataTransformation(u'height peak picking', u'MS:1000802', u'Spectral peak processing conducted on the acquired data to convert profile data to centroided data. The maximum intensity of all raw data points that belong to the peak is reported.',
                       'data transformation', [u'peak picking', u'data processing action', u'data transformation']),
    DataTransformation(u'top hat baseline reduction', u'MS:1001994', u"Top-hat morphological filter based on the basic morphological operations 'erosion' and 'dilatation'.",
                       'data transformation', [u'baseline reduction', u'data processing action', u'data transformation']),
    DataTransformation(u'convex hull baseline reduction', u'MS:1001995', u'Constructs the baseline by fitting multiple parabolas to the spectrum starting with the large scale structures.',
                       'data transformation', [u'baseline reduction', u'data processing action', u'data transformation']),
    DataTransformation(u'median baseline reduction', u'MS:1001996', u'The spectrum that will be baseline subtracted is divided into a number of segments.',
                       'data transformation', [u'baseline reduction', u'data processing action', u'data transformation']),
    DataTransformation(u'wavelet transformation smoothing', u'MS:1001997', u'The random noise is removed by using the undecimated wavelet transform." [DOI:10.1093/bioinformatics/btl355',       'data transformation', [
                       u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'Gaussian smoothing', u'MS:1000784', u'Reduces intensity spikes by convolving the data with a one-dimensional Gaussian function.',
                       'data transformation', [u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'moving average smoothing', u'MS:1000785', u'Reduces intensity spikes by averaging each point with two or more adjacent points. The more adjacent points that used, the stronger the smoothing effect.',
                       'data transformation', [u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'Savitzky-Golay smoothing', u'MS:1000782', u"Reduces intensity spikes by applying local polynomial regression (of degree k) on a distribution (of at least k+1 equally spaced points) to determine the smoothed value for each point. It tends to preserve features of the distribution such as relative maxima, minima and width, which are usually 'flattened' by other adjacent averaging techniques.",
                       'data transformation', [u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'LOWESS smoothing', u'MS:1000783', u'Reduces intensity spikes by applying a modelling method known as locally weighted polynomial regression. At each point in the data set a low-degree polynomial is fit to a subset of the data, with explanatory variable values near the point whose response is being estimated. The polynomial is fit using weighted least squares, giving more weight to points near the point whose response is being estimated and less weight to points further away. The value of the regression function for the point is then obtained by evaluating the local polynomial using the explanatory variable values for that data point. The LOESS fit is complete after regression function values have been computed for each of the n data points. Many of the details of this method, such as the degree of the polynomial model and the weights, are flexible.',
                       'data transformation', [u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'msPrefix precursor recalculation', u'MS:1000781', u'Recalculates one or more precursor selected ions by peak detection in the isolation windows of high accuracy MS precursor scans.',
                       'data transformation', [u'precursor recalculation', u'data processing action', u'data transformation']),
    DataTransformation(u'low intensity data point removal', u'MS:1000594', u'The removal of very low intensity data points that are likely to be spurious noise rather than real signal.',
                       'data transformation', [u'data filtering', u'data processing action', u'data transformation']),
    DataTransformation(u'high intensity data point removal', u'MS:1000746', u'The removal of very high intensity data points.',
                       'data transformation', [u'data filtering', u'data processing action', u'data transformation']),
    DataTransformation(u'area normalization', u'MS:1001999', u'Normalization of areas below the curves.',
                       'data transformation', [u'intensity normalization', u'data processing action', u'data transformation']),
    DataTransformation(u'below precursor intensity dominance charge state calculation', u'MS:1000779', u'Infers charge state as single or ambiguously multiple by determining the fraction of intensity below the precursor m/z.',
                       'data transformation', [u'charge state calculation', u'data processing action', u'data transformation']),
    DataTransformation(u'sophisticated numerical annotation procedure', u'MS:1001998', u'It searches for known patterns in the measured spectrum." [DOI:10.1021/ac951158i',       'data transformation', [
                       u'area peak picking', u'peak picking', u'data processing action', u'data transformation']),
]
# [[[end]]]


data_transformations_by_name = {c.name: c for c in data_transformations}


def data_transformation(name):
    try:
        return data_transformations_by_name[name]
    except KeyError:
        return DataTransformation(name, name, name, name, [name])


class ProcessingMethod(object):

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
        if isinstance(operation, basestring):
            operation = data_transformation(operation)
        self.operations[operation] = value

    def update(self, operations):
        if isinstance(operations, dict):
            for op, val in operations.items():
                self.add(op, val)
        else:
            for op in operations:
                if isinstance(op, dict):
                    self.add(op['name'], op['value'])
                else:
                    self.add(op)


class DataProcessingInformation(object):

    def __init__(self, methods=None, id=None):
        self.methods = methods or []
        self.id = id

    def __getitem__(self, i):
        return self.methods[i]

    def __len__(self):
        return len(self.methods)

    def __iter__(self):
        return iter(self.methods)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.methods}, id={self.id!r})"
        return template.format(self=self)


if __name__ == '__main__':
    __generate_list_code()
