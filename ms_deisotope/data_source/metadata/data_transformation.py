from collections import OrderedDict

from six import string_types as basestring

from .cv import Term, render_list


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


data_transformations = [
    DataTransformation(u'file format conversion', u'MS:1000530',
                       u'data transformation', [u'data transformation']),
    DataTransformation(u'data processing action', u'MS:1000543',
                       u'data transformation', [u'data transformation']),
    DataTransformation(u'Conversion to mzML', u'MS:1000544', u'data transformation', [
                       u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to mzMLb', u'MS:1002839', u'data transformation', [
                       u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to dta', u'MS:1000741', u'data transformation', [
                       u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to mzXML', u'MS:1000545', u'data transformation', [
                       u'file format conversion', u'data transformation']),
    DataTransformation(u'Conversion to mzData', u'MS:1000546', u'data transformation', [
                       u'file format conversion', u'data transformation']),
    DataTransformation(u'charge deconvolution', u'MS:1000034', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'peak picking', u'MS:1000035', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'deisotoping', u'MS:1000033', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'baseline reduction', u'MS:1000593', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'smoothing', u'MS:1000592', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'precursor recalculation', u'MS:1000780', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'data filtering', u'MS:1001486', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'retention time alignment', u'MS:1000745', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'intensity normalization', u'MS:1001484', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'm/z calibration', u'MS:1001485', u'data transformation',
                       [u'data processing action', u'data transformation']),
    DataTransformation(u'charge state calculation', u'MS:1000778', u'data transformation', [
                       u'data processing action', u'data transformation']),
    DataTransformation(u'area peak picking', u'MS:1000801', u'data transformation', [
                       u'peak picking', u'data processing action', u'data transformation']),
    DataTransformation(u'height peak picking', u'MS:1000802', u'data transformation', [
                       u'peak picking', u'data processing action', u'data transformation']),
    DataTransformation(u'top hat baseline reduction', u'MS:1001994', u'data transformation', [
                       u'baseline reduction', u'data processing action', u'data transformation']),
    DataTransformation(u'convex hull baseline reduction', u'MS:1001995', u'data transformation', [
                       u'baseline reduction', u'data processing action', u'data transformation']),
    DataTransformation(u'median baseline reduction', u'MS:1001996', u'data transformation', [
                       u'baseline reduction', u'data processing action', u'data transformation']),
    DataTransformation(u'wavelet transformation smoothing', u'MS:1001997', u'data transformation', [
                       u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'Gaussian smoothing', u'MS:1000784', u'data transformation', [
                       u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'moving average smoothing', u'MS:1000785', u'data transformation', [
                       u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'Savitzky-Golay smoothing', u'MS:1000782', u'data transformation',
                       [u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'LOWESS smoothing', u'MS:1000783', u'data transformation', [
                       u'smoothing', u'data processing action', u'data transformation']),
    DataTransformation(u'msPrefix precursor recalculation', u'MS:1000781', u'data transformation', [
                       u'precursor recalculation', u'data processing action', u'data transformation']),
    DataTransformation(u'low intensity data point removal', u'MS:1000594', u'data transformation', [
                       u'data filtering', u'data processing action', u'data transformation']),
    DataTransformation(u'high intensity data point removal', u'MS:1000746', u'data transformation', [
                       u'data filtering', u'data processing action', u'data transformation']),
    DataTransformation(u'area normalization', u'MS:1001999', u'data transformation', [
                       u'intensity normalization', u'data processing action', u'data transformation']),
    DataTransformation(u'below precursor intensity dominance charge state calculation', u'MS:1000779', u'data transformation', [
                       u'charge state calculation', u'data processing action', u'data transformation']),
    DataTransformation(u'sophisticated numerical annotation procedure', u'MS:1001998', u'data transformation', [
                       u'area peak picking', u'peak picking', u'data processing action', u'data transformation']),
]


data_transformations_by_name = {c.name: c for c in data_transformations}


def data_transformation(name):
    try:
        return data_transformations_by_name[name]
    except KeyError:
        return DataTransformation(name, name, name, [name])


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


class Software(object):
    def __init__(self, name, id=None, version=None, options=None):
        self.name = name
        self.id = id or name
        self.version = version
        self.options = options

    def __str__(self):
        return self.name

    def __repr__(self):
        template = '{self.__class__.__name__}({self.name}, {self.id}, {self.version})'
        return template.format(self=self)

    def __eq__(self, other):
        try:
            return self.name == other.name
        except AttributeError:
            return str(self) == str(other)


if __name__ == '__main__':
    __generate_list_code()
