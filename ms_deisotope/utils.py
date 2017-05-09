from __future__ import print_function
import operator
import random
from datetime import datetime

try:
    from matplotlib import pyplot as plt
    from ms_peak_picker.utils import draw_peaklist, draw_raw
    has_plot = True

except ImportError:
    has_plot = False

try:
    range = xrange
except NameError:
    range = range


try:  # pragma: no cover
    i128 = long
    basestring = basestring
except NameError:  # pragma: no cover
    i128 = int
    basestring = (bytes, str)


def printer(message):
    print(datetime.now().isoformat(' ') + ' ' + str(message))


def debug_printer(message):
    print("DEBUG:" + datetime.now().isoformat(' ') + ' ' + str(message))


# From six
def add_metaclass(metaclass):
    """Class decorator for creating a class with a metaclass."""
    def wrapper(cls):
        orig_vars = cls.__dict__.copy()
        slots = orig_vars.get('__slots__')
        if slots is not None:
            if isinstance(slots, str):
                slots = [slots]
            for slots_var in slots:
                orig_vars.pop(slots_var)
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper


def simple_repr(self):  # pragma: no cover
    template = "{self.__class__.__name__}({d})"

    def formatvalue(v):
        if isinstance(v, float):
            return "%0.4f" % v
        else:
            return str(v)

    if not hasattr(self, "__slots__"):
        d = [
            "%s=%s" % (k, formatvalue(v)) if v is not self else "(...)" for k, v in sorted(
                self.__dict__.items(), key=lambda x: x[0])
            if (not k.startswith("_") and not callable(v)) and not (v is None)]
    else:
        d = [
            "%s=%s" % (k, formatvalue(v)) if v is not self else "(...)" for k, v in sorted(
                [(name, getattr(self, name)) for name in self.__slots__], key=lambda x: x[0])
            if (not k.startswith("_") and not callable(v)) and not (v is None)]

    return template.format(self=self, d=', '.join(d))


class Base(object):
    __repr__ = simple_repr


class Constant(object):
    """A convenience type meant to be used to instantiate singletons for signaling
    specific states in return values.

    Attributes
    ----------
    name: str
        The name of the constant
    """
    def __init__(self, name, is_true=True):
        self.name = name
        self.is_true = is_true

    def __eq__(self, other):
        return self.name == str(other)

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return str(self.name)

    def __nonzero__(self):
        return self.is_true

    def __bool__(self):
        return self.is_true


class DeconvolutionProcessResult(object):
    def __init__(self, deconvoluter, peak_set, priorities, errors=None):
        self.deconvoluter = deconvoluter
        self.peak_set = peak_set
        self.priorities = priorities
        self.errors = errors

    def __getitem__(self, i):
        if i == 0:
            return self.peak_set
        elif i == 1:
            return self.priorities
        else:
            raise IndexError(i)

    def __iter__(self):
        yield self.peak_set
        yield self.priorities

    def __repr__(self):
        return "DeconvolutionProcessResult(%s, %s)" % tuple(self)


class TargetedDeconvolutionResultBase(Base):
    """Base class to store all of the necessary information to retrieve the optimal
    solution for a single peak.

    Attributes
    ----------
    deconvoluter : DeconvoluterBase
        The deconvoluter to use to look up the result
    """
    def __init__(self, deconvoluter, *args, **kwargs):
        self.deconvoluter = deconvoluter

    def get(self):
        """Fetch the optimal solution after the computation has finished.

        Returns
        -------
        DeconvolutedPeak
        """
        raise NotImplementedError()


class TrivialTargetedDeconvolutionResult(TargetedDeconvolutionResultBase):
    """Stores the necessary information to retrieve the local optimal solution for
    a single peak for a deconvolution algorithm where the local optimum is the best
    solution.

    Attributes
    ----------
    query_peak : FittedPeak
        The peak queried with
    solution_peak : DeconvolutedPeak
        The optimal solution peak
    """
    def __init__(self, deconvoluter, solution_peak, query_peak, *args, **kwargs):
        super(TrivialTargetedDeconvolutionResult, self).__init__(deconvoluter, *args, **kwargs)
        self.solution_peak = solution_peak
        self.query_peak = query_peak

    def get(self):
        return self.solution_peak


def ppm_error(x, y):
    return (x - y) / y


def dict_proxy(attribute):
    """Return a decorator for a class to give it a `dict`-like API proxied
    from one of its attributes

    Parameters
    ----------
    attribute : str
        The string corresponding to the attribute which will
        be used

    Returns
    -------
    function
    """
    getter = operator.attrgetter(attribute)

    def wrap(cls):

        def __getitem__(self, key):
            return getter(self)[key]

        def __setitem__(self, key, value):
            d = getter(self)
            d[key] = value

        def keys(self):
            return getter(self).keys()

        def __iter__(self):
            return iter(getter(self))

        def items(self):
            return getter(self).items()

        cls.__getitem__ = __getitem__
        cls.__setitem__ = __setitem__
        cls.keys = keys
        cls.items = items
        cls.__iter__ = __iter__

        return cls
    return wrap


def uid(n=128):
    int_ = random.getrandbits(n)
    return int_
