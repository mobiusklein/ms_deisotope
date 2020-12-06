from __future__ import print_function
import operator
import random

from datetime import datetime
from collections import OrderedDict

from six import add_metaclass


try:  # pragma: no cover
    from .plot import (
        draw_peaklist, draw_raw, annotate_scan, annotate_scan_single,
        has_plot)
except (RuntimeError, ImportError):  # pragma: no cover
    has_plot = False

    def annotate_scan(scan, products, nperrow=4, ax=None):
        raise ImportError('matplotlib')

    def annotate_scan_single(*args, **kwargs):
        raise ImportError('matplotlib')

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


def simple_repr(self):  # pragma: no cover
    '''A convenient function for automatically generating a ``__repr__``-like
    string for arbitrary objects.

    Returns
    -------
    str
    '''
    template = "{self.__class__.__name__}({d})"

    def formatvalue(v):
        if isinstance(v, float):
            return "%0.4f" % v
        else:
            return str(v)

    if not hasattr(self, "__slots__") or len(self.__slots__) == 0 or hasattr(self, '__dict__'):
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
    '''A convenience base class for non-critical code to provide types
    with automatic :meth:`__repr__` methods using :func:`simple_repr`
    '''
    __slots__ = ()
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

    def __contains__(self, i):
        return i in self.name

    def __getitem__(self, i):
        return self.name[i]

    def __iter__(self):
        return iter(self.name)

    def __len__(self):
        return len(self.name)


class DeconvolutionProcessResult(object):
    """Hold information from a multi-stage deconvolution process.

    Emulates an interface matching a tuple of (:attr:`peak_set`, :attr:`priorities`)

    Attributes
    ----------
    deconvoluter : :class:`~.DeconvoluterBase`
        The deconvoluter used.
    errors : list, optional
        A list of :class:`~.Exception` instances that may have been
        thrown during any stage of deconvolution.
    peak_set : :class:`~.DeconvolutedPeakSet`
        The resulting deconvoluted peaks
    priorities : list
        The extracted results from the targeted deconvolution list
    """

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
    deconvoluter : :class:`~.DeconvoluterBase`
        The deconvoluter to use to look up the result
    """
    def __init__(self, deconvoluter, *args, **kwargs):
        self.deconvoluter = deconvoluter

    def get(self):
        """Fetch the optimal solution after the computation has finished.

        Returns
        -------
        :class:`~.DeconvolutedPeak`
        """
        raise NotImplementedError()


class TrivialTargetedDeconvolutionResult(TargetedDeconvolutionResultBase):
    """Stores the necessary information to retrieve the local optimal solution for
    a single peak for a deconvolution algorithm where the local optimum is the best
    solution.

    Attributes
    ----------
    query_peak : :class:`~.FittedPeak`
        The peak queried with
    solution_peak : :class:`~.DeconvolutedPeak`
        The optimal solution peak
    """
    def __init__(self, deconvoluter, solution_peak, query_peak, *args, **kwargs):
        super(TrivialTargetedDeconvolutionResult, self).__init__(deconvoluter, *args, **kwargs)
        self.solution_peak = solution_peak
        self.query_peak = query_peak

    def get(self):
        """Fetch the optimal solution after the computation has finished.

        Returns
        -------
        :class:`~.DeconvolutedPeak`
        """
        return self.solution_peak


def ppm_error(x, y):
    return (x - y) / y


def ppm2da(mass, error_tolerance):
    x = mass + mass * error_tolerance
    return x - mass


def da2ppm(mass, error_tolerance):
    x = mass + error_tolerance
    return (x - mass) / mass


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
    '''Generate a random "universally unique" ID number with ``n``
    bits of entropy.

    Parameters
    ----------
    n: int
        The number of random bits to generate

    Returns
    -------
    int
    '''
    int_ = random.getrandbits(n)
    return int_


class LRUDict(object):  # pragma: no cover
    def __init__(self, *args, **kwargs):
        maxsize = kwargs.pop("maxsize", 24)
        self.store = OrderedDict()
        self.maxsize = maxsize
        self.purge()

    def __len__(self):
        return len(self.store)

    def popitem(self, last=True):
        return self.store.popitem(last=last)

    def pop(self, key, default=None):
        return self.store.pop(key, default)

    def purge(self):
        overflow = max(0, len(self) - self.maxsize)
        for _ in range(overflow):
            self.popitem(last=False)

    def __repr__(self):
        return "LRUDict(%r)" % (dict(self.store),)

    def __contains__(self, key):
        return key in self.store

    def __iter__(self):
        return iter(self.store)

    def keys(self):
        return self.store.keys()

    def values(self):
        return self.store.values()

    def items(self):
        return self.store.items()

    def __getitem__(self, key):
        value = self.store[key]
        self._mark_used(key)
        return value

    def __setitem__(self, key, value):
        self.store[key] = value
        self.purge()

    def _mark_used(self, key):
        value = self.store.pop(key, None)
        self.store[key] = value


def decimal_shift(x):
    i = 1.0
    while i < 10 ** 10:
        if round(x * i) > 0:
            return 1.0 / i
        i *= 10.
    return 1.0 / i


class _MappingOverAttributeProxy(object):
    '''A replacement for __dict__ for unpickling an object which once
    has __slots__ now but did not before.'''

    def __init__(self, obj):
        self.obj = obj

    def __getitem__(self, key):
        return getattr(self.obj, key)

    def __setitem__(self, key, value):
        setattr(self.obj, key, value)

    def __contains__(self, key):
        return hasattr(self.obj, key)

    def __repr__(self):
        return "{self.__class__.__name__}({self.obj})".format(self=self)
