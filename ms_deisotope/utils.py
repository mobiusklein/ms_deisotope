from __future__ import print_function
import operator
import random
import math
from datetime import datetime
from collections import OrderedDict
from six import add_metaclass


try:
    from ms_peak_picker.utils import draw_peaklist, draw_raw
    from matplotlib import pyplot as plt, gridspec
    has_plot = True

    def annotate_scan(scan, products, nperrow=4, ax=None):
        if ax is None:
            figure = plt.figure()
        else:
            figure = ax.figure
        n = len(products)
        if n > 0:
            gs = gridspec.GridSpec(1 + max(int(math.ceil(n / float(nperrow))), 1), nperrow)
        else:
            gs = gridspec.GridSpec(1 + (n / nperrow), nperrow)
        ax = figure.add_subplot(gs[0, :])
        if scan.is_profile:
            draw_raw(scan.arrays, ax=ax)
        if scan.peak_set is None:
            scan.pick_peaks()
        draw_peaklist(scan.peak_set, ax=ax, lw=0.5, alpha=0.75)
        if scan.deconvoluted_peak_set:
            draw_peaklist(scan.deconvoluted_peak_set, ax=ax, lw=0.5, alpha=0.75)
        ax.set_title(scan.id)
        k = -1
        for i in range((n / nperrow) + 1):
            for j in range(nperrow):
                k += 1
                try:
                    product_scan = products[k]
                except IndexError:
                    break
                ax = figure.add_subplot(gs[i + 1, j])
                if scan.is_profile:
                    draw_raw(scan.arrays, ax=ax, alpha=0.8)
                pinfo = product_scan.precursor_information
                lower, upper = (product_scan.isolation_window.lower_bound - 2,
                                product_scan.isolation_window.upper_bound + 2)
                peak = max(scan.peak_set.between(lower + 1.2, upper - 1.2), key=lambda x: x.intensity)
                if pinfo.extracted_charge != 0:
                    target_mz = pinfo.extracted_mz
                else:
                    target_mz = pinfo.mz

                draw_peaklist(scan.peak_set, ax=ax, alpha=0.5, lw=0.5)
                if scan.deconvoluted_peak_set:
                    draw_peaklist(
                        scan.deconvoluted_peak_set.between(
                            lower - 1.2, upper + 1.2, use_mz=True),
                        ax=ax, alpha=0.9, color='blue')

                ax.set_ylim(0, peak.intensity * 1.25)
                ax.set_xlim(lower, upper)
                upper_intensity = peak.intensity
                ax.vlines(target_mz, 0, upper_intensity * 1.5, alpha=0.50,
                          color='red', lw=1)
                if product_scan.isolation_window.lower != 0:
                    ax.vlines(product_scan.isolation_window.lower_bound, 0,
                              upper_intensity * 1.5, linestyle='--', alpha=0.5)
                if product_scan.isolation_window.upper != 0:
                    ax.vlines(product_scan.isolation_window.upper_bound, 0,
                              upper_intensity * 1.5, linestyle='--', alpha=0.5)
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                ax.set_ylabel("")
                ax.set_xlabel("")
                if pinfo.charge == "ChargeNotProvided":
                    charge = 0
                else:
                    charge = pinfo.charge
                ax.set_title("%0.3f @ %d" % (target_mz, charge))
        fig = ax.figure
        fig.set_figheight(fig.get_figheight() * 2)
        fig.tight_layout()
        return ax

except (RuntimeError, ImportError):
    has_plot = False

    def annotate_scan(scan, products, nperrow=4, ax=None):
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

    def __contains__(self, i):
        return i in self.name

    def __getitem__(self, i):
        return self.name[i]

    def __iter__(self):
        return iter(self.name)

    def __len__(self):
        return len(self.name)


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


class LRUDict(object):
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
