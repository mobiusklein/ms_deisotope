import operator

try:
    has_plot = True
    from matplotlib import pyplot as plt
    from ms_peak_picker.utils import draw_peaklist

except ImportError:
    has_plot = False

try:
    range = xrange
except:
    range = range


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


def ppm_error(x, y):
    return (x - y) / y


def dict_proxy(attribute):
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
