try:
    range = xrange
except:
    range = range


def simple_repr(self):  # pragma: no cover
    template = "{self.__class__.__name__}({d})"
    d = [
        "%s=%r" % (k, v) if v is not self else "(...)" for k, v in sorted(
            self.__dict__.items(), key=lambda x: x[0])
        if (not k.startswith("_") and not callable(v)) and not (k == "signal")]
    return template.format(self=self, d=', '.join(d))


class Base(object):
    __repr__ = simple_repr


def ppm_error(x, y):
    return (x - y) / y

try:
    has_plot = True
    from matplotlib import pyplot as plt

    def draw_peaklist(peaklist, ax=None, **kwargs):
        kwargs.setdefault("width", 0.1)
        if ax is None:
            fig, ax = plt.subplots(1)
        ax.bar([p.mz for p in peaklist], [p.intensity for p in peaklist], **kwargs)
        ax.xaxis.set_ticks_position('none')
        return ax

except ImportError:
    has_plot = False
