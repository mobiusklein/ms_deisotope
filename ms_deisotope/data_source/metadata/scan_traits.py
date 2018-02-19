from collections import namedtuple


def isclose(x, y, atol=1e-3):
    return abs(x - y) < atol


class IsolationWindow(namedtuple("IsolationWindow", ['lower', 'target', 'upper'])):
    """Describes the m/z interval a precursor ion was isolated from in the precursor scan

    Attributes
    ----------
    lower: float
        The distance from the center of the isolation window towards 0
    target: float
        The center of the isolation window in m/z
    upper: float
        The distance from the center of the isolation window towards :math:`\infty`
    lower_bound: float
        The m/z coordinate of the lower bound :attr:`target` - :attr:`lower`
    upper_bound: float
        The m/z coordinate of the upper bound :attr:`target` + :attr:`upper`
    """

    @property
    def lower_bound(self):
        return self.target - self.lower

    @property
    def upper_bound(self):
        return self.target + self.upper

    def __contains__(self, x):
        return self.spans(x, 0.1)

    def spans(self, x, tolerance=0.1):
        return (self.lower_bound - tolerance) <= x <= (self.upper_bound + tolerance)

    def is_empty(self):
        if self.lower is None:
            return self.upper is None
        return self.lower == self.upper == 0.0

    def __nonzero__(self):
        return not self.is_empty()

    def __bool__(self):
        return self.__nonzero__()

    def __eq__(self, other):
        return isclose(self.lower, other.lower) and isclose(
            self.upper, other.upper) and isclose(
            self.target, other.target)

    def __ne__(self, other):
        return not (self == other)


class ScanAcquisitionInformation(object):
    def __init__(self, combination, scan_list):
        self.combination = combination
        self.scan_list = scan_list

    def __getitem__(self, i):
        return self.scan_list[i]

    def __iter__(self):
        return iter(self.scan_list)

    def __len__(self):
        return len(self.scan_list)

    def __repr__(self):
        return "ScanAcquisitionInformation(combination=%r, scan_list=%r)" % (
            self.combination, self.scan_list)

    def __eq__(self, other):
        return self.combination == other.combination and self.scan_list == other.scan_list

    def __ne__(self, other):
        return not (self == other)


class ScanEventInformation(object):
    def __init__(self, start_time, window_list, drift_time=None):
        self.start_time = start_time
        self.window_list = window_list or []
        self.drift_time = drift_time

    def has_ion_mobility(self):
        return self.drift_time is not None and self.drift_time > 0

    def __getitem__(self, i):
        return self.window_list[i]

    def __iter__(self):
        return iter(self.window_list)

    def __len__(self):
        return len(self.window_list)

    def __repr__(self):
        template = "ScanEventInformation(start_time={}, window_list={}{})"
        if self.has_ion_mobility():
            tail = ", drift_time={}".format(self.drift_time)
        else:
            tail = ''
        form = template.format(self.start_time, self.window_list, tail)
        return form

    def __eq__(self, other):
        eq = isclose(self.start_time, other.start_time) and (
            self.window_list == other.window_list) and isclose(
            self.drift_time, other.drift_time)
        return eq

    def __ne__(self, other):
        return not (self == other)

    def total_scan_window(self):
        low = float('inf')
        high = 0
        for window in self:
            low = min(low, window.lower)
            high = max(high, window.upper)
        return ScanWindow(low, high)


class ScanWindow(namedtuple("ScanWindow", ['lower', 'upper'])):
    def __contains__(self, i):
        return self.lower <= i <= self.upper

    def is_empty(self):
        if self.lower is None:
            return self.upper is None
        return self.lower == self.upper == 0.0

    def __nonzero__(self):
        return not self.is_empty()

    def __bool__(self):
        return self.__nonzero__()

    def __eq__(self, other):
        return isclose(self.lower, other.lower) and isclose(
            self.upper, other.upper)

    def __ne__(self, other):
        return not (self == other)
