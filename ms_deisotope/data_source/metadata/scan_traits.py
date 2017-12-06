from collections import namedtuple


class IsolationWindow(namedtuple("IsolationWindow", ['lower', 'target', 'upper'])):

    @property
    def lower_bound(self):
        return self.target - self.lower

    @property
    def upper_bound(self):
        return self.target + self.upper

    def __contains__(self, x):
        return self.lower_bound <= x <= self.upper_bound

    def is_empty(self):
        if self.lower is None:
            return self.upper is None
        return self.lower == self.upper == 0.0

    def __nonzero__(self):
        return not self.is_empty()


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
