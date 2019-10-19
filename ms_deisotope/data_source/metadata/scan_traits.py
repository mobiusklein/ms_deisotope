from collections import namedtuple, MutableSequence

from ms_deisotope.utils import _MappingOverAttributeProxy
from .cv import Term, TermSet


def _isclose(x, y, atol=1e-3):
    return abs(x - y) < atol


class IsolationWindow(namedtuple("IsolationWindow", ['lower', 'target', 'upper'])):
    r"""Describes the m/z interval a precursor ion was isolated from in the precursor scan

    An :class:`IsolationWindow` instance is comparable, hashable, and orderable. It is based
    on a :class:`namedtuple`, and supports the same interface and methods as a :class:`tuple`.

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
    width: float
        The sum of :attr:`lower` and :attr:`upper`, the total m/z space spanned by
        the window
    """

    __slots__ = ()

    @classmethod
    def make_empty(cls, point):
        """A helper method for instantiating an empty isolation window centered at
        `point`

        Parameters
        ----------
        point : float
            The point to center the new empty window at.

        Returns
        -------
        :class:`IsolationWindow`
        """
        return cls(0, point, 0)

    @property
    def lower_bound(self):
        """The m/z coordinate of the lower bound :attr:`target` - :attr:`lower`
        """
        return self.target - self.lower

    @property
    def upper_bound(self):
        """The m/z coordinate of the upper bound :attr:`target` + :attr:`upper`
        """
        return self.target + self.upper

    @property
    def width(self):
        '''The sum of :attr:`lower` and :attr:`upper`, the total m/z space spanned by
        the window'''
        return self.lower + self.upper

    def __contains__(self, x):
        return self.spans(x, 0.1)

    def spans(self, x, tolerance=0.1):
        """Test whether `x` is within the interval defined by this window's bounds.

        This method permits a small amount of error (controlled by `tolerance`) when
        computing the bounds.

        Equivalent to: :code:`(self.lower_bound - tolerance) <= x <= (self.upper_bound + tolerance)`

        Parameters
        ----------
        x : float
            The number to query
        tolerance : float, optional
            The amount of error to accept when computing the bounds (the default is 0.1)

        Returns
        -------
        :class:`bool`
        """
        return (self.lower_bound - tolerance) <= x <= (self.upper_bound + tolerance)

    def is_empty(self):
        """Tests if the window is empty (e.g. its upper bound is equal to its lower bound)

        Returns
        -------
        :class:`bool`
        """
        if self.lower is None:
            return self.upper is None
        return self.lower == self.upper == 0.0

    def __nonzero__(self):
        return not self.is_empty()

    def __bool__(self):
        return self.__nonzero__()

    def __eq__(self, other):
        return _isclose(self.lower, other.lower) and _isclose(
            self.upper, other.upper) and _isclose(
                self.target, other.target)

    def __ne__(self, other):
        return not (self == other)

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)


class ScanAcquisitionInformation(MutableSequence):
    """Describes the set distinct scans along the measurable range by
    the instrument that were performed to produce the acquired data

    :class:`ScanAcquisitionInformation` objects implement the :class:`MutableSequence`
    interface, acting as a sequence of :class:`ScanEventInformation` objects.

    Attributes
    ----------
    combination : :class:`str`
        A controlled vocabulary string describing the way in which scans were combined
    scan_list : :class:`list` of :class:`ScanEventInformation`
        The list of scan events performed
    """

    __slots__ = ("combination", "scan_list")

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __init__(self, combination, scan_list):  # pylint: disable=super-init-not-called
        self.combination = combination
        self.scan_list = scan_list

    def __getitem__(self, i):
        return self.scan_list[i]

    def __setitem__(self, i, value):
        self.scan_list[i] = value

    def __delitem__(self, i):
        del self.scan_list[i]

    def insert(self, index, value):
        """Insert an event at a given position

        Parameters
        ----------
        index : int
            The position to add the new event
        value : :class:`ScanEventInformation`
            The event to insert
        """
        if not isinstance(value, ScanEventInformation):
            raise TypeError("Expected ScanEventInformation but got %r" % (type(value), ))
        self.scan_list.insert(index, value)

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

    __slots__ = ("start_time", "window_list", "drift_time", "injection_time", "traits")

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __init__(self, start_time, window_list, drift_time=None, injection_time=None, traits=None, **options):
        self.start_time = start_time
        self.window_list = window_list or []
        self.drift_time = drift_time
        self.injection_time = injection_time
        self.traits = traits or {}
        self.traits.update(options)

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
        if self.injection_time is not None:
            tail += ', injection_time={}'.format(self.injection_time)
        if self.traits:
            tail += ", traits={}".format(self.traits)
        form = template.format(self.start_time, self.window_list, tail)
        return form

    def __eq__(self, other):
        eq = _isclose(self.start_time, other.start_time) and (
            self.window_list == other.window_list)
        if not eq:
            return False
        if self.has_ion_mobility() != other.has_ion_mobility():
            return False
        if self.has_ion_mobility() and not _isclose(self.drift_time, other.drift_time):
            return False
        return True

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
    __slots__ = ()

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

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
        return _isclose(self.lower, other.lower) and _isclose(
            self.upper, other.upper)

    def __ne__(self, other):
        return not (self == other)


class ScanAttribute(Term):
    __slots__ = ()


scan_attributes = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('scan attribute', term_cls_name="ScanAttribute", writer=cog.out)
# ]]]
scan_attributes = TermSet([
    ScanAttribute(u'ion mobility attribute', u'MS:1002892',
                  (u'An attribute describing ion mobility searches.'),
                  'scan attribute',
                  [u'scan attribute', u'PSM-level attribute', u'object attribute', u'single identification result attribute', u'identification attribute', u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'mass resolution', u'MS:1000011',
                  (u'Smallest mass difference between two equal magnitude peaks'
                   u'so that the valley between them is a specified fraction of'
                   u'the peak height.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'scan start time', u'MS:1000016',
                  (u'The time that an analyzer started a scan, relative to the'
                   u'start of the MS run.'),
                  'scan attribute',
                  [u'scan attribute', u'PSM-level attribute', u'object attribute', u'single identification result attribute', u'identification attribute', u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'scan rate', u'MS:1000015',
                  (u'Rate in Th/sec for scanning analyzers.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'elution time', u'MS:1000826',
                  (u'The time of elution from all used chromatographic columns'
                   u'(one or more) in the chromatographic separation step,'
                   u'relative to the start of the chromatography.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'mass resolving power', u'MS:1000800',
                  (u'The observed mass divided by the difference between two'
                   u'masses that can be separated: m/dm. The procedure by which'
                   u'dm was obtained and the mass at which the measurement was'
                   u'made should be reported.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'analyzer scan offset', u'MS:1000803',
                  (u'Offset between two analyzers in a constant neutral loss or'
                   u'neutral gain scan. The value corresponds to the neutral loss'
                   u'or neutral gain value.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'ion injection time', u'MS:1000927',
                  (u'The length of time spent filling an ion trapping device.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'instrument specific scan attribute', u'MS:1002527',
                  (u'Instrument specific scan properties that are associated with'
                   u'a value.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'interchannel delay', u'MS:1000880',
                  (u'The duration of intervals between scanning, during which the'
                   u'instrument configuration is switched.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'dwell time', u'MS:1000502',
                  (u'The time spent gathering data across a peak.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'filter string', u'MS:1000512',
                  (u'A string unique to Thermo instrument describing instrument'
                   u'settings for the scan.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'preset scan configuration', u'MS:1000616',
                  (u'A user-defined scan configuration that specifies the'
                   u'instrumental settings in which a spectrum is acquired. An'
                   u'instrument may cycle through a list of preset scan'
                   u'configurations to acquire data. This is a more generic term'
                   u'for the Thermo \\"scan event\\", which is defined in the'
                   u'Thermo Xcalibur glossary as: a mass spectrometer scan that'
                   u'is defined by choosing the necessary scan parameter'
                   u'settings. Multiple scan events can be defined for each'
                   u'segment of time.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'second column elution time', u'MS:1002083',
                  (u'The time of elution from the second chromatographic column'
                   u'in the chromatographic separation step, relative to the'
                   u'start of the chromatography on the second column.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'first column elution time', u'MS:1002082',
                  (u'The time of elution from the first chromatographic column in'
                   u'the chromatographic separation step, relative to the start'
                   u'of chromatography on the first column.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'inverse reduced ion mobility', u'MS:1002815',
                  (u'Ion mobility measurement for an ion or spectrum of ions as'
                   u'measured in an ion mobility mass spectrometer. This might'
                   u'refer to the central value of a bin into which all ions'
                   u'within a narrow range of mobilities have been aggregated.'),
                  'scan attribute',
                  [u'ion selection attribute', u'ion mobility attribute', u'object attribute', u'scan attribute', u'PSM-level attribute', u'single identification result attribute', u'identification attribute', u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'ion mobility drift time', u'MS:1002476',
                  (u'Drift time of an ion or spectrum of ions as measured in an'
                   u'ion mobility mass spectrometer. This time might refer to the'
                   u'central value of a bin into which all ions within a narrow'
                   u'range of drift time have been aggregated.'),
                  'scan attribute',
                  [u'ion selection attribute', u'ion mobility attribute', u'object attribute', u'scan attribute', u'PSM-level attribute', u'single identification result attribute', u'identification attribute', u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'FAIMS compensation voltage', u'MS:1001581',
                  (u'The DC potential applied to the asymmetric waveform in FAIMS'
                   u'that compensates for the difference between high and low'
                   u'field mobility of an ion.'),
                  'scan attribute',
                  [u'ion mobility attribute', u'scan attribute', u'PSM-level attribute', u'object attribute', u'single identification result attribute', u'identification attribute', u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'synchronous prefilter selection', u'MS:1002528',
                  (u'Synchronous prefilter selection.'),
                  'scan attribute',
                  [u'instrument specific scan attribute', u'scan attribute', u'object attribute']),
])
# [[[end]]]

__all__ = [
    "IsolationWindow", "ScanAcquisitionInformation", "ScanEventInformation",
    "ScanWindow", "ScanAttribute", "scan_attributes"
]
