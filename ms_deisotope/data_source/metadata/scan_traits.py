"""A module defining types for describing instrument scans and windows,
as well as a :class:`~.Term` subclass for the `scan attribute` family.
"""

from collections import namedtuple

try:
    from collections.abc import MutableSequence
except ImportError:
    from collections import MutableSequence


from pyteomics.auxiliary import unitfloat  # pylint: disable=unused-import

from ms_deisotope.utils import _MappingOverAttributeProxy
from .cv import Term, TermSet


def _isclose(x, y, atol=1e-3):
    return abs(x - y) < atol


_IsolationWindowBase = namedtuple(
    "IsolationWindow", ['lower', 'target', 'upper'])


class IsolationWindow(_IsolationWindowBase):
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
    __hash__ = _IsolationWindowBase.__hash__

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

    def __reduce__(self):
        return self.__class__, (self.lower, self.target, self.upper)


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
            raise TypeError(
                "Expected ScanEventInformation but got %r" % (type(value), ))
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

    def __reduce__(self):
        return self.__class__, (self.combination, self.scan_list)

    def copy(self):
        return self.__class__(self.combination, [se.copy() for se in self.scan_list])


class _IonMobilityMixin(object):
    __slots__ = ()

    @property
    def _ion_mobility(self):
        return IonMobilityMethods(self)

    def has_ion_mobility(self):
        """Check whether the scan has ion mobility measure.

        Returns
        -------
        bool
        """
        return self._ion_mobility.has_ion_mobility()

    @property
    def drift_time(self):
        """Fetch the drift time quantity of the scan

        Returns
        -------
        value: float or :const:`None`
            The measured drift time, or :const:`None`
        """
        return self._ion_mobility.drift_time()

    @drift_time.setter
    def drift_time(self, value):
        if value is None or value == 0:
            return
        else:
            raise ValueError("Cannot set drift time directly yet.")

    @property
    def ion_mobility_type(self):
        """Fetch the ion mobility type of the scan.

        Returns
        -------
        ims_type: :class:`ScanAttribute` or :const:`None`
            The ion mobility type, or :const:`None`
        """
        return self._ion_mobility.ion_mobility_type()


class ScanEventInformation(_IonMobilityMixin):
    """Describe a single instrument scan, one or more of which compose
    an actual :class:`Scan` object which describes a spectrum.

    Attributes
    ----------
    start_time: :class:`unitfloat`
        The time since the experiment began when the scan started being acquired.
    window_list: :class:`list` of :class:`ScanWindow`
        The intervals over the m/z axis that were sampled from during this scan.
    injection_time: :class:`unitfloat`
        The time spent injecting ions into the analyzer.
    traits: :class:`dict`
        A mapping from :class:`ScanAttribute` to optional values describing this scan event.
    drift_time: :class:`unitfloat`
        The IMS drift time, if present, :const:`None` otherwise. Determined from :attr:`traits`.
    ion_mobility_type: :class:`ScanAttribute`
        The IMS mechanism type, if present, :const:`None` otherwise. Determined from :attr:`traits`.
    """

    __slots__ = ("start_time", "window_list", "injection_time", "traits")

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __init__(self, start_time, window_list, drift_time=None, injection_time=None, traits=None, **options):
        self.traits = traits or {}
        self.traits.update(options)
        self.start_time = start_time
        self.window_list = window_list or []
        self.drift_time = drift_time
        self.injection_time = injection_time

    @property
    def scan_configuration(self):
        return self.traits.get('preset scan configuration')

    def copy(self):
        return self.__class__(self.start_time, self.window_list[:], None, self.injection_time, self.traits.copy())

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
        """Create a :class:`ScanWindow` spanning from the lowest m/z of the lowest window
        to the highest m/z of the highest window.

        Returns
        -------
        :class:`ScanWindow`
        """
        low = float('inf')
        high = 0
        for window in self:
            low = min(low, window.lower)
            high = max(high, window.upper)
        return ScanWindow(low, high)

    def __reduce__(self):
        return self.__class__, (self.start_time, self.window_list, None, self.injection_time, self.traits)


class ScanWindow(namedtuple("ScanWindow", ['lower', 'upper'])):
    """Represent a single contiguous m/z interval that was sampled from during a scan.

    This type supports spanning testing using the :meth:`__contains__` method, and is
    equality comparable and hashable.

    Attributes
    ----------
    lower: float
        The lower bound of the window
    upper: float
        The upper bound of the window
    """
    __slots__ = ()

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __reduce__(self):
        return self.__class__, (self.lower, self.upper)

    def __contains__(self, i):
        return self.lower <= i <= self.upper

    def is_empty(self):
        """Test whether the window's bounds are both equal to :const:`None`
        or `0`.

        Returns
        -------
        bool
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
            self.upper, other.upper)

    def __hash__(self):
        return super(ScanWindow, self).__hash__()

    def __ne__(self, other):
        return not (self == other)


class ScanAttribute(Term):
    """Describes a single trait or attribute belonging to a scan,
    such as injection time, filter string, or instrument configuration.
    """
    __slots__ = ()


scan_attributes = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('scan attribute', term_cls_name="ScanAttribute", writer=cog.out)
# ]]]
scan_attributes = TermSet([
    ScanAttribute(u'mass resolution', u'MS:1000011',
                  (u'Smallest mass difference between two equal magnitude peaks'
                   u'so that the valley between them is a specified fraction of'
                   u'the peak height.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'scan rate', u'MS:1000015',
                  (u'Rate in Th/sec for scanning analyzers.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'scan start time', u'MS:1000016',
                  (u'The time that an analyzer started a scan, relative to the'
                   u'start of the MS run.'),
                  'scan attribute',
                  [u'scan attribute', u'PSM-level attribute', u'object attribute',
                   u'single identification result attribute', u'identification attribute',
                   u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'zoom scan', u'MS:1000497',
                  (u'Special scan mode where data with improved resolution is'
                   u'acquired. This is typically achieved by scanning a more'
                   u'narrow m/z window or scanning with a lower scan rate.'),
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
    ScanAttribute(u'elution time', u'MS:1000826',
                  (u'The time of elution from all used chromatographic columns'
                   u'(one or more) in the chromatographic separation step,'
                   u'relative to the start of the chromatography.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'interchannel delay', u'MS:1000880',
                  (u'The duration of intervals between scanning, during which the'
                   u'instrument configuration is switched.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'ion injection time', u'MS:1000927',
                  (u'The length of time spent filling an ion trapping device.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'first column elution time', u'MS:1002082',
                  (u'The time of elution from the first chromatographic column in'
                   u'the chromatographic separation step, relative to the start'
                   u'of chromatography on the first column.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'second column elution time', u'MS:1002083',
                  (u'The time of elution from the second chromatographic column'
                   u'in the chromatographic separation step, relative to the'
                   u'start of the chromatography on the second column.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'instrument specific scan attribute', u'MS:1002527',
                  (u'Instrument specific scan properties that are associated with'
                   u'a value.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'ion mobility attribute', u'MS:1002892',
                  (u'An attribute describing ion mobility searches.'),
                  'scan attribute',
                  [u'scan attribute', u'PSM-level attribute', u'object attribute',
                   u'single identification result attribute', u'identification attribute',
                   u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'scan number', u'MS:1003057',
                  (u'Ordinal number of the scan indicating its order of'
                   u'acquisition within a mass spectrometry acquisition run.'),
                  'scan attribute',
                  [u'scan attribute', u'object attribute']),
    ScanAttribute(u'synchronous prefilter selection', u'MS:1002528',
                  (u'Synchronous prefilter selection.'),
                  'scan attribute',
                  [u'instrument specific scan attribute', u'scan attribute', u'object attribute']),
    ScanAttribute(u'FAIMS compensation voltage', u'MS:1001581',
                  (u'The DC potential applied to the asymmetric waveform in FAIMS'
                   u'that compensates for the difference between high and low'
                   u'field mobility of an ion.'),
                  'scan attribute',
                  [u'ion mobility attribute', u'scan attribute', u'PSM-level attribute', u'object attribute',
                   u'single identification result attribute', u'identification attribute', u'analysis attribute',
                   u'spectrum identification result details']),
    ScanAttribute(u'ion mobility drift time', u'MS:1002476',
                  (u'Drift time of an ion or spectrum of ions as measured in an'
                   u'ion mobility mass spectrometer. This time might refer to the'
                   u'central value of a bin into which all ions within a narrow'
                   u'range of drift time have been aggregated.'),
                  'scan attribute',
                  [u'ion selection attribute', u'ion mobility attribute', u'object attribute', u'scan attribute',
                   u'PSM-level attribute', u'single identification result attribute', u'identification attribute',
                   u'analysis attribute', u'spectrum identification result details']),
    ScanAttribute(u'inverse reduced ion mobility', u'MS:1002815',
                  (u'Ion mobility measurement for an ion or spectrum of ions as'
                   u'measured in an ion mobility mass spectrometer. This might'
                   u'refer to the central value of a bin into which all ions'
                   u'within a narrow range of mobilities have been aggregated.'),
                  'scan attribute',
                  [u'ion selection attribute', u'ion mobility attribute', u'object attribute', u'scan attribute',
                   u'PSM-level attribute', u'single identification result attribute', u'identification attribute',
                   u'analysis attribute', u'spectrum identification result details']),
])
# [[[end]]]


FAIMS_compensation_voltage = scan_attributes['MS:1001581']
ion_mobility_drift_time = scan_attributes['MS:1002476']
inverse_reduced_ion_mobility = scan_attributes['MS:1002815']
ion_mobility_attribute = scan_attributes['MS:1002892']

ION_MOBILITY_TYPES = {
    FAIMS_compensation_voltage,
    ion_mobility_drift_time,
    inverse_reduced_ion_mobility,
    ion_mobility_attribute,
}


class IonMobilityMethods(object):
    """Determine the ion mobility measure of a scan event.

    This interface attempts to find FAIMS compensation voltage,
    ion mobility drift time, and inverse reduced ion mobility.

    Attributes
    ----------
    scan_event: :class:`ScanEventInformation`
        The scan to interpret.
    """

    def __init__(self, scan_event):
        self.scan_event = scan_event

    def _get_traits(self):
        return self.scan_event.traits

    def get(self):
        """Find the first ion mobility trait of the scan.

        Returns
        -------
        ims_type: :class:`ScanAttribute`
            The ion mobility type.
        value: float
            The drift time value.
        """
        for ims_type in ION_MOBILITY_TYPES:
            result = self._get_traits().get(ims_type.name)
            if result is not None:
                return (ims_type, result)
        return None

    def has_ion_mobility(self):
        """Check whether the scan has ion mobility measure.

        Returns
        -------
        bool
        """
        return self.get() is not None

    def drift_time(self):
        """Fetch the drift time quantity of the scan

        Returns
        -------
        value: float or :const:`None`
            The measured drift time, or :const:`None`
        """
        _ims_type_value = self.get()
        if _ims_type_value is None:
            return None
        _ims_type, value = _ims_type_value
        return value

    def add_ion_mobility(self, ion_mobility_type, drift_time):
        """Set the drift time value for a specific type of ion mobility on
        the scan event.

        Parameters
        ----------
        ion_mobility_type : :class:`unitfloat`
            The type of ion mobility to set the drift time for
        drift_time : :class:`unitfloat`
            The drift time, with appropriate units.
        """
        self._get_traits()[ion_mobility_type] = drift_time

    def remove_ion_mobility_type(self, ion_mobility_type):
        """Remove a specific type of ion mobility from the scan.

        Parameters
        ----------
        ion_mobility_type : :class:`ScanAttribute`
            The type of ion mobility to remove.
        """
        self._get_traits().pop(ion_mobility_type)

    def ion_mobility_type(self):
        """Fetch the ion mobility type of the scan.

        Returns
        -------
        ims_type: :class:`ScanAttribute` or :const:`None`
            The ion mobility type, or :const:`None`
        """
        _ims_type_value = self.get()
        if _ims_type_value is None:
            return None
        ims_type, _value = _ims_type_value
        return ims_type


binary_data_arrays = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('binary data array', writer=cog.out)
# ]]]
binary_data_arrays = TermSet([
    Term(u'm/z array', u'MS:1000514',
         (u'A data array of m/z values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'intensity array', u'MS:1000515',
         (u'A data array of intensity values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'charge array', u'MS:1000516',
         (u'A data array of charge values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'signal to noise array', u'MS:1000517',
         (u'A data array of signal-to-noise values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'time array', u'MS:1000595',
         (u'A data array of relative time offset values from a reference'
          u'time.'),
         'binary data array',
         [u'binary data array']),
    Term(u'wavelength array', u'MS:1000617',
         (u'A data array of electromagnetic radiation wavelength values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'non-standard data array', u'MS:1000786',
         (u'A data array that contains data not covered by any other'
          u'term in this group. Please do not use this term, if the'
          u'binary data array type might be commonly used - contact the'
          u'PSI-MS working group in order to have another CV term added.'),
         'binary data array',
         [u'binary data array']),
    Term(u'flow rate array', u'MS:1000820',
         (u'A data array of flow rate measurements.'),
         'binary data array',
         [u'binary data array']),
    Term(u'pressure array', u'MS:1000821',
         (u'A data array of pressure measurements.'),
         'binary data array',
         [u'binary data array']),
    Term(u'temperature array', u'MS:1000822',
         (u'A data array of temperature measurements.'),
         'binary data array',
         [u'binary data array']),
    Term(u'mean charge array', u'MS:1002478',
         (u'Array of mean charge values where the mean charge is'
          u'calculated as a weighted mean of the charges of individual'
          u'peaks that are aggregated into a processed spectrum.'),
         'binary data array',
         [u'binary data array']),
    Term(u'resolution array', u'MS:1002529',
         (u'A data array of resolution values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'baseline array', u'MS:1002530',
         (u'A data array of signal baseline values (the signal in the'
          u'absence of analytes).'),
         'binary data array',
         [u'binary data array']),
    Term(u'noise array', u'MS:1002742',
         (u'A data array of noise values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'sampled noise m/z array', u'MS:1002743',
         (u'A data array of parallel, independent m/z values for a'
          u'sampling of noise across a spectrum (typically much smaller'
          u'than MS:1000514, the m/z array).'),
         'binary data array',
         [u'binary data array']),
    Term(u'sampled noise intensity array', u'MS:1002744',
         (u'A data array of intensity values for the amplitude of noise'
          u'variation superposed on the baseline (MS:1002745) across a'
          u'spectrum (for use with MS:1002743, sampled noise m/z array).'),
         'binary data array',
         [u'binary data array']),
    Term(u'sampled noise baseline array', u'MS:1002745',
         (u'A data array of baseline intensity values (the intensity in'
          u'the absence of analytes) for a sampling of noise across a'
          u'spectrum (for use with MS:1002743, sampled noise m/z array).'),
         'binary data array',
         [u'binary data array']),
    Term(u'ion mobility array', u'MS:1002893',
         (u'An array of ion mobility data.'),
         'binary data array',
         [u'binary data array']),
    Term(u'mass array', u'MS:1003143',
         (u'A data array of mass values.'),
         'binary data array',
         [u'binary data array']),
    Term(u'mean drift time array', u'MS:1002477',
         (u'Array of drift times, averaged from a matrix of binned m/z'
          u'and drift time values, corresponding to spectrum of'
          u'individual peaks encoded with an m/z array.'),
         'binary data array',
         [u'ion mobility array', u'binary data array']),
    Term(u'mean ion mobility array', u'MS:1002816',
         (u'Array of drift times, averaged from a matrix of binned m/z'
          u'and ion mobility values, corresponding to a spectrum of'
          u'individual peaks encoded with an m/z array.'),
         'binary data array',
         [u'ion mobility array', u'binary data array']),
    Term(u'mean inverse reduced ion mobility array', u'MS:1003006',
         (u'Array of inverse reduced ion mobilities, averaged from a'
          u'matrix of binned m/z and ion mobility values, corresponding'
          u'to a spectrum of individual peaks encoded with an m/z array.'),
         'binary data array',
         [u'ion mobility array', u'binary data array']),
    Term(u'raw ion mobility array', u'MS:1003007',
         (u'Array of raw drift times.'),
         'binary data array',
         [u'ion mobility array', u'binary data array']),
    Term(u'raw inverse reduced ion mobility array', u'MS:1003008',
         (u'Array of raw inverse reduced ion mobilities.'),
         'binary data array',
         [u'ion mobility array', u'binary data array']),
])
# [[[end]]]


__all__ = [
    "IsolationWindow", "ScanAcquisitionInformation", "ScanEventInformation",
    "ScanWindow", "ScanAttribute", "scan_attributes", "array_types"
]
