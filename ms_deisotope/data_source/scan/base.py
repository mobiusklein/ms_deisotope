'''Represent the basic structures of a mass spectrum and its processed contents,
and provide an interface for manipulating that data.
'''
import warnings

from collections import namedtuple

import numpy as np

from ms_peak_picker import average_signal
from ms_deisotope.averagine import neutral_mass, mass_charge_ratio


try:
    from ms_deisotope.plot import annotate_scan as _annotate_precursors, draw_raw
except ImportError:
    def _missing_matplotlib(*args, **kwargs):
        raise ImportError(
            "This method requires matplotlib. Please install it.")
    _annotate_precursors = _missing_matplotlib
    draw_raw = _missing_matplotlib

from ms_deisotope.utils import Constant


DEFAULT_CHARGE_WHEN_NOT_RESOLVED = 1
ChargeNotProvided = Constant("ChargeNotProvided")


class ScanBunch(namedtuple("ScanBunch", ["precursor", "products"])):
    """Represents a single MS1 scan and all MSn scans derived from it,
    or a collection of related MSn scans.

    Attributes
    ----------
    precursor: :class:`~.ScanBase`
        A single MS1 scan which may have undergone MSn
    products: list
        A list of 0 or more :class:`~.ScanBase` objects which were derived
        from :attr:`precursor` or another element of this list derived
        from it.
    """

    def __new__(cls, *args, **kwargs):  # pylint: disable=super-on-old-class
        inst = super(ScanBunch, cls).__new__(cls, *args, **kwargs)
        inst._id_map = {}
        if inst.precursor is not None:
            inst._id_map[inst.precursor.id] = inst.precursor
        for scan in inst.products:
            inst._id_map[scan.id] = scan
        return inst

    def precursor_for(self, scan):
        """Find the precursor :class:`~.ScanBase` instance
        for the given scan object

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The MSn scan to look for the MSn-1 scan for

        Returns
        -------
        :class:`~.ScanBase`
        """
        if scan.precursor_information is not None:
            scan_id = scan.precursor_information.precursor_scan_id
            return self.get_scan_by_id(scan_id)
        return None

    def get_scan_by_id(self, scan_id):
        """Retrieve the scan object for the specified scan id from this
        group in memory.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        :class:`~.ScanBase`
        """
        return self._id_map[scan_id]

    def annotate_precursors(self, nperrow=4, ax=None):
        '''Plot the spectra in this group as a grid, with the full
        MS1 spectrum in profile in the top row, and each MSn spectrum's
        precursor ion revealed in a grid panel below, with isolation
        window and selected ion/monoisotopic peak annotated.

        Parameters
        ----------
        nperrow: :class:`int`
            The number of precursors to annotate per row
            in the grid.
        ax: :class:`matplotlib._axes.Axes`, optional
            The axis to draw on. If not provided, a new figure
            will be created, along with a new axis.

        Returns
        -------
        :class:`matplotlib._axes.Axes`
        '''
        return _annotate_precursors(
            self.precursor, self.products, nperrow=nperrow, ax=ax)

    def _repr_pretty_(self, p, cycle):  # pragma: no cover
        if cycle:
            p.text("ScanBunch(...)")
            return
        p.text("ScanBunch(\n")
        with p.group(2):
            with p.group(4, "precursor=\n"):
                p.pretty(self.precursor)
            with p.group(4, ",\nproducts=\n"):
                p.pretty(self.products)
        p.text(")")

    def pack(self):
        '''Build a new :class:`ScanBunch` where each scan in it is returned by calling
        :meth:`~.Scan.pack`

        Returns
        -------
        :class:`ScanBunch`
        '''
        return self.__class__(self.precursor.pack(), [
            p.pack() for p in self.products
        ])


class RawDataArrays(namedtuple("RawDataArrays", ['mz', 'intensity'])):
    """Represent the m/z and intensity arrays associated with a raw
    mass spectrum.

    Supports scaling and summing, as well as low level m/z search.

    Thin wrapper around a ``namedtuple``, so this object supports
    the same interfaces as a tuple.

    Attributes
    ----------
    mz: np.ndarray
        The m/z axis of a mass spectrum
    intensity: np.ndarray
        The intensity measured at the corresponding m/z of a mass spectrum
    """

    def plot(self, *args, **kwargs):
        """Draw the profile spectrum described by the
        contained arrays.

        Parameters
        ----------
        ax: :class:`matplotlib._axes.Axes`
            The figure axes onto which to draw the plot. If not provided,
            this will default to the current figure interactively.
        **kwargs
            All keywords are forwarded to :meth:`plot` on ``ax``.

        Returns
        -------
        :class:`matplotlib._axes.Axes`
            The axes drawn on
        """
        ax = draw_raw(self, *args, **kwargs)
        return ax

    def __eq__(self, other):
        try:
            return np.allclose(
                self[0], other[0]) and np.allclose(
                    self[1], other[1])
        except ValueError:
            return False

    def __ne__(self, other):
        return not (self == other)

    def __mul__(self, i):
        return self.__class__(self.mz, self.intensity * i)

    def __div__(self, d):
        return self.__class__(self.mz, self.intensity / d)

    def __add__(self, other):
        if len(self.mz) == len(other.mz) and np.allclose(self.mz, other.mz):
            return self.__class__(self.mz, self.intensity + other.intensity)
        else:
            return self.__class__(*average_signal([self, other])) * 2

    def find_mz(self, mz):
        """Find the nearest index to the query ``mz``

        Parameters
        ----------
        mz : float
            The m/z value to search for

        Returns
        -------
        int
            The index nearest to the query m/z
        """
        n = len(self.mz)
        lo = 0
        hi = n

        while hi != lo:
            mid = int((hi + lo) // 2)
            y = self.mz[mid]
            err = y - mz
            if abs(err) < 0.1:
                best_index = mid
                best_err = abs(err)
                i = mid
                while i >= 0:
                    y = self.mz[i]
                    err = y - mz
                    if err <= -0.1:
                        break
                    abs_err = abs(err)
                    if abs_err < best_err:
                        best_err = abs_err
                        best_index = i
                    i -= 1
                i = mid
                while i < n:
                    y = self.mz[i]
                    err = y - mz
                    if err >= 0.1:
                        break
                    abs_err = abs(err)
                    if abs_err < best_err:
                        best_err = abs_err
                        best_index = i
                    i += 1
                return best_index
            elif hi - lo == 1:
                return mid
            elif err > 0:
                hi = mid
            else:
                lo = mid
        return 0

    def between_mz(self, low, high):
        """Returns a slice of the arrays between ``low`` and ``high``
        m/z

        Parameters
        ----------
        low : float
            The lower bound m/z
        high : float
            The upper bound m/z

        Returns
        -------
        :class:`.RawDataArrays`
        """
        i = self.find_mz(low)
        j = self.find_mz(high) + 1
        if not (low <= self.mz[i] <= high):
            i += 1
        return self.__class__(self.mz[i:j], self.intensity[i:j])


class ScanBase(object):
    '''Abstract base class for Scan-like objects
    '''

    def has_ion_mobility(self):
        '''Check whether this scan has drift time information associated with
        it.

        If this scan has been aggregated, it will only check the first scan in
        the aggregate.
        '''
        acq = self.acquisition_information
        if acq is None:
            return False
        scan_event = acq[0]
        return scan_event.has_ion_mobility()

    @property
    def drift_time(self):
        '''A convenience method to access the first
        scan event to retrieve its drift time.

        Returns
        -------
        float or None
        '''
        acq = self.acquisition_information
        if acq is None:
            return None
        scan_event = acq[0]
        return scan_event.drift_time

    @property
    def scan_id(self):
        '''An alias for :attr:`id`
        '''
        return self.id

    def copy(self, deep=True):
        """Return a deep copy of the :class:`Scan` object
        wrapping the same reference data.

        Returns
        -------
        :class:`Scan`
        """
        return self.clone(deep)

    def __copy__(self):
        return self.clone()

    def __eq__(self, other):
        if other is None:
            return False
        if not isinstance(other, ScanBase):
            return False
        try:
            eq = (self.scan_id == other.scan_id) and (
                abs(self.scan_time - other.scan_time) < 1e-3) and (
                    self.index == other.index) and (
                        self.ms_level == other.ms_level)
            if not eq:
                return False
        except AttributeError:
            return False
        try:
            eq = self.arrays == other.arrays
            if not eq:
                return False
        except AttributeError:
            # ProcessedScan doesn't have an arrays attribute
            pass
        try:
            eq = self.peak_set == other.peak_set
            if not eq:
                return False
        except AttributeError:
            if ((self.peak_set is None and other.peak_set is not None) or (
                    self.peak_set is not None and other.peak_set is None)):
                pass
            else:
                return False

        try:
            eq = self.deconvoluted_peak_set == other.deconvoluted_peak_set
            if not eq:
                return False
        except AttributeError:
            if ((self.deconvoluted_peak_set is None and other.deconvoluted_peak_set is not None) or (
                    self.deconvoluted_peak_set is not None and other.deconvoluted_peak_set is None)):
                pass
            else:
                return False

        eq = self.precursor_information == other.precursor_information
        if not eq:
            return False
        eq = self.isolation_window == other.isolation_window
        if not eq:
            return False
        try:
            a = self.acquisition_information
            b = other.acquisition_information
            if a is not None and b is not None:
                eq = a == b
            else:
                eq = True
            if not eq:
                return False
        except AttributeError:
            pass
        try:
            a = self.activation
            b = other.activation
            if a is not None and b is not None:
                eq = a == b
            else:
                eq = True
            if not eq:
                return False
        except AttributeError:
            pass

        return True

    def __ne__(self, other):
        return not (self == other)

    def bind(self, source):
        '''Attach this object and its other referent members
        to ``source``, letting them load information.
        '''
        if self.precursor_information is not None:
            self.precursor_information.bind(source)
        return self

    def unbind(self):
        '''Detattch this object and its other referent members
        from their currently bound :attr:`source`.

        This may cause errors if more information is requested but is not
        cached, or if requesting another :class:`ScanBase` be loaded.
        '''
        if self.precursor_information is not None:
            self.precursor_information.unbind()
        return self


class PrecursorInformation(object):
    """Store information relating a tandem MS scan to its precursor MS scan.

    .. note:
        The attributes prefixed with `extracted_` refer to the quantities estimated
        from the data, while those unprefixed are the values read directly from the
        data source. These values regularly do not agree. When available, the extracted
        values should be more accurate.

    Attributes
    ----------
    charge : int
        The charge reported in the source metadata
    defaulted : bool
        Whether the information in the extracted fields reflects empirical
        information or fell back on the vendor-reported values.
    extracted_charge : int
        The charge estimated from the source data
    extracted_intensity : float
        The sum of the peak heights of the extracted isotopic pattern
    extracted_neutral_mass : float
        The monoisotopic neutral mass estimated from the source data
    extracted_peak : :class:`.DeconvolutedPeak`
        The deconvoluted peak summarizing the precursor ion
    intensity : float
        The abundance reported in the source metadata
    mz : float
        The m/z reported in the source metadata
    orphan : bool
        Whether there was an isotopic pattern to extract in the precursor scan. Usually
        paired with `defaulted`
    peak : :class:`.FittedPeak`
        The peak nearest :attr:`mz`, and the starting point for estimating information
        about the precursor ion
    precursor_scan_id : str
        The id string for the precursor scan
    source : :class:`ScanIterator`
        Any object implementing the :class:`ScanIterator` interface to be used to look up
        the precursor scan with :attr:`precursor_scan_id`
    """

    def __init__(self, mz, intensity, charge, precursor_scan_id=None, source=None,
                 extracted_neutral_mass=0, extracted_charge=0, extracted_intensity=0,
                 peak=None, extracted_peak=None, defaulted=False, orphan=False,
                 product_scan_id=None, annotations=None, coisolation=None):
        try:
            charge = int(charge)
        except Exception:
            pass
        try:
            extracted_charge = int(extracted_charge)
        except Exception:
            pass
        if not annotations:
            annotations = {}
        if not coisolation:
            coisolation = []

        self.mz = mz
        self.intensity = intensity
        self.charge = charge

        self.precursor_scan_id = precursor_scan_id
        self.source = source

        self.extracted_neutral_mass = extracted_neutral_mass
        self.extracted_charge = extracted_charge
        self.extracted_intensity = extracted_intensity

        self.peak = peak
        self.extracted_peak = extracted_peak
        self.defaulted = defaulted
        self.orphan = orphan
        self.product_scan_id = product_scan_id

        self.annotations = annotations
        self.coisolation = coisolation

    def __repr__(self):
        return "PrecursorInformation(mz=%0.4f/%0.4f, intensity=%0.4f/%0.4f, charge=%r/%r, scan_id=%r)" % (
            self.mz,
            self.extracted_mz if self.extracted_neutral_mass != 0. else 0.,
            self.intensity or 0., self.extracted_intensity or 0., self.charge,
            self.extracted_charge or 0., self.precursor_scan_id)

    def __reduce__(self):
        return self.__class__, (0, 0, 0), self.__getstate__()

    def __getstate__(self):
        # explicitly do not propagate :attr:`source` when serializing.
        return (self.mz, self.intensity, self.charge, self.precursor_scan_id, None, self.extracted_neutral_mass,
                self.extracted_charge, self.extracted_intensity, self.peak, self.extracted_peak,
                self.defaulted, self.orphan, self.product_scan_id, self.annotations, self.coisolation)

    def __setstate__(self, state):
        (self.mz, self.intensity, self.charge, self.precursor_scan_id, self.source, self.extracted_neutral_mass,
         self.extracted_charge, self.extracted_intensity, self.peak, self.extracted_peak,
         self.defaulted, self.orphan, self.product_scan_id) = state[:13]
        if len(state) > 13:
            self.annotations = state[13]
        if len(state) > 14:
            self.coisolation = list(state[14])

    def __eq__(self, other):
        if other is None:
            return False
        eq = self.precursor_scan_id == other.precursor_scan_id
        if not eq:
            return False
        eq = self.product_scan_id == other.product_scan_id
        if not eq:
            return False
        self_fit = self.extracted_neutral_mass != 0
        other_fit = other.extracted_neutral_mass != 0
        self_mass = self.extracted_mz if self_fit else self.mz
        other_mass = other.extracted_mz if other_fit else other.mz
        eq = np.isclose(self_mass, other_mass)
        if not eq:
            return False
        self_charge = self.extracted_charge if self_fit else self.charge
        other_charge = other.extracted_charge if other_fit else other.charge
        eq = self_charge == other_charge
        if not eq:
            return False
        return True

    def __ne__(self, other):
        return not (self == other)

    def bind(self, source):
        '''Attach this object and its other referent members
        to ``source``, letting them load information.
        '''
        self.source = source
        return self

    def unbind(self):
        self.source = None

    def extract(self, peak, override_charge=None):
        '''Populate the extracted attributes of this object from the attributes
        of a :class:`~.DeconvolutedPeak` instance.

        Parameters
        ----------
        peak: :class:`~.DeconvolutedPeak`
            The peak to copy attributes from
        override_charge: :class:`int`, optional
            If provided, this charge will be used instead of the charge of ``peak``
        '''
        self.extracted_neutral_mass = peak.neutral_mass
        self.extracted_charge = int(
            peak.charge) if override_charge is None else override_charge
        self.extracted_intensity = peak.intensity
        self.extracted_peak = peak

    def default(self, orphan=False):
        '''Populate the extracted attributes of this object from the matching
        original attributes.

        This usually reflects a failure to find an acceptable deconvolution solution,
        and may indicate that there was no peak at the specified location when ``orphan``
        is :const:`True`

        Parameters
        ----------
        orphan: :class:`bool`
            Whether or not to set :attr:`orphan` to :const:`True`, indicating no peak was
            found near :attr:`mz`.
        '''
        if self.charge == ChargeNotProvided:
            warnings.warn(
                "A precursor has been defaulted with an unknown charge state.")
            self.extracted_charge = ChargeNotProvided
            self.extracted_neutral_mass = neutral_mass(
                self.mz, DEFAULT_CHARGE_WHEN_NOT_RESOLVED)
            self.extracted_intensity = self.intensity
            self.defaulted = True
        else:
            self.extracted_charge = int(self.charge)
            self.extracted_neutral_mass = self.neutral_mass
            self.extracted_intensity = self.intensity
            self.defaulted = True
        if orphan:
            self.orphan = True

    @property
    def neutral_mass(self):
        if self.charge == ChargeNotProvided:
            warnings.warn(
                "A precursor with an unknown charge state was used to compute a neutral mass.")
            return neutral_mass(self.mz, DEFAULT_CHARGE_WHEN_NOT_RESOLVED)
        return neutral_mass(self.mz, self.charge)

    @property
    def extracted_mz(self):
        if self.extracted_charge == ChargeNotProvided or (
                self.extracted_charge == 0 and self.charge == ChargeNotProvided):
            warnings.warn(
                "A precursor with an unknown charge state was used to compute a m/z.")
            return mass_charge_ratio(self.mz, DEFAULT_CHARGE_WHEN_NOT_RESOLVED)
        return mass_charge_ratio(self.extracted_neutral_mass, self.extracted_charge)

    @property
    def precursor(self):
        if self.precursor_scan_id is None:
            return None
        return self.source.get_scan_by_id(self.precursor_scan_id)

    @property
    def product(self):
        if self.product_scan_id is None:
            return None
        return self.source.get_scan_by_id(self.product_scan_id)

    def copy(self):
        dup = self.__class__(
            self.mz, self.intensity, self.charge, self.precursor_scan_id, self.source,
            self.extracted_neutral_mass, self.extracted_charge, self.extracted_intensity,
            self.peak, self.extracted_peak, self.defaulted, self.orphan,
            self.product_scan_id, self.annotations, self.coisolation)
        return dup

    def clone(self):
        return self.copy()