'''Represent the basic structures of a mass spectrum and its processed contents,
and provide an interface for manipulating that data.
'''
import warnings
try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence

import numpy as np

from ms_peak_picker import (
    pick_peaks, reprofile, average_signal,
    scan_filter, PeakIndex, PeakSet)

from ms_deisotope.utils import decimal_shift
from ms_deisotope.deconvolution import deconvolute_peaks

from ms_deisotope.data_source.metadata.scan_traits import (
    IsolationWindow,
    ScanAcquisitionInformation)

from ms_deisotope.data_source.metadata.activation import ActivationInformation
from ms_deisotope.data_source.metadata.instrument_components import InstrumentInformation


from .base import (ScanBase, RawDataArrays, PrecursorInformation)


class Scan(ScanBase):
    """Container for mass spectral data and associated descriptive information.

    A :class:`Scan` object is a generic object intended to be created by a :class:`ScanDataSource` and describes
    a mass spectrum at each level of processing (Profile --> Peak Fitted --> Deconvoluted). The raw object
    provided by the source is wrapped and queried lazily when an attribute is requested, delegated through
    :attr:`source`.

    To access the raw mass spectral signal, see :attr:`arrays`. This can be converted to a centroided peak
    list using :meth:`pick_peaks`, populating :attr:`peak_set`. After picking peaks, :meth:`deconvolute`
    can be used to deconvolute the spectrum, populating :attr:`deconvoluted_peak_set`.

    To apply signal processing filters to the spectrum, :meth:`reprofile` can create a copy of the object
    replacing :attr:`arrays` with a synthetic profile created from the fitted peaks in :attr:`peak_set`.
    Multiple spectra may be averaged together to form a new synthetic profile spectrum using :meth:`average_with`
    (when the scans are already chosen) or :meth:`average` (when selecting nearby scans automatically, MS1-only).
    The :meth:`denoise` method can remove local unstructured noise. :meth:`transform` can apply any filter
    from :mod:`ms_peak_picker.scan_filter`, returning a copy.

    Attributes
    ----------
    deconvoluted_peak_set : :class:`ms_deisotope.DeconvolutedPeakSet` or None
        Deconvoluted peaks resulting from charge state deconvolution and deisotoping. Will
        be `None` if deconvolution has not been done.
    peak_set : :class:`ms_peak_picker.PeakSet` or None
        Picked peaks and (possibly) associated raw data points as produced by :meth:`pick_peaks`.
        Will be `None` if peak picking has not been done.
    product_scans : list
        A list of :class:`Scan` instances which were produced by fragmenting ions from this one.
        This attribute is not guaranteed to be populated depending upon how the scan is loaded.
    source : :class:`ScanDataSource`
        The object which produced this scan and which defines the methods for retrieving common
        attributes from the underlying data structures.
    precursor_information: :class:`PrecursorInformation` or None
        Descriptive metadata for the ion which was chosen for fragmentation, and a reference to
        the precursor scan
    arrays: :class:`RawDataArrays`
        A pair of :class:`numpy.ndarray` objects corresponding to the raw m/z and intensity data points
    id: str
        The unique identifier for this scan as given by the source
    title: str
        The human-readable display string for this scan as shown in some external software
    ms_level: int
        The degree of fragmentation performed. 1 corresponds to a MS1 or "Survey" scan, 2 corresponds
        to MS/MS, and so on. If :attr:`ms_level` > 1, the scan is considered a "tandem scan" or "MS^n" scan
    scan_time: float
        The time the scan was acquired during data acquisition. The unit of time will always be minutes.
    drift_time: float or None
        The time measured by the ion mobility spectrometer for this scan or frame. This quantity is None
        if the scan does not have ion mobility information associated with it, which is usually recorded
        in :attr:`acquisition_information`
    index: int
        The integer number indicating how many scans were acquired prior to this scan.
    is_profile: bool
        Whether this scan's raw data points corresponds to a profile scan or whether the raw data was
        pre-centroided.
    polarity: int
        If the scan was acquired in positive mode, the value ``+1``.  If the scan was acquired in negative
        mode, the value ``-1``. May be used to indicating how to calibrate charge state determination methods.
    activation: :class:`.ActivationInformation` or None
        If this scan is an MS^n scan, this attribute will contain information about the process
        used to produce it from its parent ion.
    instrument_configuration: :class:`~.InstrumentInformation`
        The instrument configuration used to acquire this scan.
    acquisition_information: :class:`.ScanAcquisitionInformation` or None
        Describes the type of event that produced this scan, as well as the scanning method
        used.
    isolation_window: :class:`.IsolationWindow` or None
        Describes the range of m/z that were isolated from a parent scan to create this scan
    annotations: dict
        A set of key-value pairs describing the scan not part of the standard interface
    """

    def __init__(self, data, source, peak_set=None, deconvoluted_peak_set=None, product_scans=None, annotations=None):
        if product_scans is None:
            product_scans = []
        if annotations is None:
            annotations = dict()
        self.source = source
        self.peak_set = peak_set
        self.deconvoluted_peak_set = deconvoluted_peak_set

        self._data = data

        self._arrays = None
        self._id = None
        self._title = None
        self._ms_level = None
        self._scan_time = None
        self._precursor_information = None
        self._index = None
        self._is_profile = None
        self._polarity = None
        self._activation = None
        self._acquisition_information = None
        self._isolation_window = None
        self._instrument_configuration = None

        self._annotations = None
        self._external_annotations = annotations

        self.product_scans = product_scans

    def clone(self, deep=True):
        """Return a copy of the :class:`Scan` object
        wrapping the same reference data, potentially a deep
        one

        Parameters
        ----------
        deep: :class:`bool`

        Returns
        -------
        :class:`Scan`
        """
        dup = self.__class__(
            self._data, self.source,

            self.peak_set.clone() if self.peak_set is not None else None
            if deep else self.peak_set,

            self.deconvoluted_peak_set.clone() if self.deconvoluted_peak_set is not None else None
            if deep else self.deconvoluted_peak_set,

            [s.clone(deep=deep) for s in self.product_scans], self._external_annotations.copy())
        return dup

    def _load(self):
        _ = self.arrays
        _ = self.id
        _ = self.title
        _ = self.ms_level
        _ = self.scan_time
        _ = self.index
        _ = self.polarity
        _ = self.precursor_information
        _ = self.activation
        _ = self.acquisition_information
        _ = self.isolation_window
        _ = self.is_profile
        _ = self.instrument_configuration
        _ = self.annotations
        _ = None
        del _
        return self

    def _unload(self):
        self._arrays = None
        self._id = None
        self._title = None
        self._ms_level = None
        self._scan_time = None
        self._precursor_information = None
        self._index = None
        self._is_profile = None
        self._polarity = None
        self._activation = None
        self._acquisition_information = None
        self._isolation_window = None
        self._instrument_configuration = None

    def clear(self, full=False):
        """Releases all associated in-memory data and clears the cached
        attributes.

        The data reference attribute :attr:`_data` is retained
        and unchanged.

        Parameters
        ----------
        full: bool
            Whether to clear more attributes to aggressively free memory.
        """
        if self.source is not None:
            self.source._scan_cleared(self)
        self._unload()
        if full:
            self.peak_set = None
            self.deconvoluted_peak_set = None
            self.product_scans = []
            self._external_annotations = {}

    @property
    def ms_level(self):
        '''The degree of fragmentation performed. 1 corresponds to a MS1 or "Survey" scan, 2 corresponds
        to MS/MS, and so on. If :attr:`ms_level` > 1, the scan is considered a "tandem scan" or "MS^n" scan

        Returns
        -------
        :class:`int`
        '''
        if self._ms_level is None:
            self._ms_level = self.source._ms_level(self._data)
        return self._ms_level

    @ms_level.setter
    def ms_level(self, value):
        self._ms_level = int(value)

    @property
    def is_profile(self):
        '''Whether this scan's raw data points corresponds to a profile scan or whether the raw data was
        pre-centroided.

        Returns
        -------
        :class:`bool`
        '''
        if self._is_profile is None:
            self._is_profile = self.source._is_profile(self._data)
        return self._is_profile

    @is_profile.setter
    def is_profile(self, value):
        self._is_profile = bool(value)

    @property
    def polarity(self):
        '''If the scan was acquired in positive mode, the value ``+1``.  If the scan was acquired in negative
        mode, the value ``-1``. May be used to indicating how to calibrate charge state determination methods.

        Returns
        -------
        :class:`int`
        '''
        if self._polarity is None:
            self._polarity = self.source._polarity(self._data)
        return self._polarity

    @polarity.setter
    def polarity(self, value):
        self._polarity = int(value)

    @property
    def scan_time(self):
        '''The time the scan was acquired during data acquisition. The unit of time will always
        be minutes.

        Returns
        -------
        :class:`float`
        '''
        if self._scan_time is None:
            self._scan_time = self.source._scan_time(self._data)
        return self._scan_time

    @scan_time.setter
    def scan_time(self, value):
        self._scan_time = float(value)

    @property
    def arrays(self):
        '''A pair of :class:`numpy.ndarray` objects corresponding to the raw m/z and
        intensity data points.

        These arrays are wrapped in a :class:`~.RawDataArrays` instance, which provides
        additional methods.

        Returns
        -------
        :class:`~.RawDataArrays`
        '''
        if self._arrays is None:
            self._arrays = RawDataArrays(*self.source._scan_arrays(self._data))
        return self._arrays

    @arrays.setter
    def arrays(self, value):
        if isinstance(value, RawDataArrays) or value is None:
            self._arrays = value
        elif isinstance(value, Sequence):
            if len(value) == 2:
                self._arrays = RawDataArrays(*map(np.asanyarray, value))
            elif len(value) == 3:
                self._arrays = RawDataArrays(*map(np.asanyarray, value[:2]), arrays=dict(value[2]))
            else:
                raise ValueError("Too many values to convert. Please provide two arrays, "
                                 "or two arrays and a dictionary of additional arrays.")
        else:
            raise TypeError(
                "arrays must be an instance of RawDataArrays or a pair of numpy arrays")

    @property
    def title(self):
        '''The human-readable display string for this scan as shown in some external software.

        Returns
        -------
        :class:`str`
        '''
        if self._title is None:
            self._title = self.source._scan_title(self._data)
        return self._title

    @title.setter
    def title(self, value):
        self._title = value

    @property
    def id(self):
        '''The within run unique scan identifier.

        Returns
        -------
        :class:`str`
        '''
        if self._id is None:
            self._id = self.source._scan_id(self._data)
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    scan_id = id

    @property
    def index(self):
        '''The integer number indicating how many scans were acquired prior to this scan.

        Returns
        -------
        :class:`int`
        '''
        if self._index is None:
            self._index = self.source._scan_index(self._data)
        return self._index

    @index.setter
    def index(self, value):
        self._index = int(value)

    @property
    def precursor_information(self):
        '''Descriptive metadata for the ion which was chosen for fragmentation, and a reference to
        the precursor scan.

        Returns
        -------
        :class:`~.PrecursorInformation`'''
        if self.ms_level < 2:
            return None
        if self._precursor_information is None:
            self._precursor_information = self.source._precursor_information(
                self._data)
        return self._precursor_information

    @precursor_information.setter
    def precursor_information(self, value):
        if not isinstance(value, PrecursorInformation) and value is not None:
            raise TypeError("precursor_information must be a %r instance" % (
                PrecursorInformation, ))
        self._precursor_information = value

    @property
    def activation(self):
        '''If this scan is an MS^n scan, this attribute will contain information about the process
        used to produce it from its parent ion.

        Returns
        -------
        :class:`~.ActivationInformation`
        '''
        if self.ms_level < 2:
            return None
        if self._activation is None:
            self._activation = self.source._activation(self._data)
        return self._activation

    @activation.setter
    def activation(self, value):
        if not isinstance(value, ActivationInformation) and value is not None:
            raise TypeError(
                "activation must be an %r instance" % (ActivationInformation, ))
        self._activation = value

    @property
    def isolation_window(self):
        '''Describes the range of m/z that were isolated from a parent scan to create this scan.

        Returns
        -------
        :class:`~.IsolationWindow`
        '''
        if self.ms_level < 2:
            return None
        if self._isolation_window is None:
            self._isolation_window = self.source._isolation_window(self._data)
        return self._isolation_window

    @isolation_window.setter
    def isolation_window(self, value):
        if isinstance(value, IsolationWindow) or value is None:
            self._isolation_window = value
        elif isinstance(value, Sequence):
            if len(value) == 2:
                lo, hi = value
                width = (hi - lo) / 2.
                center = lo + width
                self._isolation_window = IsolationWindow(center, width, width)
            elif len(value) == 3:
                lo, center, hi = value
                self._isolation_window = IsolationWindow(lo, center, hi)
            else:
                raise ValueError("Could not convert %r to an %r" %
                                 (value, IsolationWindow))
        else:
            raise TypeError(
                "isolation_window must be an either an %r instance or a sequence of two or three elements" % (
                    IsolationWindow))

    @property
    def acquisition_information(self):
        '''Describes the type of event that produced this scan, as well as the scanning method
        used.'''
        if self._acquisition_information is None:
            self._acquisition_information = self.source._acquisition_information(
                self._data)
        return self._acquisition_information

    @acquisition_information.setter
    def acquisition_information(self, value):
        if not isinstance(value, ScanAcquisitionInformation) and value is not None:
            raise TypeError("acquisition_information must be an instance of %r" % (
                ScanAcquisitionInformation, ))
        self._acquisition_information = value

    @property
    def instrument_configuration(self):
        '''The instrument configuration used to acquire this scan.'''
        if self._instrument_configuration is None:
            self._instrument_configuration = self.source._instrument_configuration(
                self._data)
        return self._instrument_configuration

    @instrument_configuration.setter
    def instrument_configuration(self, value):
        if not isinstance(value, InstrumentInformation) and value is not None:
            raise TypeError("instrument_configuration must be an instance of %r" % (
                InstrumentInformation, ))
        self._instrument_configuration = value

    @property
    def annotations(self):
        '''A set of key-value pairs describing the scan not part of the standard interface'''
        if self._annotations is None:
            self._annotations = self.source._annotations(self._data)
            self._annotations.update(self._external_annotations)
        return self._annotations

    @annotations.setter
    def annotations(self, value):
        self._external_annotations = dict(value)
        self._annotations = self._external_annotations.copy()

    def bind(self, source):
        super(Scan, self).bind(source)
        self.source = source
        return self

    def unbind(self):
        super(Scan, self).unbind()
        self.source = None
        return self

    def __repr__(self):
        try:
            precursor_information = self.precursor_information
        except KeyError:
            precursor_information = None
        return "Scan(%r, index=%d, time=%0.4f, ms_level=%r%s)" % (
            self.id, (self.index if self.index is not None else -1), (
                self.scan_time if self.scan_time is not None else -1), self.ms_level,
            ", " + repr(precursor_information) if precursor_information else '')

    # peak manipulation

    def __iter__(self):
        if self.peak_set is None:
            raise ValueError("Cannot iterate over peaks in a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        return iter(self.peak_set)

    def __getitem__(self, i):
        if self.peak_set is None:
            raise ValueError("Cannot retrieve peaks in a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        return self.peak_set[i]

    def __len__(self):
        if self.peak_set is None:
            raise ValueError("Cannot count peaks in a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        return len(self.peak_set)

    def __bool__(self):
        return bool(self._data)

    def __nonzero__(self):
        return self.__bool__()

    def has_peak(self, *args, **kwargs):
        """A wrapper around :meth:`ms_peak_picker.PeakSet.has_peak` to query the
        :class:`ms_peak_picker.FittedPeak` objects picked for this scan.

        Parameters
        ----------
        mz: float
            The m/z to search for
        error_tolerance: float
            The parts per million mass error tolerance to use

        Returns
        -------
        :class:`ms_peak_picker.FittedPeak` or None
            The peak closest to the query m/z within the error tolerance window or None
            if there are no peaks satisfying the requirements

        Raises
        ------
        ValueError:
            If the scan has not yet had peaks picked yet

        See Also
        --------
        :meth:`.Scan.pick_peaks`
        """
        if self.peak_set is None:
            raise ValueError("Cannot search for peaks in a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        return self.peak_set.has_peak(*args, **kwargs)

    def pick_peaks(self, *args, **kwargs):
        """A wrapper around :func:`ms_peak_picker.pick_peaks` which will populate the
        :attr:`peak_set` attribute of this scan.

        Parameters
        ----------
        fit_type : str, optional
            The name of the peak model to use. One of "quadratic", "gaussian", "lorentzian", or "apex"
        signal_to_noise_threshold : int, optional
            Minimum signal-to-noise measurement to accept a peak
        intensity_threshold : int, optional
            Minimum intensity measurement to accept a peak
        threshold_data : bool, optional
            Whether to apply thresholds to the data
        target_envelopes : list, optional
            A sequence of (start m/z, end m/z) pairs, limiting peak picking to only those intervals
        transforms : list, optional
            A list of :class:`scan_filter.FilterBase` instances or callable that
            accepts (mz_array, intensity_array) and returns (mz_array, intensity_array) or
            `str` matching one of the premade names in `scan_filter.filter_register`
        verbose : bool, optional
            Whether to log extra information while picking peaks
        start_mz : float, optional
            A minimum m/z value to start picking peaks from
        stop_mz : float, optional
            A maximum m/z value to stop picking peaks after
        *args :
            Passed along to :func:`ms_peak_picker.pick_peaks`
        **kwargs :
            Passed along to :func:`ms_peak_picker.pick_peaks`

        Returns
        -------
        Scan
            Returns self
        """
        # Check to see if the user requested one of the ms_peak_picker fits or wanted
        # to use the vendor peak picker if provided.
        fit_type_k = kwargs.get("fit_type")
        if len(args) > 0:
            fit_type_a = args[0]
        else:
            fit_type_a = None
        if fit_type_k == 'vendor' or fit_type_a == 'vendor':
            try:
                peaks = self.source._pick_peaks_vendor(self._data, *args, **kwargs)
                self.peak_set = peaks
                return self
            except NotImplementedError:
                pass
        # Prepare the peak picking parameters
        mzs, intensities = self.arrays
        if len(mzs) == 0:
            self.peak_set = PeakIndex(mzs, intensities, PeakSet([]))
            return self
        if self.is_profile:
            peak_mode = 'profile'
        else:
            peak_mode = 'centroid'

        kwargs.setdefault('peak_mode', peak_mode)

        self.peak_set = pick_peaks(mzs, intensities, *args, **kwargs)
        return self

    def deconvolute(self, *args, **kwargs):
        """A wrapper around :func:`ms_deisotope.deconvolution.deconvolute_peaks`.

        The scan must have had its peaks picked before it can be deconvoluted.

        Parameters
        ----------
        decon_config : dict, optional
            Parameters to use to initialize the deconvoluter instance produced by
            ``deconvoluter_type``
        charge_range : tuple of integers, optional
            The range of charge states to consider.
        error_tolerance : float, optional
            PPM error tolerance to use to match experimental to theoretical peaks
        priority_list : list, optional
            The set of peaks to target for deconvolution to be able to enforce external
            constraints on, such as selected precursors for fragmentation.
        left_search_limit : int, optional
            The maximum number of neutron shifts to search to the left  (decrease) from
            each query peak
        right_search_limit : int, optional
            The maximum number of neutron shifts to search to the right (increase) from
            each query peak
        left_search_limit_for_priorities : int, optional
            The maximum number of neutron shifts to search to the left (decrease) from
            each query peak for priority targets
        right_search_limit_for_priorities : None, optional
            The maximum number of neutron shifts to search to the right (increase) from
            each query peak for priority targets
        charge_carrier : float, optional
            The mass of the charge carrier. Defaults to PROTON
        truncate_after : float, optional
            The percentage of the isotopic pattern to include. Defaults to TRUNCATE_AFTER
        deconvoluter_type : type or callable, optional
            A callable returning a deconvoluter. Defaults to :class:`~.AveraginePeakDependenceGraphDeconvoluter`
        **kwargs
            Additional keywords passed to :func:`~.deconvolute_peaks`

        Returns
        -------
        Scan
            Returns self

        Raises
        ------
        ValueError
            If :attr:`peak_set` is None, a :class:`ValueError` will be raised
            indicating that a scan must be centroided before it can be deconvoluted

        See Also
        --------
        :func:`~.deconvolute_peaks`
        """
        if self.peak_set is None:
            raise ValueError("Cannot deconvolute a scan that has not been "
                             "centroided. Call `pick_peaks` first.")
        charge_range = kwargs.get("charge_range", (1, 8))
        if self.polarity < 0 and max(charge_range) > 0:
            charge_range = tuple(c * self.polarity for c in charge_range)
        kwargs['charge_range'] = charge_range
        decon_results = deconvolute_peaks(self.peak_set, *args, **kwargs)
        self.deconvoluted_peak_set = decon_results.peak_set
        return self

    def pack(self, bind=False):
        '''Pack the (dispersed) representation of the data in this :class:`Scan`
        into a packed :class:`ProcessedScan` object.

        .. note::
            A reference to :attr:`precursor_information` is passed to the returned
            :class:`ProcessedScan`, so both objects share it. Because the :attr:`~.PrecursorInformation.product`
            property works by looking up the scan in :attr:`source`, it's not possible to
            retrieve the :class:`ProcessedScan` this way.

        Parameters
        ----------
        bind: bool
            Whether or not the :class:`ProcessedScan` object should also be bound
            to :attr:`source`

        Returns
        -------
        :class:`ProcessedScan`
        '''
        precursor_info = self.precursor_information
        scan = ProcessedScan(
            self.id, self.title, precursor_info,
            self.ms_level, self.scan_time, self.index,
            self.peak_set.pack() if self.peak_set is not None else None,
            self.deconvoluted_peak_set,
            self.polarity,
            self.activation,
            self.acquisition_information,
            self.isolation_window,
            self.instrument_configuration,
            self.product_scans,
            self.annotations,
            source=self.source if bind else None)
        return scan

    # signal transformation

    def reprofile(self, max_fwhm=0.2, dx=0.005, model_cls=None, override_fwhm=None):
        """Use the picked peaks in :attr:`peak_set` to create a new
        profile mass spectrum using a peak shape model.

        Parameters
        ----------
        max_fwhm : float, optional
            Maximum peak width above which peaks will be ignored
        dx : float, optional
            The distance between each new point in m/z space in the
            reprofiled spectrum
        model_cls : ms_peak_picker.peak_statistics.PeakShapeModel, optional
            The peak shape model to use to generate the profile data from
            the centroided peaks. Defaults a Gaussian model

        Returns
        -------
        Scan
            A shallow copy of this scan with its :attr:`arrays` replaced with
            the new reprofiled arrays

        Raises
        ------
        ValueError
            A scan that has not been centroided and is already in profile mode
            must have its peaks picked before it can be reprofiled.
        """
        if self.peak_set is None and self.is_profile:
            raise ValueError(
                "Cannot reprofile a scan that has not been centroided")
        elif self.peak_set is None and not self.is_profile:
            self.pick_peaks()
        if not self.peak_set:
            arrays = (np.array([], dtype=float), np.array([], dtype=float))
        else:
            arrays = reprofile(self.peak_set, max_fwhm, dx,
                               model_cls, override_fwhm=override_fwhm)
        scan = WrappedScan(
            self._data, self.source, arrays,
            list(self.product_scans), is_profile=True,
            annotations=self._external_annotations)
        return scan

    def denoise(self, scale=5.0, window_length=2.0, region_width=10):
        """Create a shallow copy of the scan with a noise reduction
        transformation applied.

        This method uses the scan filter :class:`ms_peak_picker.scan_filter.FTICRBaselineRemoval`
        which uses the MasSpike noise reduction algorithm.

        Parameters
        ----------
        scale : float, optional
            The multiplier of the local noise window to remove
        window_length : float, optional
            The width (in m/z) of each window
        region_width : int, optional
            The width (in m/z) of each region of windows

        Returns
        -------
        Scan
            The denoised version of this scan
        """
        mzs, intensities = self.arrays
        mzs = mzs.astype(float)
        intensities = intensities.astype(float)
        transform = scan_filter.FTICRBaselineRemoval(
            window_length=window_length, scale=scale, region_width=region_width)
        mzs, intensities = transform(mzs, intensities)
        return WrappedScan(self._data, self.source,
                           (mzs, intensities), list(self.product_scans),
                           is_profile=self.is_profile,
                           annotations=self._external_annotations)

    def transform(self, filters=None):
        '''Applies a series of :class:`ms_peak_picker.scan_filter.FilterBase`,
        or strings that are recognized by :func:`ms_peak_picker.scan_filter.transform`

        Arguments
        ---------
        filters: :class:`Iterable`
            An iterable of transformations of :class:`ms_peak_picker.scan_filter.FilterBase`
            or strings.

        Returns
        -------
        :class:`WrappedScan`
        '''
        mzs, intensities = self.arrays
        mzs = mzs.astype(float)
        intensities = intensities.astype(float)
        mzs, intensities = scan_filter.transform(
            mzs, intensities, filters=filters)
        return WrappedScan(self._data, self.source,
                           (mzs, intensities), list(self.product_scans),
                           is_profile=self.is_profile,
                           annotations=self._external_annotations)

    def average_with(self, scans, dx=None, weight_sigma=None):
        r"""Average together multiple scans' raw data arrays to create a composite intensity
        profile for a common m/z axis.

        Parameters
        ----------
        scans: list:
            A list of :class:`Scan` objects
        dx : float, optional
            The distance between each point in the generated common m/z axis.
        weight_sigma : float, optional
            When this value is not None, scans are weighted according to a
            gaussian distribution with a $\sigma$ equal to this value

        Returns
        -------
        :class:`AveragedScan`
            A shallow copy of this scan with its :attr:`arrays` attribute replaced
            with the averaged array
        """
        if dx is None:
            dx = 0.01
            default_dx = True
        else:
            default_dx = False
        scans = [self] + list(scans)
        arrays = []
        for scan in scans:
            if scan.is_profile:
                if scan.arrays.mz.size == 0:
                    continue
                arrays.append(scan.arrays)
            else:
                scan_arrays = scan.reprofile(dx=dx).arrays
                if scan_arrays.mz.size == 0:
                    continue
                arrays.append(scan_arrays)
        if weight_sigma:
            if weight_sigma == 1:
                weight_sigma = 0.025
            weights = self._compute_smoothing_weights(
                scans, mean=self.scan_time, sigma=weight_sigma)
        else:
            weights = None
        if default_dx and arrays:
            if len(arrays) > 2:
                reference = arrays[len(arrays) // 2 + 1]
            else:
                reference = arrays[0]
            empirical_dx = decimal_shift(2 * np.median(np.diff(reference.mz)))
            dx = min(dx, empirical_dx)
        new_arrays = average_signal(arrays, dx=dx, weights=weights)
        indices = [scan.index for scan in scans]
        return AveragedScan(
            self._data, self.source, new_arrays,
            indices, list(self.product_scans), is_profile=True,
            annotations=self._external_annotations)

    def _get_adjacent_scans(self, index_interval=None, rt_interval=None):
        if index_interval is None and rt_interval is None:
            raise ValueError(
                "One of `index_interval` or `rt_interval` must be provided")
        if self.ms_level > 1:
            raise ValueError("Cannot average MSn scans at this time")
        if not self.source:
            raise ValueError("Can't average an unbound scan")
        before = []
        after = []
        if index_interval is not None:
            before = []
            current_index = self.index
            for _ in range(index_interval):
                next_scan = self.source.find_previous_ms1(current_index)
                if next_scan is None:
                    break
                before.append(next_scan)
                current_index = next_scan.index
            before = before[::-1]
            after = []
            current_index = self.index
            for _ in range(index_interval):
                try:
                    next_scan = self.source.find_next_ms1(current_index)
                except ValueError:
                    break
                if next_scan is None:
                    break
                after.append(next_scan)
                current_index = next_scan.index
        elif rt_interval is not None:
            reference_time = self.scan_time
            before = []
            current_index = self.index
            current_time = self.scan_time
            while abs(reference_time - current_time) < rt_interval and current_index > 0:
                next_scan = self.source.find_previous_ms1(current_index)
                if next_scan is None:
                    break
                before.append(next_scan)
                current_index = next_scan.index
                current_time = next_scan.scan_time

            before = before[::-1]

            after = []
            current_index = self.index
            current_time = self.scan_time
            while abs(reference_time - current_time) < rt_interval and current_index > 0:
                try:
                    next_scan = self.source.find_next_ms1(current_index)
                except ValueError:
                    break
                if next_scan is None:
                    break
                after.append(next_scan)
                current_index = next_scan.index
                current_time = next_scan.scan_time
        else:
            raise ValueError(
                "One of `index_interval` or `rt_interval` must be provided")
        return before, after

    def _compute_smoothing_weights(self, scans, mean, sigma=0.025):
        sigma_sqrd_2 = (2 * sigma ** 2)
        time_array = np.array([s.scan_time for s in scans])
        weights = np.exp((-(time_array - mean) ** 2) / sigma_sqrd_2)
        return weights

    def average(self, index_interval=None, rt_interval=None, dx=None, weight_sigma=None):
        r"""Average together multiple scans' raw data arrays to create a composite intensity
        profile for a common m/z axis.

        Only MS1 scans will be averaged with this method

        Either an absolute number of scans before and after can be specified using
        ``index_interval`` or a time window may be specified using ``rt_interval``.

        Parameters
        ----------
        index_interval : int, optional
            The number of scans preceding and proceding to average with.
        rt_interval : float, optional
            The range of time (in minutes) preceding and proceding to
            look for other scans to average with.
        dx : float, optional
            The distance between each point in the generated common m/z axis.
        weight_sigma : float, optional
            When this value is not None, scans are weighted according to a
            gaussian distribution with a $\sigma$ equal to this value

        Returns
        -------
        :class:`AveragedScan`
            A shallow copy of this scan with its :attr:`arrays` attribute replaced
            with the averaged array
        """
        if dx is None:
            dx = 0.01
            default_dx = True
        else:
            default_dx = False
        before, after = self._get_adjacent_scans(index_interval, rt_interval)
        scans = before + [self] + after
        arrays = []
        for scan in scans:
            if scan.is_profile:
                scan_arrays = scan.arrays
            else:
                scan_arrays = scan.reprofile(dx=dx).arrays
            if len(scan_arrays.mz) > 0:
                arrays.append(scan_arrays)
        if weight_sigma:
            if weight_sigma == 1:
                weight_sigma = 0.025
            weights = self._compute_smoothing_weights(
                scans, mean=self.scan_time, sigma=weight_sigma)
        else:
            weights = None
        if default_dx:
            if len(arrays) > 2:
                reference = arrays[len(arrays) // 2 + 1]
            else:
                reference = arrays[0]
            empirical_dx = decimal_shift(2 * np.median(np.diff(reference.mz)))
            dx = min(dx, empirical_dx)

        new_arrays = average_signal(arrays, dx=dx, weights=weights)
        indices = [scan.index for scan in scans]
        return AveragedScan(
            self._data, self.source, new_arrays,
            indices, list(self.product_scans), is_profile=True,
            annotations=self._external_annotations)


class WrappedScan(Scan):
    '''A wrapper around a :class:`Scan` object with one or more attributes overridden.
    '''
    overridable_keys = [
        "_arrays",
        "_id",
        "_title",
        "_ms_level",
        "_scan_time",
        "_precursor_information",
        "_index",
        "_is_profile",
        "_polarity",
        "_activation",
        "_acquisition_information",
        "_isolation_window",
        "_instrument_configuration"
    ]

    def __init__(self, data, source, array_data, product_scans=None, annotations=None, **overrides):
        super(WrappedScan, self).__init__(
            data, source, peak_set=None,
            deconvoluted_peak_set=None,
            annotations=annotations,
            product_scans=product_scans)
        self._arrays = RawDataArrays(*array_data)
        self._overrides = overrides
        for key, value in overrides.items():
            if not key.startswith("_"):
                key = "_" + key
            if key in self.overridable_keys:
                setattr(self, key, value)
            else:
                warnings.warn("Cannot override attribute %s" % (key,))

    def clone(self, deep=True):
        dup = self.__class__(
            self._data, self.source, self.arrays,
            [s.clone(deep=deep) for s in self.product_scans],
            annotations=self._external_annotations,
            **self._overrides)
        if deep:
            dup.peak_set = self.peak_set.clone() if self.peak_set is not None else None
            dup.deconvoluted_peak_set = self.deconvoluted_peak_set.clone()\
                if self.deconvoluted_peak_set is not None else None
        else:
            dup.peak_set = self.peak_set
            dup.deconvoluted_peak_set = self.deconvoluted_peak_set
        return dup


class AveragedScan(WrappedScan):
    '''An averaged :class:`Scan` object, storing additional information for retrieving
    the spectra that were averaged together.

    Attributes
    ----------
    scan_indices: list
        The :attr:`index` values for all of the scans that were averaged together
    '''

    def __init__(self, data, source, array_data, scan_indices, product_scans=None, annotations=None, **overrides):
        super(AveragedScan, self).__init__(
            data, source, array_data,
            product_scans=product_scans, annotations=annotations, **overrides)
        self.scan_indices = scan_indices

    def clone(self, deep=True):
        dup = self.__class__(
            self._data, self.source, self.arrays,
            self.scan_indices,
            [s.clone(deep=deep) for s in self.product_scans],
            annotations=self._external_annotations,
            **self._overrides)
        if deep:
            dup.peak_set = self.peak_set.clone() if self.peak_set is not None else None
            dup.deconvoluted_peak_set = self.deconvoluted_peak_set.clone()\
                if self.deconvoluted_peak_set is not None else None
        else:
            dup.peak_set = self.peak_set
            dup.deconvoluted_peak_set = self.deconvoluted_peak_set
        return dup


class ProcessedScan(ScanBase):
    """Container for mass spectral data and associated descriptive information that has been
    processed and reduced.

    A :class:`ProcessedScan` object has been processed and compacted. It does not carry a
    :class:`~.RawDataArrays` :attr:`arrays` attribute, only at least one of :attr:`peak_set`
    and :attr:`deconvoluted_peak_set`. The rest of its attributes are all loaded in memory
    directly without needing to reconstruct them from a specific format lazily through a
    :class:`~.ScanDataSource`.

    Attributes
    ----------
    deconvoluted_peak_set : :class:`~.DeconvolutedPeakSet` or None
        Deconvoluted peaks resulting from charge state deconvolution and deisotoping. Will
        be `None` if deconvolution has not been done.
    peak_set : :class:`ms_peak_picker.PeakSet` or None
        Picked peaks and (possibly) associated raw data points as produced by :meth:`pick_peaks`.
        Will be `None` if peak picking has not been done.
    product_scans : list
        A list of :class:`Scan` instances which were produced by fragmenting ions from this one.
        This attribute is not guaranteed to be populated depending upon how the scan is loaded.
    source : :class:`ScanDataSource`
        The object which produced this scan and which defines the methods for retrieving common
        attributes from the underlying data structures.
    precursor_information: :class:`PrecursorInformation` or None
        Descriptive metadata for the ion which was chosen for fragmentation, and a reference to
        the precursor scan
    id: str
        The unique identifier for this scan as given by the source
    title: str
        The human-readable display string for this scan as shown in some external software
    ms_level: int
        The degree of fragmentation performed. 1 corresponds to a MS1 or "Survey" scan, 2 corresponds
        to MS/MS, and so on. If :attr:`ms_level` > 1, the scan is considered a "tandem scan" or "MS^n" scan
    scan_time: float
        The time the scan was acquired during data acquisition. The unit of time will always be minutes.
    drift_time: float or None
        The time measured by the ion mobility spectrometer for this scan or frame. This quantity is None
        if the scan does not have ion mobility information associated with it, which is usually recorded
        in :attr:`acquisition_information`
    index: int
        The integer number indicating how many scans were acquired prior to this scan.
    is_profile: bool
        Whether this scan's raw data points corresponds to a profile scan or whether the raw data was
        pre-centroided.
    polarity: int
        If the scan was acquired in positive mode, the value ``+1``.  If the scan was acquired in negative
        mode, the value ``-1``. May be used to indicating how to calibrate charge state determination methods.
    activation: :class:`.ActivationInformation` or None
        If this scan is an MS^n scan, this attribute will contain information about the process
        used to produce it from its parent ion.
    instrument_configuration: :class:`~.InstrumentInformation`
        The instrument configuration used to acquire this scan.
    acquisition_information: :class:`.ScanAcquisitionInformation` or None
        Describes the type of event that produced this scan, as well as the scanning method
        used.
    isolation_window: :class:`.IsolationWindow` or None
        Describes the range of m/z that were isolated from a parent scan to create this scan
    annotations: dict
        A set of key-value pairs describing the scan not part of the standard interface
    """

    def __init__(self, id, title, precursor_information,
                 ms_level, scan_time, index, peak_set,
                 deconvoluted_peak_set, polarity=None, activation=None,
                 acquisition_information=None, isolation_window=None,
                 instrument_configuration=None, product_scans=None,
                 annotations=None, source=None):
        if product_scans is None:
            product_scans = []
        if annotations is None:
            annotations = {}
        self.id = id
        self.title = title
        self.precursor_information = precursor_information
        self.ms_level = ms_level
        self.scan_time = scan_time
        self.index = index
        self.peak_set = peak_set
        self.deconvoluted_peak_set = deconvoluted_peak_set
        self.polarity = polarity
        self.activation = activation
        self.acquisition_information = acquisition_information
        self.isolation_window = isolation_window
        self.instrument_configuration = instrument_configuration
        self.product_scans = product_scans
        self.annotations = annotations
        self.source = source

    def clear(self, full=False):
        '''Clear storage-heavy attribute values

        Parameters
        ----------
        full: bool
            Whether to clear attributes more aggressively to free up space.
        '''
        self.peak_set = None
        self.deconvoluted_peak_set = None
        self.activation = None
        self.acquisition_information = None
        self.isolation_window = None
        self.instrument_configuration = None
        self.product_scans = None

    @property
    def scan_id(self):
        return self.id

    @property
    def is_profile(self):
        '''Whether this scan's raw data points corresponds to a profile scan or whether the raw data was
        pre-centroided.
        '''
        return False

    def _resolve_peaks(self):
        if self.deconvoluted_peak_set is not None:
            return self.deconvoluted_peak_set
        elif self.peak_set is not None:
            return self.peak_set
        else:
            return []

    def __iter__(self):
        return iter(self._resolve_peaks())

    def __getitem__(self, index):
        return self._resolve_peaks()[index]

    def __len__(self):
        return len(self._resolve_peaks())

    def has_peak(self, mass, error_tolerance=2e-5):
        """A wrapper around :meth:`~.DeconvolutedPeakSet.has_peak` to query the
        :class:`~.DeconvolutedPeak` objects picked for this scan. If no deconvoluted
        peaks are available, but centroided peaks are, this method will instead
        behave like :class:`Scan.has_peak`

        Parameters
        ----------
        mass: float
            The mass to search for
        error_tolerance: float
            The parts per million mass error tolerance to use

        Returns
        -------
        :class:`~.PeakBase`
        """
        if self.deconvoluted_peak_set is not None:
            return self.deconvoluted_peak_set.has_peak(mass, error_tolerance)
        elif self.peak_set is not None:
            return self.peak_set.has_peak(mass, error_tolerance)
        else:
            raise ValueError("No peaks available")

    def __repr__(self):
        if self.deconvoluted_peak_set is not None:
            peaks = self.deconvoluted_peak_set
        elif self.peak_set is not None:
            peaks = self.peak_set
        else:
            peaks = []

        pinfo = self.precursor_information
        if pinfo:
            pinfo_string = ", %s" % pinfo
        else:
            pinfo_string = ""

        return "ProcessedScan(id=%s, ms_level=%d, %d peaks%s)" % (
            self.id, self.ms_level, len(peaks), pinfo_string)

    def bind(self, source):
        super(ProcessedScan, self).bind(source)
        self.source = source
        return self

    def unbind(self):
        super(ProcessedScan, self).unbind()
        self.source = None
        return self

    def clone(self, deep=True):
        """Return a copy of the :class:`ProcessedScan` object, potentially a deep
        one

        Parameters
        ----------
        deep: :class:`bool`

        Returns
        -------
        :class:`Scan`
        """
        dup = self.__class__(
            self.id, self.title, self.precursor_information, self.ms_level,
            self.scan_time, self.index,

            self.peak_set.clone() if self.peak_set is not None else None
            if deep else self.peak_set,

            self.deconvoluted_peak_set.clone() if self.deconvoluted_peak_set is not None else None
            if deep else self.deconvoluted_peak_set,

            self.polarity, self.activation, self.acquisition_information,
            self.isolation_window, self.instrument_configuration,
            list(self.product_scans), self.annotations.copy())
        return dup
