from ms_deisotope.feature_map.feature_map import DeconvolutedLCMSFeatureForest
from ms_deisotope.data_source.metadata.scan_traits import IsolationWindow

class IsolationWindowDemultiplexingDeconvolutedLCMSFeatureForestCollection(object):
    '''Aggregate MS2 features per isolation window over the course of a
    deconvoluted LC-MS/MS run.

    Attributes
    ----------
    reader : ProcessedMzMLDeserializer
        A reader over a deconvoluted LC-MS/MS run
    minimum_mass : float
        The mass below which MS2 peaks will be ignored
    minimum_intensity : float
        The intensity below which MS2 peaks will be ignored
    error_tolerance : float
        The mass error tolerance (PPM) when aggregating peaks into features
    delta_rt : float
        The maximum gap size to tolerate in the time dimension

    '''

    def __init__(self, reader, minimum_mass=80.0, minimum_intensity=10.0, error_tolerance=1e-5, delta_rt=1.0):
        self.reader = reader
        self.error_tolerance = error_tolerance
        self.minimum_mass = minimum_mass
        self.minimum_intensity = minimum_intensity
        self.delta_rt = delta_rt
        self.forest_for_isolation_window = dict()

    def forest_for(self, isolation_window):
        '''Get the :class:`~.DeconvolutedLCMSFeatureForest` for a specific isolation window,
        or if one does not exist, create it.

        Parameters
        ----------
        isolation_window : :class:`~.IsolationWindow`
            The isolation window to select with

        Returns
        -------
        :class:`~.DeconvolutedLCMSFeatureForest`
        '''
        if isolation_window is None:
            raise ValueError("An isolation window None!")
        try:
            forest = self.forest_for_isolation_window[isolation_window]
        except KeyError:
            forest = DeconvolutedLCMSFeatureForest(
                error_tolerance=self.error_tolerance)
            self.forest_for_isolation_window[isolation_window] = forest
        return forest

    def all_forests_overlapping_window(self, isolation_window):
        '''Get all :class:`~.DeconvolutedLCMSFeatureForest` that are from isolation windows
        which overlap with `isolation_window`.

        The result list is ordered by isolation window lower bound.

        Parameters
        ----------
        isolation_window : :class:`~.IsolationWindow`
            The isolation window to query for overlaps with

        Returns
        -------
        list of (:class:`~.IsolationWindow`, :class:`~.DeconvolutedLCMSFeatureForest`) pairs.
            All feature
        '''
        forests = []
        lb = isolation_window.lower_bound
        ub = isolation_window.upper_bound
        for window, forest in self.forest_for_isolation_window.items():
            if window.spans(lb) or window.spans(ub):
                forests.append((window, forest))
        if len(forests) > 1:
            forests.sort(key=lambda x: x[0].lower_bound)
        return forests

    def handle_scan(self, scan):
        '''Add `scan`'s peaks to the appropriate feature forest by its isolation window.

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The deconvoluted scan to process
        '''
        isolation_window = scan.isolation_window
        forest = self.forest_for(isolation_window)
        time = scan.scan_time
        for peak in scan.deconvoluted_peak_set:
            if peak.neutral_mass >= self.minimum_mass and peak.intensity >= self.minimum_intensity:
                forest.handle_peak(peak, time)

    def complete(self):
        '''Apply any last post-processing to the feature forests, smoothing
        overlapping features and breaking gaps.
        '''
        for key, forest in self.forest_for_isolation_window.items():
            forest.smooth_overlaps()
            forest.split_sparse(self.delta_rt)

    def aggregate_peaks(self):
        '''Run the aggregation process to completion, and post-process the results,
        producing a fully-fledged set of :class:`~.DeconvolutedLCMSFeatureForest`
        '''
        self.reader.reset()
        self.reader.make_iterator(grouped=False)
        n = len(self.reader)
        for i, scan in enumerate(self.reader):
            if i % 1000 == 0 and i:
                print("Aggregated %d/%d Scans (%0.2f%%)" %
                      (i, n, i * 100.0 / n))
            if scan.ms_level > 1:
                self.handle_scan(scan)
        self.complete()

    @classmethod
    def demultiplex(cls, reader, minimum_mass=80.0, minimum_intensity=10.0, error_tolerance=1e-5, delta_rt=1.0):
        '''Aggregate MS2 features per isolation window over the course of a deconvoluted LC-MS/MS run.

        Parameters
        ----------
        reader : ProcessedMzMLDeserializer
            A reader over a deconvoluted LC-MS/MS run
        minimum_mass : float
            The mass below which MS2 peaks will be ignored
        minimum_intensity : float
            The intensity below which MS2 peaks will be ignored
        error_tolerance : float
            The mass error tolerance (PPM) when aggregating peaks into features
        delta_rt : float
            The maximum gap size to tolerate in the time dimension

        Returns
        -------
        IsolationWindowDemultiplexingDeconvolutedLCMSFeatureForestCollection
        '''
        task = cls(reader, minimum_mass, minimum_intensity,
                   error_tolerance, delta_rt)
        task.aggregate_peaks()
        return task


demultiplex_ms2_features = IsolationWindowDemultiplexingDeconvolutedLCMSFeatureForestCollection.demultiplex
