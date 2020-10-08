from ms_deisotope.peak_set import DeconvolutedPeakSet
from ms_deisotope.spectrum_graph import (PathFinder, amino_acids, MassWrapper)

from brainpy import SimpleComposition


H2O = SimpleComposition({"H": 2, "O": 1})
NH3 = SimpleComposition({"H": 3, "N": 1})

default_losses = [MassWrapper("-H2O", -H2O.mass()), MassWrapper("-NH3", -NH3.mass())]


class GraphBasedDenoiser(object):
    def __init__(self, components=None, precursor_error_tolerance=2e-5, product_error_tolerance=2e-5,
                 precursor_losses=None, product_losses=None):
        if components is None:
            components = amino_acids
        if precursor_losses is None:
            precursor_losses = default_losses[:]
        if product_losses is None:
            product_losses = default_losses[:]
        self.precursor_error_tolerance = precursor_error_tolerance
        self.product_error_tolerance = product_error_tolerance
        self.path_finder = PathFinder(components, product_error_tolerance)
        self.components = self.path_finder.components
        self.precursor_losses = list(precursor_losses)
        self.product_losses = list(product_losses)

    def initialize_graph(self, peak_set, precursor_mass):
        """Construct a spectrum graph which connects peaks by :attr:`component` masses

        This method relies on :attr:`path_finder` for most of its heavy lifting through
        :meth:`~.PathFinder.find_edges` and :meth:`~.PathFinder.find_complements`.

        Parameters
        ----------
        peak_set : :class:`~.DeconvolutedPeakSet`
            The peak set to search
        precursor_mass : float
            The intact mass of the precursor ion

        Returns
        -------
        :class:`~.SpectrumGraph`
        """
        graph = self.path_finder.find_edges(peak_set)
        p = MassWrapper("Precursor", precursor_mass)
        self.path_finder.find_complements(peak_set, graph, p)
        return graph

    def find_precursor_peak(self, peak_set, precursor_mass, losses=None):
        """Find peaks derived from intact- and neutral losses of the precursor ion

        Parameters
        ----------
        peak_set : :class:`~.DeconvolutedPeakSet`
            The peak set to search
        precursor_mass : float
            The intact mass of the precursor ion
        losses : :class:`~.Iterable` of :class:`~.MassWrapper`, optional
            An iterable of losses expressed as :class:`~.MassWrapper` instances. If not
            provided, -H2O and -NH3 will be assumed.

        Returns
        -------
        :class:`set`:
            The indices of the peaks derived from the intact precursor to retain
        """
        if losses is None:
            losses = self.precursor_losses
        peaks = peak_set.all_peaks_for(precursor_mass, self.product_error_tolerance)
        keep = set()
        for peak in peaks:
            keep.add(peak.index.neutral_mass)
        for loss in losses:
            peaks = peak_set.all_peaks_for(
                precursor_mass + loss.mass, self.product_error_tolerance)
            for peak in peaks:
                keep.add(peak.index.neutral_mass)
        return keep

    def find_loss_peaks(self, peak_set, kept_peak_indices, losses=None):
        """Find neutral loss peaks from the set of peaks that are already set to be retained.

        Parameters
        ----------
        peak_set : :class:`~.DeconvolutedPeakSet`
            The peak set to search
        kept_peak_indices : :class:`set`
            The indices of the peaks already retained to find losses for.
        losses : :class:`~.Iterable` of :class:`~.MassWrapper`, optional
            An iterable of losses expressed as :class:`~.MassWrapper` instances. If not
            provided, -H2O and -NH3 will be assumed.

        Returns
        -------
        :class:`set`:
            `kept_peak_indices` is updated and returned
        """
        if losses is None:
            losses = self.product_losses
        for index in tuple(kept_peak_indices):
            peak = peak_set[index]
            for loss in losses:
                for loss_peak in peak_set.all_peaks_for(
                        peak.neutral_mass + loss.mass, self.product_error_tolerance):
                    kept_peak_indices.add(loss_peak.index.neutral_mass)
        return kept_peak_indices

    def walk_graph(self, graph):
        """Traverse `graph` along its edges, and save the peak indices at either end of each
        edge for use later.

        [description]

        Parameters
        ----------
        graph : :class:`~.SpectrumGraph`
            The graph built by :meth:`initialize_graph`

        Returns
        -------
        :class:`set`:
            The set of peaks which were joined by edges.
        """
        keep = set()
        for e in graph.transitions:
            keep.add(e.start.index)
            keep.add(e.end.index)
        return keep

    def select_peaks(self, peak_set, indices):
        """Given a :class:`~.DeconvolutedPeakSet` and a set of indices into it, construct a
        new :class:`~.DeconvolutedPeakSet` with copies of the peaks at those indices

        Parameters
        ----------
        peak_set : :class:`~.DeconvolutedPeakSet`
            The peak set to search
        indices : :class:`set`
            The indices of the peaks to retain

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`
        """
        keep = []
        for i in indices:
            keep.append(peak_set[i].clone())
        keep = DeconvolutedPeakSet(keep)
        keep.reindex()
        return keep

    def apply_filter(self, peak_set, precursor_mass):
        """Does the actual work of filtering the peak set.

        This method is designed to decouple the :class:`~.ScanBase`
        consuming interface from the actual filtering process.

        Parameters
        ----------
        peak_set : :class:`~.DeconvolutedPeakSet`
            The peak set to search
        precursor_mass : float
            The intact mass of the precursor ion

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`

        See Also
        --------
        :meth:`filter`
        :meth:`__call__`
        """
        graph = self.initialize_graph(peak_set, precursor_mass)
        indices = self.walk_graph(graph)
        self.find_loss_peaks(peak_set, indices)
        indices.update(self.find_precursor_peak(peak_set, precursor_mass))
        peak_set = self.select_peaks(peak_set, indices)
        return peak_set

    def filter(self, scan):
        """Filter an MSn scan's deconvoluted peak set according to graph-connectedness.

        This noise filter will remove any peak that is not connected to another peak by a
        mass shift from :attr:`components` or connected to such a peak via a loss mass shift
        given by :attr:`product_losses` or :attr:`precursor_losses`

        This method is also used via the :meth:`__call__` slot.

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The MSn scan to filter.

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`

        See Also
        --------
        :meth:`apply_filter`
        """
        precursor_mass = scan.precursor_information.extracted_neutral_mass
        if precursor_mass == 0:
            precursor_mass = scan.precursor_information.neutral_mass
        peak_set = scan.deconvoluted_peak_set
        if peak_set is None:
            raise ValueError("scan.deconvoluted_peak_set cannot be None!")
        filtered_peak_set = self.apply_filter(peak_set, precursor_mass)
        return filtered_peak_set

    def __call__(self, scan):
        """Filter an MSn scan's deconvoluted peak set according to graph-connectedness.

        This is a wrapper around the :meth:`filter` method.

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The MSn scan to filter.

        Returns
        -------
        :class:`~.DeconvolutedPeakSet`

        See Also
        --------
        :meth:`filter`
        """
        return self.filter(scan)
