'''Implements algorithms for clustering scans by
spectral similarity.
'''
import operator
import functools
import bisect
import json

from collections import deque

import numpy as np
from scipy.special import comb

from ms_deisotope.task.log_utils import LogUtilsMixin
from ms_deisotope.data_source.dispatch import (
    SubsequenceMappingProxy as _SubsequenceMappingProxy,
    DynamicallyLoadingProxyResolver as _DynamicallyLoadingResolver)

from .similarity_methods import peak_set_similarity

try:
    basestring
except NameError:
    from six import string_types as basestring

peak_set_getter = operator.attrgetter("peak_set")
deconvoluted_peak_set_getter = operator.attrgetter("deconvoluted_peak_set")


def _ppm_error(x, y):
    return (x - y) / y


@functools.total_ordering
class SpectrumCluster(object):
    '''A collection of similar spectra.

    A :class:`SpectrumCluster` is compare-able by :attr:`neutral_mass`
    and acts as a :class:`collections.Sequence` over its members stored
    in :attr:`scans`.

    Attributes
    ----------
    scans: :class:`list`
        A list of spectra which are similar to each other.
    neutral_mass: :class:`float`
        The neutral mass of the representative spectrum's precursor.
    '''
    def __init__(self, scans=None, neutral_mass=None, average_similarity=None, annotations=None):
        if scans is None:
            scans = []
        if annotations is None:
            annotations = dict()
        self.scans = scans
        self.neutral_mass = None
        if neutral_mass:
            self.neutral_mass = neutral_mass
        elif self.scans:
            self.neutral_mass = self.scans[0].precursor_information.neutral_mass
        else:
            self.neutral_mass = 0.0
        if average_similarity is None and len(self) == 1:
            average_similarity = 1.0
        self._average_similarity = average_similarity
        self.annotations = annotations

    def __lt__(self, other):
        return self.neutral_mass < other.neutral_mass

    def __gt__(self, other):
        return self.neutral_mass > other.neutral_mass

    def __eq__(self, other):
        return (abs(self.neutral_mass - other.neutral_mass) / other.neutral_mass) < 1e-6

    def __repr__(self):
        return "SpectrumCluster(%f, %d)" % (self.neutral_mass, len(self))

    def __iter__(self):
        return iter(self.scans)

    def __len__(self):
        return len(self.scans)

    def __getitem__(self, i):
        return self.scans[i]

    def _invalidate(self):
        self._average_similarity = None

    def append(self, item, incremental_similarity=False):
        '''Add a new spectrum to the cluster.

        Parameters
        ----------
        item: :class:`~.ScanBase`
        '''
        if not self.neutral_mass:
            self.neutral_mass = item.precursor_information.neutral_mass
        if incremental_similarity:
            self._incremental_similarity(item)
        else:
            self._invalidate()
        self.scans.append(item)

    def _calculate_similarity_with(self, scan, from_ix=0, to_ix=None, *args, **kwargs):
        if from_ix is None:
            from_ix = 0
        if to_ix is None:
            to_ix = len(self)
        acc = []
        for member in self[from_ix:to_ix]:
            acc.append(peak_set_similarity(scan, member, *args, **kwargs))
        return acc

    def _incremental_similarity(self, scan, *args, **kwargs):
        new_sims = self._calculate_similarity_with(scan, *args, **kwargs)
        aggregate_size = comb(len(self), 2)
        n = (aggregate_size + len(new_sims))
        if n == 0:
            n = 1
        self._average_similarity = (aggregate_size * self.average_similarity() + sum(new_sims)
                                   ) / n

    def _full_similarity(self, *args, **kwargs):
        similarity_method = kwargs.pop(
            "similarity_method", peak_set_similarity)
        ratings = []
        n = len(self)
        for i in range(n):
            scan_i = self[i]
            for j in range(i + 1, n):
                scan_j = self[j]
                ratings.append(similarity_method(
                    scan_i, scan_j, *args, **kwargs))
        self._average_similarity = sum(ratings) / len(ratings)

    def average_similarity(self, *args, **kwargs):
        '''Calculate the within-cluster similarity among all cluster members
        and returns the average.

        All arguments are forwarded to :func:`~.peak_set_similarity`.

        If the cluster is a singleton or smaller, by definition its similarity
        is ``1.0``.

        Parameters
        ----------
        similarity_method: Callable, optional
            The peak set similarity method to use. Defaults to :func:`~.peak_set_similarity`
            All other arguments are forwarded to this function.

        Returns
        -------
        :class:`float`
        '''
        n = len(self)
        if n <= 1:
            return 1.0
        if self._average_similarity is not None:
            return self._average_similarity
        self._full_similarity(*args, **kwargs)
        return self._average_similarity

    def similarity_matrix(self, *args, **kwargs):
        '''Compute a symmetric similarity matrix comparing each
        spectrum within the cluster to each other spectrum.

        Parameters
        ----------
        similarity_method: Callable, optional
            The peak set similarity method to use. Defaults to :func:`~.peak_set_similarity`
            All other arguments are forwarded to this function.

        Returns
        -------
        :class:`np.ndarray`: m
            The value at position ``m[i,j]`` is the simiarlity of ``self[i]`` with ``self[j]``
        '''
        similarity_method = kwargs.pop("similarity_method", peak_set_similarity)
        n = len(self)
        mat = np.identity(n)
        for i in range(n):
            scan_i = self[i]
            for j in range(i + 1, n):
                scan_j = self[j]
                mat[i, j] = mat[j, i] = (similarity_method(scan_i, scan_j, *args, **kwargs))
        return mat

    def to_dict(self):
        '''Convert the cluster to a JSON-safe :class:`dict`

        Returns
        -------
        :class:`dict`
        '''
        d = {}
        d['neutral_mass'] = self.neutral_mass
        d['size'] = len(self)
        d['average_similarity'] = self.average_similarity()
        scans = []
        for scan in self:
            scan_source = scan.source
            if scan_source is not None:
                source_name = scan_source.source_file
                if not isinstance(source_name, basestring):
                    if hasattr(source_name, 'name'):
                        source_name = source_name.name
            else:
                source_name = ":detatched:"
            scans.append({
                "id": scan.id,
                "source": source_name,
                "neutral_mass": scan.precursor_information.neutral_mass,
            })
        d['scans'] = scans
        return d

    def split_on_charge(self):
        scans = sorted(self, key=lambda x: x.tic(), reverse=True)
        index = {}
        for scan in scans:
            key = scan.precursor_information.charge
            try:
                index[key].append(scan)
            except KeyError:
                index[key] = SpectrumCluster([scan], scan.precursor_information.neutral_mass)
        return index


class SpectrumClusterCollection(object):
    '''A sorted :class:`~.Sequence` of :class:`SpectrumCluster` instances
    that supports searching by precursor mass.
    '''
    def __init__(self, clusters=None):
        if clusters is None:
            clusters = []
        self.clusters = list(clusters)

    def add(self, cluster):
        '''Add a new :class:`SpectrumCluster` to the collection,
        preserving sorted order.

        Parameters
        ----------
        cluster: :class:`SpectrumCluster`
            The cluster to add
        '''
        bisect.insort(self.clusters, cluster)

    def __getitem__(self, i):
        return self.clusters[i]

    def __setitem__(self, i, value):
        self.clusters[i] = value

    def __len__(self):
        return len(self.clusters)

    def __iter__(self):
        return iter(self.clusters)

    def __repr__(self):
        template = "{self.__class__.__name__}({size})"
        size = len(self)
        return template.format(self=self, size=size)

    def _binary_search(self, mass, error_tolerance=1e-5):
        array = self.clusters
        n = len(array)
        lo = 0
        hi = n

        while hi != lo:
            mid = (hi + lo) // 2
            y = array[mid].neutral_mass
            err = (y - mass) / mass
            if hi - lo == 1:
                best_index = mid
                best_error = abs(err)
                i = mid - 1
                while i >= 0:
                    x = array[i]
                    err = abs((x.neutral_mass - mass) / mass)
                    if err < best_error:
                        best_error = err
                        best_index = i
                    elif err > error_tolerance:
                        break
                    i -= 1
                lo_index = i + 1
                i = mid + 1
                while i < n:
                    x = array[i]
                    err = abs((x.neutral_mass - mass) / mass)
                    if err < best_error:
                        best_error = err
                        best_index = i
                    elif err > error_tolerance:
                        break
                    i += 1
                hi_index = i
                return best_index, lo_index, hi_index
            elif err > 0:
                hi = mid
            else:
                lo = mid
        return 0, 0, 0

    def find(self, mass, error_tolerance=1e-5):
        '''Finds the cluster whose precursor mass is closest to
        ``mass`` within ``error_tolerance`` ppm error.

        Parameters
        ----------
        mass: :class:`float`
            The mass to search for.
        error_tolerance: :class:`float`
            The PPM error tolerance to use. Defaults to 1e-5.

        Returns
        -------
        :class:`SpectrumCluster` or :const:`None`
        '''
        target_ix, _, _ = self._binary_search(mass, error_tolerance)
        target = self[target_ix]
        if abs(target.neutral_mass - mass) / mass > error_tolerance:
            return None
        return target

    def find_all(self, mass, error_tolerance=1e-5):
        '''Finds all clusters whose precursor mass is within ``error_tolerance``
        ppm of ``mass``.

        Parameters
        ----------
        mass: :class:`float`
            The mass to search for.
        error_tolerance: :class:`float`
            The PPM error tolerance to use. Defaults to 1e-5.

        Returns
        -------
        :class:`list`
        '''
        _, lo_ix, hi_ix = self._binary_search(mass, error_tolerance)
        result = [
            target for target in self[lo_ix:hi_ix]
            if abs(target.neutral_mass - mass) / mass <= error_tolerance
        ]
        return result


class _PeakGetterStrategyBase(object):
    def __call__(self, scan):
        return self.peaks(scan)


class _FittedPeakAccessorStrategy(_PeakGetterStrategyBase):
    def peaks(self, scan):
        return scan.peak_set

    def tic(self, scan):
        return scan.tic.centroided()


class _DeconvolutedPeakAccessorStrategy(_PeakGetterStrategyBase):
    def peaks(self, scan):
        return scan.deconvoluted_peak_set

    def tic(self, scan):
        return scan.tic.deconvoluted()


class _DynamicPeakAccessorStrategy(_PeakGetterStrategyBase):

    def _ensure_peak_set(self, scan):
        if scan.peak_set is None:
            scan.pick_peaks()
        return scan.peak_set

    def peaks(self, scan):
        if scan.deconvoluted_peak_set is not None:
            return scan.deconvoluted_peak_set
        else:
            return self._ensure_peak_set(scan)

    def tic(self, scan):
        if scan.deconvoluted_peak_set is not None:
            return scan.tic.deconvoluted()
        else:
            self._ensure_peak_set(scan)
            return scan.tic.centroided()


class ScanClusterBuilder(LogUtilsMixin):
    """Clusters spectra based upon peak pattern similarity

    Attributes
    ----------
    clusters : :class:`SpectrumClusterCollection`
        The clusters built so far
    minimum_similarity : float
        The minimum similarity score needed to consider two spectra
        similar enough to form a cluster
    peak_getter : :class:`Callable`
        A function to call on each spectrum to retrieve the peak list
        to cluster over
    precursor_error_tolerance : float
        The maximum precursor mass error (in PPM) to permit between
        two spectra to consider comparing them
    track_incremental_similarity : bool
        Whether to incrementally update a cluster's similarity when it
        grows.
    """

    @classmethod
    def _guess_peak_getter(cls, getter):
        if getter is None:
            return _DynamicPeakAccessorStrategy()
        if callable(getter):
            return getter
        if isinstance(getter, basestring):
            if getter == "d":
                return _DeconvolutedPeakAccessorStrategy()
            if getter in ('p', 'c'):
                return _FittedPeakAccessorStrategy()
            else:
                raise ValueError("Cannot infer peak getter strategy from %r" % (getter, ))

        raise ValueError(
            "Cannot infer peak set getter strategy from %r" % (getter, ))

    def __init__(self, clusters=None, precursor_error_tolerance=1e-5, minimum_similarity=0.1,
                 peak_getter=None, track_incremental_similarity=True):
        peak_getter = self._guess_peak_getter(peak_getter)
        if clusters is None:
            clusters = []
        self.clusters = SpectrumClusterCollection(clusters)
        self.precursor_error_tolerance = precursor_error_tolerance
        self.minimum_similarity = minimum_similarity
        self.peak_getter = peak_getter
        self.track_incremental_similarity = track_incremental_similarity

    def _binsearch_simple(self, x):
        n = len(self)
        lo = 0
        hi = n

        while hi != lo:
            mid = (hi + lo) // 2
            y = self[mid].neutral_mass
            err = y - x
            # Do refinement in `find_best_cluster_for_scan`
            if hi - lo == 1:
                return mid
            elif err > 0:
                hi = mid
            else:
                lo = mid
        return 0

    def peak_set_similarity(self, scan_i, scan_j):
        '''Calculate the similarity between the peaks of
        ``scan_i`` and ``scan_j``, where the peaks are
        those points given by :attr:`peak_getter`

        Parameters
        ----------
        scan_i: :class`~.ScanBase`
        scan_j: :class`~.ScanBase`

        Returns
        -------
        :class:`float`
        '''
        peak_set_a = self.peak_getter(scan_i)
        peak_set_b = self.peak_getter(scan_j)
        return peak_set_similarity(
            peak_set_a, peak_set_b)

    def find_best_cluster_for_scan(self, scan):
        '''Locate the best cluster to add ``scan`` to according to
        precursor mass and peak set similarity.

        Parameters
        ----------
        scan: :class:`~.ScanBase`

        Returns
        -------
        :class:`SpectrumCluster`
        '''
        best_cluster = None
        best_similarity = 0.0
        n = len(self.clusters)
        if n == 0:
            return best_cluster

        center_i = self._binsearch_simple(scan.precursor_information.neutral_mass)
        i = center_i

        while i >= 0:
            cluster = self.clusters[i]
            if abs(_ppm_error(scan.precursor_information.neutral_mass,
                              cluster.neutral_mass)) > self.precursor_error_tolerance:
                break
            similarity = self.peak_set_similarity(scan, cluster[0])
            i -= 1
            if similarity > best_similarity and similarity > self.minimum_similarity:
                best_similarity = similarity
                best_cluster = cluster

        i = center_i + 1
        while i < n:
            cluster = self.clusters[i]
            if abs(_ppm_error(scan.precursor_information.neutral_mass,
                              cluster.neutral_mass)) > self.precursor_error_tolerance:
                break
            similarity = self.peak_set_similarity(scan, cluster[0])
            i += 1
            if similarity > best_similarity and similarity > self.minimum_similarity:
                best_similarity = similarity
                best_cluster = cluster
        return best_cluster

    def add_scan(self, scan):
        '''Add ``scan`` to the cluster collection, adding it to the best
        matching cluster, or starting a new cluster around it if no good
        match can be found.

        Parameters
        ----------
        scan: :class:`~.ScanBase`
        '''
        best_cluster = self.find_best_cluster_for_scan(scan)
        if best_cluster:
            best_cluster.append(scan, incremental_similarity=self.track_incremental_similarity)
        else:
            self.clusters.add(SpectrumCluster([scan]))

    def __iter__(self):
        return iter(self.clusters)

    def __len__(self):
        return len(self.clusters)

    def __getitem__(self, i):
        return self.clusters[i]

    def _get_tic(self, scan):
        try:
            return self.peak_getter.tic(scan)
        except AttributeError:
            return sum(p.intensity for p in self.peak_getter(scan))

    @classmethod
    def cluster_scans(cls, scans, precursor_error_tolerance=1e-5, minimum_similarity=0.1,
                      peak_getter=None, sort=True, track_incremental_similarity=False):
        '''Cluster scans by precursor mass and peak set similarity.

        Parameters
        ----------
        scans: :class:`Iterable`
            An iterable of :class:`Scan`-like objects
        precursor_error_tolerance: :class:`float`
            The PPM mass accuracy threshold for precursor mass differences to
            tolerate when deciding whether to compare spectra. Defaults to 1e-5.
        minimum_similarity: :class:`float`
            The minimum peak set similarity required to consider adding a spectrum
            to a cluster. Defaults to 0.1
        peak_getter: :class:`Callable`
            A callable object used to get peaks from elements of ``scans``.
        sort: :class:`bool`
            Whether or not to sort spectra by their total ion current before clustering.
        '''
        self = cls([], precursor_error_tolerance, minimum_similarity,
                   peak_getter, track_incremental_similarity)
        if sort:
            scans = self._sort_by_tic(scans)
        if len(scans) > 100:
            self.log("Clustering %d Scans" % (len(scans), ))
        n = len(scans)
        report_interval = max(min(n // 10, 1000), 50)
        for i, scan in enumerate(scans):
            if i % report_interval == 0 and i:
                self.log("... Handled %d Scans (%0.2f%%)" % (i, i * 100.0 / n))
            self.add_scan(scan)
        return self.clusters

    def _sort_by_tic(self, scans):
        should_log = len(scans) > 100
        if should_log:
            self.log("Sorting Scans By TIC")
        augmented = []
        n = len(scans)
        for i, scan in enumerate(scans):
            if i % 1000 == 0 and i > 0:
                self.log("... Loaded TIC for %d Scans (%0.2f%%)" % (i, i * 100.0 / n))
            augmented.append((self._get_tic(scan), scan))
        augmented.sort(key=lambda x: x[0], reverse=True)
        scans = [a[1] for a in augmented]
        return scans

    @classmethod
    def iterative_clustering(cls, scans, precursor_error_tolerance=1e-5, similarity_thresholds=None,
                             peak_getter=None):
        '''Cluster scans by precursor mass and peak set similarity, iteratively refining
        clusters with increasing similarity threshold requirements.

        Parameters
        ----------
        scans: :class:`Iterable`
            An iterable of :class:`Scan`-like objects
        precursor_error_tolerance: :class:`float`
            The PPM mass accuracy threshold for precursor mass differences to
            tolerate when deciding whether to compare spectra. Defaults to 1e-5.
        similarity_thresholds: :class:`Sequence` of :class:`float`
            A series of similarity thresholds to apply as spectra are added to clusters
            and as clusters are iteratively refined.
        peak_getter: :class:`Callable`
            A callable object used to get peaks from elements of ``scans``.
        '''
        peak_getter = cls._guess_peak_getter(peak_getter)
        if similarity_thresholds is None:
            similarity_thresholds = [0.1, .4, 0.7]
        singletons = []
        to_bisect = [scans]
        logger = LogUtilsMixin()
        track_similarity = [False] * (len(similarity_thresholds) - 1)
        track_similarity.append(True)
        for similarity_threshold, track in zip(similarity_thresholds, track_similarity):
            logger.log("Clustering with Threshold %0.2f" % (similarity_threshold, ))
            next_to_bisect = []
            if len(to_bisect) > 1:
                logger.log("Refining %d Clusters" % (len(to_bisect)))
            elif to_bisect:
                logger.log("Clustering %d Scans" % (len(to_bisect[0])))
            else:
                logger.log("Nothing to cluster...")
                break
            n = len(to_bisect)
            report_interval = max(min(n // 10, 1000), 1)
            for i, group in enumerate(to_bisect):
                if i % report_interval == 0:
                    logger.log("... Handling Batch %d (%d Scans)" % (i, len(group)))
                clusters = cls.cluster_scans(
                    group, precursor_error_tolerance,
                    minimum_similarity=similarity_threshold,
                    peak_getter=peak_getter, sort=True, track_incremental_similarity=track)
                for cluster in clusters:
                    if len(cluster) == 1:
                        singletons.append(cluster)
                    else:
                        next_to_bisect.append(cluster)
            logger.log("%d Singletons and %d Groups" % (len(singletons), len(next_to_bisect)))
            to_bisect = next_to_bisect
        return SpectrumClusterCollection(sorted(list(singletons) + list(to_bisect)))


cluster_scans = ScanClusterBuilder.cluster_scans
iterative_clustering = ScanClusterBuilder.iterative_clustering


class ScanClusterWriterBase(object):
    '''A base class for writing :class:`ScanCluster` objects to an
    I/O stream like a file.

    Attributes
    ----------
    stream: :class:`io.IOBase`
        The stream to write the clusters to
    metadata: :class:`dict`
        A set of key-value pairs that describe this
        collection.
    '''

    def __init__(self, stream, metadata=None):
        self.stream = stream
        self.metadata = metadata or {}
        self._wrote_metadata = False

    def _write(self, data):
        self.stream.write(data)

    def save(self, cluster):
        '''Write ``cluster`` to the output stream, recording its
        members and calculating it's average similarity.

        Parameters
        ----------
        cluster: :class:`SpectrumCluster`
            The spectrum cluster to write out
        '''
        if not self._wrote_metadata:
            self.write_metadata()
            self._wrote_metadata = True
        self._save(cluster)

    def _save(self, cluster):
        raise NotImplementedError()

    def save_all(self, clusters):
        '''Write each :class:`SpectrumCluster` in ``clusters`` out,
        calling :meth:`save` on each one.

        Parameters
        ----------
        clusters: :class:`collections.Iterable` of :class:`SpectrumCluster`
            The spectrum clusters to write out
        '''
        raise NotImplementedError()

    def add_metadata(self, key, value):
        '''Add metadata to the writer. That metadata will be flushed
        out upon starting to write clusters out.

        Parameters
        ----------
        key: :class:`str`
            The metadata element's name
        value: :class:`str`, :class:`float`
            The metadata element's value
        '''
        if self.wrote_metadata:
            raise TypeError(
                "Cannot add additional metadata, the metadata has already been written")
        self.metadata[key] = value

    def write_metadata(self):
        '''Write the accumulated metadata out in a format-appropriate
        manner at the top of the file.
        '''
        if self._wrote_metadata:
            raise TypeError("Already wrote metadata!")


class ScanClusterWriter(ScanClusterWriterBase):
    '''Writes :class:`ScanCluster` objects to a hierarchical text stream

    Parameters
    ----------
    stream: :class:`io.IOBase`
        The stream to write the clusters to
    '''

    def __init__(self, stream, metadata=None):
        super(ScanClusterWriter, self).__init__(stream, metadata)

    def write_metadata(self):
        for key, value in self.metadata.items():
            self._write("#%s = %s\n" % (key, value))
        self._write("\n")

    def _save(self, cluster):
        '''Write ``cluster`` as a tab delimited tree, recording its
        members and calculating it's average similarity.

        Parameters
        ----------
        cluster: :class:`SpectrumCluster`
            The spectrum cluster to write out
        '''
        self._write("%f\t%d\t%f\n" % (cluster.neutral_mass, len(cluster), cluster.average_similarity()))
        for member in cluster:
            member_source = member.source
            if member_source is not None:
                source_name = member_source.source_file
                if not isinstance(source_name, basestring):
                    if hasattr(source_name, 'name'):
                        source_name = source_name.name
            else:
                source_name = ":detatched:"
            self._write("\t%s\t%s\n" % (source_name, member.id))
        self._write('\n')

    def save_all(self, clusters):
        '''Write each :class:`SpectrumCluster` in ``clusters`` out,
        calling :meth:`save` on each one.

        Parameters
        ----------
        clusters: :class:`collections.Iterable` of :class:`SpectrumCluster`
            The spectrum clusters to write out
        '''
        for cluster in clusters:
            self.save(cluster)


class JSONScanClusterWriter(ScanClusterWriterBase):
    '''A :class:`ScanClusterWriterBase` that uses JSON
    lines.
    '''

    def write_metadata(self):
        json.dump(self.metadata, self.stream)
        self._write("\n")

    def _save(self, cluster):
        json.dump(cluster.to_dict(), self.stream)
        self._write("\n")

    def save_all(self, clusters):
        for cluster in clusters:
            self.save(cluster)


class ScanClusterReaderBase(object):
    '''Base class for reading spectrum clusters from disk.

    Attributes
    ----------
    stream: :class:`io.IOBase`
        The stream to read the clusters from
    resolver_map: :class:`dict`
        A mapping from scan source name to a :class:`Callable`
        which will return a :class:`~.ScanBase` object representing
        that spectrum.
    metadata: :class:`dict`
        A set of key-value pairs that describe this
        collection.
    clusters: :class:`list`
        The read clusters.
    '''
    def __init__(self, stream, resolver_map):
        self.stream = stream
        self.resolver_map = resolver_map
        self.clusters = list()
        self.metadata = {}
        self._generator = None

    def _resolve(self, source, scan_id):
        resolver = self.resolver_map[source]
        try:
            return resolver.get_scan_by_id(scan_id)
        except AttributeError:
            return resolver(scan_id)

    def _parse(self):
        '''Parse the cluster collection from :attr:`stream`
        '''
        self._load_metadata()
        return self._load_clusters()

    def _load_metadata(self):
        '''Read the metadata header from :attr:`stream`.
        '''
        raise NotImplementedError()

    def _load_clusters(self):
        '''Read the data describing :class:`SpectrumCluster` objects from
        :attr:`stream`.
        '''
        raise NotImplementedError()

    def __iter__(self):
        return self

    def __next__(self):
        '''Advance the iterator, retrieving the next :class:`SpectrumCluster`

        Returns
        -------
        :class:`SpectrumCluster`
        '''
        if self._generator is None:
            self._generator = self._parse()
        return next(self._generator)

    def next(self):
        '''Advance the iterator, retrieving the next :class:`SpectrumCluster`

        Returns
        -------
        :class:`SpectrumCluster`
        '''
        return self.__next__()


class ScanClusterReader(ScanClusterReaderBase):
    '''Reads :class:`SpectrumCluster` objects from hierarchical text files written by
    :class:`ScanClusterWriter`.
    '''
    def __init__(self, stream, resolver_map):
        super(ScanClusterReader, self).__init__(stream, resolver_map)
        self._line_buffer = deque()

    def _next_line(self):
        if self._line_buffer:
            return self._line_buffer.popleft()
        return self.stream.readline()

    def _return_line(self, line):
        self._line_buffer.append(line)

    def _stream_lines(self):
        line = self._next_line()
        while line:
            yield line
            line = self._next_line()

    def _load_metadata(self):
        line = self._next_line()
        while line.startswith("#"):
            key, value = line.strip().split(" = ", 1)
            try:
                value = float(value)
            except ValueError:
                value = str(value)
            self.metadata[key] = value
        self._return_line(line)

    def _load_clusters(self):
        current_cluster = []
        mass = None
        similarity = None
        for line in self._stream_lines():
            line = line.rstrip()
            if not line or not line.startswith('\t'):
                if current_cluster:
                    yield SpectrumCluster(current_cluster, mass, similarity)
                    current_cluster = []
                    mass = None
                    similarity = None
                if line:
                    tokens = line.split('\t')
                    mass = float(tokens[0])
                    # size = int(tokens[1])
                    similarity = float(tokens[2])
            elif line.startswith("\t"):
                tokens = line.split("\t")
                source = tokens[1]
                scan_id = tokens[2]
                scan = self._resolve(source, scan_id)
                current_cluster.append(scan)
        if current_cluster:
            yield SpectrumCluster(current_cluster, mass, similarity)


class JSONScanClusterReader(ScanClusterReader):
    '''Reads :class:`SpectrumCluster` objects from hierarchical text files written by
    :class:`JSONScanClusterWriter`.
    '''
    def _load_metadata(self):
        return json.loads(self.stream.readline())

    def _load_clusters(self):
        line = self.stream.readline()
        while line:
            data = json.loads(line)
            mass = data['neutral_mass']
            similarity = data['average_similarity']
            scans = []
            for scan_bundle in data['scans']:
                scans.append(self._resolve(scan_bundle['source'], scan_bundle['id']))
            yield SpectrumCluster(scans, mass, similarity)
            line = self.stream.readline()
