import operator
import functools
import bisect
import json

from .similarity_methods import peak_set_similarity

from ms_deisotope.task import LogUtilsMixin

peak_set_getter = operator.attrgetter("peak_set")
deconvoluted_peak_set_getter = operator.attrgetter("deconvoluted_peak_set")


def ppm_error(x, y):
    return (x - y) / y


@functools.total_ordering
class SpectrumCluster(object):
    def __init__(self, scans=None):
        if scans is None:
            scans = []
        self.scans = scans

    @property
    def neutral_mass(self):
        return self.scans[0].precursor_information.neutral_mass

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

    def append(self, item):
        self.scans.append(item)

    def average_similarity(self, *args, **kwargs):
        n = len(self)
        if n == 1:
            return 1.0
        ratings = []
        for i in range(n):
            scan_i = self[i]
            for j in range(i + 1, n):
                scan_j = self[j]
                ratings.append(peak_set_similarity(scan_i, scan_j, *args, **kwargs))
        return sum(ratings) / len(ratings)

    def to_dict(self):
        d = {}
        d['neutral_mass'] = self.neutral_mass
        d['size'] = len(self)
        d['average_similarity'] = self.average_similarity()
        scans = []
        for scan in self:
            scans.append({
                "id": scan.id,
                "source": scan.source.source_file if scan.source is not None else ":detatched:",
                "neutral_mass": scan.precursor_information.neutral_mass,
            })
        d['scans'] = scans
        return d


class SpectrumClusterCollection(object):
    def __init__(self, clusters=None):
        if clusters is None:
            clusters = []
        self.clusters = list(clusters)

    def add(self, cluster):
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
                    x = self.structures[i]
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
                    x = self.structures[i]
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
        target_ix, lo_ix, hi_ix = self._binary_search(mass, error_tolerance)
        target = self[target_ix]
        if abs(target.neutral_mass - mass) / mass > error_tolerance:
            return None
        return target

    def find_all(self, mass, error_tolerance=1e-5):
        target_ix, lo_ix, hi_ix = self._binary_search(mass, error_tolerance)
        result = [
            target for target in self[lo_ix:hi_ix]
            if abs(target.neutral_mass - mass) / mass <= error_tolerance
        ]
        return result


def binsearch_simple(array, x):
    n = len(array)
    lo = 0
    hi = n

    while hi != lo:
        mid = (hi + lo) // 2
        y = array[mid].neutral_mass
        err = y - x
        # Do refinement in :func:`cluster_scans`
        if hi - lo == 1:
            return mid
        elif err > 0:
            hi = mid
        else:
            lo = mid
    return 0


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
    """

    @classmethod
    def _guess_peak_getter(cls, getter):
        if getter is None:
            return list
        if callable(getter):
            return getter
        if isinstance(getter, basestring):
            if getter == "d":
                return deconvoluted_peak_set_getter
            if getter in ('p', 'c'):
                return peak_set_getter
            else:
                return operator.attrgetter(getter)
        raise ValueError("Cannot infer peak set getter from %r" % (getter, ))

    def __init__(self, clusters=None, precursor_error_tolerance=1e-5, minimum_similarity=0.1,
                 peak_getter=None):
        peak_getter = self._guess_peak_getter(peak_getter)
        if clusters is None:
            clusters = []
        self.clusters = SpectrumClusterCollection(clusters)
        self.precursor_error_tolerance = precursor_error_tolerance
        self.minimum_similarity = minimum_similarity
        self.peak_getter = peak_getter

    def peak_set_similarity(self, scan_i, scan_j):
        return peak_set_similarity(
            self.peak_getter(scan_i),
            self.peak_getter(scan_j))

    def find_best_cluster_for_scan(self, scan):
        best_cluster = None
        best_similarity = 0.0

        if len(self) == 0:
            return best_cluster

        center_i = binsearch_simple(self.clusters, scan.precursor_information.neutral_mass)
        i = center_i

        while i >= 0:
            cluster = self.clusters[i]
            if abs(ppm_error(scan.precursor_information.neutral_mass,
                             cluster.neutral_mass)) > self.precursor_error_tolerance:
                break
            similarity = self.peak_set_similarity(scan, cluster[0])
            i -= 1
            if similarity > best_similarity and similarity > self.minimum_similarity:
                best_similarity = similarity
                best_cluster = cluster

        n = len(self.clusters)
        i = center_i + 1
        while i < n:
            cluster = self.clusters[i]
            if abs(ppm_error(scan.precursor_information.neutral_mass,
                             cluster.neutral_mass)) > self.precursor_error_tolerance:
                break
            similarity = self.peak_set_similarity(scan, cluster[0])
            i += 1
            if similarity > best_similarity and similarity > self.minimum_similarity:
                best_similarity = similarity
                best_cluster = cluster
        return best_cluster

    def add_scan(self, scan):
        best_cluster = self.find_best_cluster_for_scan(scan)
        if best_cluster:
            best_cluster.append(scan)
        else:
            self.clusters.add(SpectrumCluster([scan]))

    def __iter__(self):
        return iter(self.clusters)

    def __len__(self):
        return len(self.clusters)

    def __getitem__(self, i):
        return self.clusters[i]

    @classmethod
    def cluster_scans(cls, scans, precursor_error_tolerance=1e-5, minimum_similarity=0.1,
                      peak_getter=None):
        self = cls([], precursor_error_tolerance, minimum_similarity, peak_getter)
        scans = sorted(scans, key=lambda x: sum(p.intensity for p in self.peak_getter(x)), reverse=True)
        for scan in scans:
            self.add_scan(scan)
        return self.clusters

    @classmethod
    def iterative_clustering(cls, scans, precursor_error_tolerance=1e-5, similarity_thresholds=None,
                             peak_getter=None):
        peak_getter = cls._guess_peak_getter(peak_getter)
        if similarity_thresholds is None:
            similarity_thresholds = [0.1, .4, 0.7]
        singletons = []
        to_bisect = [scans]
        for similarity_threshold in similarity_thresholds:
            next_to_bisect = []
            for group in to_bisect:
                clusters = cls.cluster_scans(
                    group, precursor_error_tolerance,
                    minimum_similarity=similarity_threshold,
                    peak_getter=peak_getter)
                for cluster in clusters:
                    if len(cluster) == 1:
                        singletons.append(cluster)
                    else:
                        next_to_bisect.append(cluster)
            to_bisect = next_to_bisect
        return SpectrumClusterCollection(sorted(list(singletons) + list(to_bisect)))


cluster_scans = ScanClusterBuilder.cluster_scans
iterative_clustering = ScanClusterBuilder.iterative_clustering


class ScanClusterWriter(object):
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)

    def save(self, cluster):
        self.write("%f\t%d\t%f\n" % (cluster.neutral_mass, len(cluster), cluster.average_similarity()))
        for member in cluster:
            member_source = member.source
            if member_source is not None:
                source_name = member_source.source_file
            else:
                source_name = ":detatched:"
            self.write("\t%s\t%s\n" % (source_name, member.id))
        self.write('\n')

    def save_all(self, clusters):
        for cluster in clusters:
            self.save(cluster)


class JSONScanClusterWriter(object):
    def __init__(self, stream):
        self.stream = stream

    def save(self, cluster):
        json.dump(self.stream, cluster.to_dict())

    def save_all(self, clusters):
        for cluster in clusters:
            self.save(cluster)
