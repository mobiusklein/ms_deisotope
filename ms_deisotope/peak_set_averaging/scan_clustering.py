from .peak_set_similarity import peak_set_similarity


def ppm_error(x, y):
    return (x - y) / y


class SpectrumCluster(object):
    def __init__(self, scans=None):
        if scans is None:
            scans = []
        self.scans = scans

    @property
    def neutral_mass(self):
        return self.scans[0].precursor_information.neutral_mass

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


def cluster_scans(scans, precursor_error_tolerance=1e-5):
    scans = sorted(scans, key=lambda x: sum(p.intensity for p in x.peak_set), reverse=True)
    clusters = []
    for scan in scans:
        best_cluster = None
        best_similarity = 0.0
        for cluster in clusters:
            if abs(ppm_error(scan.precursor_information.neutral_mass,
                             cluster.neutral_mass)) > precursor_error_tolerance:
                continue
            similarity = peak_set_similarity(scan, cluster[0])
            if similarity > best_similarity:
                best_similarity = similarity
                best_cluster = cluster
        if best_cluster is None:
            cluster = SpectrumCluster([scan])
            clusters.append(cluster)
        else:
            best_cluster.append(scan)
    return clusters
