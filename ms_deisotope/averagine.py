from collections import defaultdict
from brainpy import calculate_mass, neutral_mass, PROTON, isotopic_variants, mass_charge_ratio

from .utils import dict_proxy


def shift_isotopic_pattern(mz, cluster):
    first_peak = cluster[0]
    for peak in cluster[1:]:
        delta = (peak.mz - first_peak.mz)
        peak.mz = mz + delta

    cluster[0].mz = mz
    return cluster


class TheoreticalIsotopicPattern(object):
    def __init__(self, base_tid, truncated_tid=None):
        if truncated_tid is None:
            truncated_tid = base_tid
        self.base_tid = base_tid
        self.truncated_tid = truncated_tid

    def __getitem__(self, i):
        return self.truncated_tid[i]

    def __iter__(self):
        return iter(self.truncated_tid)

    def __len__(self):
        return len(self.truncated_tid)

    def get(self, i):
        return self.truncated_tid[i]

    def clone(self):
        return self.__class__([p.clone() for p in self.base_tid],
                              [p.clone() for p in self.truncated_tid])

    @property
    def monoisotopic_mz(self):
        return self.base_tid[0].mz

    def shift(self, mz):
        first_peak = self.base_tid[0]
        for peak in self.base_tid[1:]:
            delta = peak.mz - first_peak.mz
            peak.mz = mz + delta
        first_peak.mz = mz
        return self

    def truncate_after(self, truncate_after=0.0):
        cumsum = 0
        result = []
        for peak in self.base_tid:
            cumsum += peak.intensity
            result.append(peak)
            if cumsum >= truncate_after:
                break
        for peak in result:
            peak.intensity *= 1. / cumsum
        self.truncated_tid = result
        return self

    def ignore_below(self, ignore_below=0.0):
        total = 0
        kept_tid = []
        for i, p in enumerate(self.truncated_tid):
            if p.intensity < ignore_below and i > 1:
                continue
            else:
                total += p.intensity
                kept_tid.append(p.clone())
        for p in kept_tid:
            p.intensity /= total
        self.truncated_tid = kept_tid
        return self

    def __repr__(self):
        return "TheoreticalIsotopicPattern(%0.4f, charge=%d, (%s))" % (
            self.base_tid[0].mz,
            self.base_tid[0].charge,
            ', '.join("%0.3f" % p.intensity for p in self.base_tid))


@dict_proxy("base_composition")
class Averagine(object):
    def __init__(self, base_composition):
        self.base_composition = dict(base_composition)
        self.base_mass = calculate_mass(self.base_composition)

    def scale(self, mz, charge=1, charge_carrier=PROTON):
        neutral = neutral_mass(mz, charge, charge_carrier)

        scale = neutral / self.base_mass
        scaled = {}
        for elem, count in self.base_composition.items():
            scaled[elem] = round(count * scale)

        scaled_mass = calculate_mass(scaled)
        delta_hydrogen = round(scaled_mass - neutral)
        H = scaled["H"]
        if H > delta_hydrogen:
            scaled["H"] = H - delta_hydrogen
        else:
            scaled["H"] = 0

        return scaled

    def isotopic_cluster(self, mz, charge=1, charge_carrier=PROTON, truncate_after=0.95, ignore_below=0.0):
        composition = self.scale(mz, charge, charge_carrier)
        cumsum = 0
        result = []
        for peak in isotopic_variants(composition, charge=charge):
            cumsum += peak.intensity
            result.append(peak)
            if cumsum >= truncate_after:
                break
        for peak in result:
            peak.intensity *= 1. / cumsum
        result = shift_isotopic_pattern(mz, result)
        result = TheoreticalIsotopicPattern(result)
        if ignore_below > 0:
            result.ignore_below(ignore_below)
        return result.truncated_tid

    def __repr__(self):
        return "Averagine(%r)" % self.base_composition

    def __eq__(self, other):
        return self.base_composition == other.base_composition

    def __hash__(self):
        return hash(self.base_composition.items())


def average_compositions(compositions):
    n = 0
    result = defaultdict(float)
    for comp in compositions:
        n += 1
        for k, v in comp.items():
            result[k] += v
    for k, v in list(result.items()):
        result[k] = v / n
    return dict(result)


def add_compositions(a, b):
    a = defaultdict(float, **a)
    for k, v in b.items():
        a[k] += v
    return dict(a)


try:
    _Averagine = Averagine
    _TheoreticalIsotopicPattern = TheoreticalIsotopicPattern
    try:
        from ms_deisotope._c.averagine import TheoreticalIsotopicPattern
    except Exception:
        pass
    from ms_deisotope._c.averagine import Averagine
except ImportError as e:
    print(e, "averagine")


peptide = Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})
glycopeptide = Averagine({"C": 10.93, "H": 15.75, "N": 1.6577, "O": 6.4773, "S": 0.02054})
glycan = Averagine({'C': 7.0, 'H': 11.8333, 'N': 0.5, 'O': 5.16666})
permethylated_glycan = Averagine({'C': 12.0, 'H': 21.8333, 'N': 0.5, 'O': 5.16666})
heparin = Averagine({'H': 10.5, 'C': 6, 'S': 0.5, 'O': 5.5, 'N': 0.5})


_neutron_shift = calculate_mass({"C[13]": 1}) - calculate_mass({"C[12]": 1})


def isotopic_shift(charge=1):
    return _neutron_shift / float(charge)


@dict_proxy("averagine")
class AveragineCache(object):
    def __init__(self, averagine, backend=None, cache_truncation=1.0):
        if backend is None:
            backend = {}
        self.backend = backend
        self.averagine = Averagine(averagine)
        self.cache_truncation = cache_truncation

    def has_mz_charge_pair(self, mz, charge=1, charge_carrier=PROTON, truncate_after=0.95, ignore_below=0.0):
        if self.cache_truncation == 0.0:
            key_mz = mz
        else:
            key_mz = round(mz / self.cache_truncation) * self.cache_truncation
        if (key_mz, charge, charge_carrier) in self.backend:
            return shift_isotopic_pattern(
                mz, [p.clone() for p in self.backend[key_mz, charge, charge_carrier]])
        else:
            tid = self.averagine.isotopic_cluster(
                mz, charge, charge_carrier, truncate_after, ignore_below)
            self.backend[key_mz, charge, charge_carrier] = [p.clone() for p in tid]
            return tid

    isotopic_cluster = has_mz_charge_pair

    def __repr__(self):
        return "AveragineCache(%r)" % self.averagine

    def clear(self):
        self.backend.clear()


try:
    _AveragineCache = AveragineCache
    _isotopic_shift = isotopic_shift
    from ms_deisotope._c.averagine import AveragineCache, isotopic_shift
except ImportError:
    pass
