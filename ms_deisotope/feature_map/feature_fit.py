from collections import namedtuple

from brainpy import neutral_mass as calc_neutral_mass
from .lcms_feature import (
    LCMSFeature,
    LCMSFeatureTreeNode,
    FeatureSetIterator)


class map_coord(namedtuple("map_coord", ("mz", 'time'))):
    def __repr__(self):
        return "(%0.3f, %0.3f)" % self


class LCMSFeatureSetFit(object):
    def __init__(self, features, theoretical, score, charge,
                 missing_features=0, supporters=None, data=None,
                 neutral_mass=None):
        if supporters is None:
            supporters = []
        self.features = features
        self.theoretical = theoretical
        self.score = score
        self.charge = charge
        self.data = data
        self.missing_features = missing_features
        self.monoisotopic_feature = features[0]
        self.supporters = supporters
        if neutral_mass is None:
            neutral_mass = neutral_mass(self.monoisotopic_feature.mz, self.charge)
        self.neutral_mass = neutral_mass
        self.mz = self.monoisotopic_feature.mz

    def clone(self):
        return self.__class__(
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.supporters, self.data, self.neutral_mass)

    def __reduce__(self):
        return self.__class__, (
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.supporters, self.data, self.neutral_mass)

    def __eq__(self, other):
        val = (self.score == other.score and
               self.charge == other.charge and
               self.features == other.features and
               self.theoretical == other.theoretical)
        if self.data is not None or other.data is not None:
            val = val and (self.data == other.data)
        return val

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.score < other.score

    def __gt__(self, other):
        return self.score > other.score

    def __hash__(self):
        return hash((self.monoisotopic_feature.mz, self.charge))

    def __iter__(self):
        return iter(self.features)

    def __len__(self):
        return len(self.features)

    @property
    def npeaks(self):
        return len(self)

    def __repr__(self):
        return "LCMSFeatureSetFit(score=%0.5f, charge=%d, size=%d, monoisotopic_mz=%0.5f, %0.2f-%0.2f)" % (
            self.score, self.charge, len(self), self.monoisotopic_feature.mz,
            self.start.time, self.end.time)

    @property
    def start(self):
        first = self.features[0]
        if first is None:
            raise Exception()
        return map_coord(first.mz, first.start_time)

    @property
    def end(self):
        for last in reversed(self.features):
            if last is None:
                continue
            return map_coord(last.mz, last.end_time)


class DeconvolutedLCMSFeatureTreeNode(LCMSFeatureTreeNode):
    __slots__ = ["_neutral_mass", "charge", "precursor_information"]

    def __init__(self, retention_time=None, members=None, precursor_information=None):
        if precursor_information is None:
            precursor_information = []
        self._neutral_mass = 0
        self.charge = 0
        super(DeconvolutedLCMSFeatureTreeNode, self).__init__(retention_time, members)
        self.precursor_information = precursor_information

    def _recalculate(self):
        self._calculate_most_abundant_member()
        self._mz = self._most_abundant_member.mz
        self._neutral_mass = self._most_abundant_member.neutral_mass
        self.charge = self._most_abundant_member.charge

    @property
    def neutral_mass(self):
        if self._neutral_mass is None:
            if self._most_abundant_member is not None:
                self._neutral_mass = self._most_abundant_member.neutral_mass
        return self._neutral_mass


class DeconvolutedLCMSFeature(LCMSFeature):
    def __init__(self, nodes=None, charge=None, adducts=None, used_as_adduct=None, score=0.0,
                 n_features=0, feature_id=None, supporters=None):
        if supporters is None:
            supporters = []
        self.charge = charge
        self.score = score
        self._neutral_mass = None
        self._last_neutral_mass = None
        self._precursor_information = None
        self.n_features = n_features
        self.supporters = supporters
        super(DeconvolutedLCMSFeature, self).__init__(nodes, adducts, used_as_adduct, feature_id=feature_id)

    @property
    def precursor_information(self):
        if self._precursor_information is None:
            pinfo = []
            for node in self:
                pinfo.extend(node.precursor_information)
            self._precursor_information = tuple(pinfo)
        return self._precursor_information

    def clone(self):
        return DeconvolutedLCMSFeature(
            self.nodes, self.charge, self.adducts, self.used_as_adduct, self.score,
            self.n_features, self.feature_id, list(self.supporters))

    def _invalidate(self):
        self._last_neutral_mass = self._neutral_mass if self._neutral_mass is not None else 0.
        self._neutral_mass = None
        self._precursor_information = None
        super(DeconvolutedLCMSFeature, self)._invalidate()

    @property
    def neutral_mass(self):
        if self._neutral_mass is None:
            best_neutral_mass = 0
            maximum_intensity = 0
            for node in self.nodes:
                intensity = node.max_intensity()
                if intensity > maximum_intensity:
                    maximum_intensity = intensity
                    best_neutral_mass = node.neutral_mass
            self._last_neutral_mass = self._neutral_mass = best_neutral_mass
        return self._neutral_mass

    def __repr__(self):
        return "%s(%0.4f, %d, %0.2f, %0.2f, %0.2f)" % (
            self.__class__.__name__, self.neutral_mass,
            self.charge, self.score,
            self.start_time, self.end_time)


try:
    has_c = True
    _map_coord = map_coord
    _LCMSFeatureSetFit = LCMSFeatureSetFit
    from ms_deisotope._c.feature_map.feature_fit import (LCMSFeatureSetFit, map_coord)
except ImportError:
    has_c = False
