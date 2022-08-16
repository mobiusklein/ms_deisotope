from ms_deisotope.utils import Base
from ms_deisotope.feature_map.feature_map import NeutralMassIndex


def ppm_error(x, y):
    return abs(x - y) / y


def binsearch(array, value):
    lo = 0
    hi = len(array)
    while hi != lo:
        mid = (hi + lo) / 2
        point = array[mid]
        if value == point:
            return mid
        elif hi - lo == 1:
            return mid
        elif point > value:
            hi = mid
        else:
            lo = mid


class FeatureRelationshipBase(object):
    def __repr__(self):
        return "{self.__class__.__name__}()".format(self=self)

    def __hash__(self):
        return hash(repr(self))


class ChargeRelationship(FeatureRelationshipBase):
    def __init__(self, charge_diff, error_tolerance=2e-5):
        self.charge_diff = charge_diff
        self.error_tolerance = error_tolerance

    def __repr__(self):
        return "{self.__class__.__name__}(charge_diff={self.charge_diff})".format(self=self)

    def test(self, a, b):
        if ppm_error(a.neutral_mass, b.neutral_mass) < self.error_tolerance:
            if (a.charge - b.charge) == self.charge_diff:
                return True
        else:
            return False

    def find(self, feature, index):
        hits = index.find_all(feature.neutral_mass, self.error_tolerance)
        for hit in hits:
            if self.test(feature, hit):
                relation = Relation(feature, hit, self)
                feature.supporters.append(relation)


class MassShiftRelationship(FeatureRelationshipBase):
    def __init__(self, mass_shift, error_tolerance=2e-5):
        self.mass_shift = mass_shift
        self.error_tolerance = error_tolerance

    def __repr__(self):
        return "{self.__class__.__name__}(mass_shift={self.mass_shift})".format(self=self)

    def test(self, a, b):
        return ppm_error(
            a.neutral_mass + self.mass_shift, b.neutral_mass) < self.error_tolerance

    def find(self, feature, index):
        hits = index.find_all(feature.neutral_mass + self.mass_shift, self.error_tolerance)
        for hit in hits:
            if self.test(feature, hit):
                relation = Relation(feature, hit, self)
                feature.supporters.append(relation)


class FittedRelationshipFunction(object):
    def __init__(self, relationship, coefficient):
        self.relationship = relationship
        self.coefficient = coefficient

    def test(self, a, b):
        return self.relationship.test(a, b)


class Relation(Base):
    def __init__(self, base, reference, relation):
        self.base = base
        self.reference = reference
        self.relation = relation


class FeatureRelationshipFitter(object):
    def __init__(self, feature_relations=None):
        if feature_relations is None:
            feature_relations = []
        self.feature_relations = feature_relations
        self.history = []

    def add_history(self, hist):
        self.history.append(hist)

    def fit(self, features, past_reference=None, minimum_mass=500):
        relations = []
        features = NeutralMassIndex(features)
        if len(features) == 0:
            return relations
        start_ix = binsearch([f.neutral_mass for f in features], minimum_mass)
        for base in features[start_ix:]:
            for feat_rel_rule in self.feature_relations:
                feat_rel_rule.find(base, features)

        if past_reference is not None and len(past_reference) > 0:
            past_reference = NeutralMassIndex(past_reference)
            for base in features[start_ix:]:
                for feat_rel_rule in self.feature_relations:
                    feat_rel_rule.find(base, past_reference)

    def predict(self, features):
        for feature in features:
            feature.score *= 1.0 + (len(feature.supporters) * 0.05)
