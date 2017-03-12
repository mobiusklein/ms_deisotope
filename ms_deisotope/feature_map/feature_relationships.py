from ms_deisotope.utils import Base


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


class MassShiftRelationship(FeatureRelationshipBase):
    def __init__(self, mass_shift, error_tolerance=2e-5):
        self.mass_shift = mass_shift
        self.error_tolerance = error_tolerance

    def __repr__(self):
        return "{self.__class__.__name__}(mass_shift={self.mass_shift})".format(self=self)

    def test(self, a, b):
        return ppm_error(
            a.neutral_mass + self.mass_shift, b.neutral_mass) < self.error_tolerance


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
        features = sorted(features, key=lambda x: x.neutral_mass)
        if len(features) == 0:
            return relations
        start_ix = binsearch([f.neutral_mass for f in features], 500.)
        for base in features[start_ix:]:
            for reference in features[start_ix:]:
                if base is reference:
                    continue
                for feat_rel_rule in self.feature_relations:
                    if feat_rel_rule.test(base, reference):
                        relation = Relation(base, reference, feat_rel_rule)
                        base.supporters.append(relation)
                        relations.append(relation)

        if past_reference is not None:
            for base in features[start_ix:]:
                for reference in past_reference:
                    for feat_rel_rule in self.feature_relations:
                        if feat_rel_rule.test(base, reference):
                            relation = Relation(base, reference, feat_rel_rule)
                            base.supporters.append(relation)
                            relations.append(relation)
        return relations
