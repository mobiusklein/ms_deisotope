
from functools import total_ordering
from collections import namedtuple, Counter


def isclose(a, b):
    return abs(a - b) < 1e-3


CHARGE_UNKNOWN = 0
OUT_OF_RANGE_INT = 999


def intensity_ratio_function(peak1, peak2):
    ratio = peak1.intensity / float(peak2.intensity)
    if ratio >= 5:
        return -4
    elif 2.5 <= ratio < 5:
        return -3
    elif 1.7 <= ratio < 2.5:
        return -2
    elif 1.3 <= ratio < 1.7:
        return -1
    elif 1.0 <= ratio < 1.3:
        return 0
    elif 0.8 <= ratio < 1.0:
        return 1
    elif 0.6 <= ratio < 0.8:
        return 2
    elif 0.4 <= ratio < 0.6:
        return 3
    elif 0.2 <= ratio < 0.4:
        return 4
    elif 0. <= ratio < 0.2:
        return 5


class DeltaFeature(object):
    def __init__(self, offset, tolerance, intensity_ratio, charge_from, charge_to, name=None):
        self.offset = offset
        self.tolerance = tolerance
        self.intensity_ratio = intensity_ratio
        self.charge_from = charge_from
        self.charge_to = charge_to
        self._hash = hash((self.offset, self.intensity_ratio, self.from_charge,
                           self.to_charge))

    def test(self, peak1, peak2):
        if (self.intensity_ratio == OUT_OF_RANGE_INT or
                intensity_ratio_function(peak1, peak2) == self.intensity_ratio) and\
            ((self.from_charge == OUT_OF_RANGE_INT and self.to_charge == OUT_OF_RANGE_INT) or
             (self.from_charge == peak1.charge and self.to_charge == peak2.charge)):

            return abs((peak1.neutral_mass + self.offset - peak2.neutral_mass) / peak2.neutral_mass) <= self.tolerance
        return False

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        v = abs(self.offset - other.offset) <= 1e-3
        if not v:
            return v
        v = self.intensity_ratio == other.intensity_ratio
        if not v:
            return v
        v = self.from_charge == other.from_charge
        if not v:
            return v
        v = self.to_charge == other.to_charge
        if not v:
            return v
        return True

    def __lt__(self, other):
        eq_count = 0
        v = self.intensity_ratio <= other.intensity_ratio
        eq_count += self.intensity_ratio == other.intensity_ratio
        if not v:
            return v
        v = self.from_charge <= other.from_charge
        eq_count += self.from_charge == other.from_charge
        if not v:
            return v
        v = self.to_charge <= other.to_charge
        eq_count += self.to_charge == other.to_charge
        if not v:
            return v
        if eq_count == 3:
            return False
        return True

    def __gt__(self, other):
        eq_count = 0
        v = self.intensity_ratio >= other.intensity_ratio
        eq_count += self.intensity_ratio == other.intensity_ratio
        if not v:
            return v
        v = self.from_charge >= other.from_charge
        eq_count += self.from_charge == other.from_charge
        if not v:
            return v
        v = self.to_charge >= other.to_charge
        eq_count += self.to_charge == other.to_charge
        if not v:
            return v
        if eq_count == 3:
            return False
        return True

    def __ne__(self, other):
        return not (self == other)

    def __call__(self, peak1, peak2, structure=None):
        return self.test(peak1, peak2)

    def find_matches(self, peak, peak_list):
        matches = []
        for peak2 in peak_list.all_peaks_for(peak.neutral_mass + self.offset, self.tolerance):
            if self(peak, peak2) and peak is not peak2:
                matches.append(peak2)
        return matches

    def is_valid_match(self, peak, match, labeled_peaks):
        return match in labeled_peaks

    def specialize(self, from_charge, to_charge, intensity_ratio):
        return self.__class__(
            self.offset, self.tolerance, self.name, intensity_ratio,
            from_charge, to_charge)

    def unspecialize(self):
        return self.__class__(
            self.offset, self.tolerance, self.name, OUT_OF_RANGE_INT,
            OUT_OF_RANGE_INT, OUT_OF_RANGE_INT)

    def _get_display_fields(self):
        fields = {}
        fields['offset'] = self.offset
        if self.from_charge != OUT_OF_RANGE_INT:
            fields["from_charge"] = self.from_charge
        if self.to_charge != OUT_OF_RANGE_INT:
            fields['to_charge'] = self.to_charge
        if self.intensity_ratio != OUT_OF_RANGE_INT:
            fields["intensity_ratio"] = self.intensity_ratio
        terms = []
        for k, v in fields.items():
            if isinstance(v, int):
                terms.append("%s=%d" % (k, v))
            elif isinstance(v, float):
                terms.append("%s=%0.4f" % (k, v))
            else:
                terms.append("%s=%r" % (k, v))
        return terms

    def __repr__(self):
        terms = self._get_display_fields()
        return "{}(name={!r}, {})".format(
            self.__class__.__name__, self.name, ", ".join(terms))


def feature_function_estimator(labeled_scans, feature_function, tolerance=2e-5, track_relations=True, verbose=False):
    total_on_series_satisfied = 0.
    total_off_series_satisfied = 0.
    total_on_series = 0.
    total_off_series = 0.
    peak_relations = []
    for i_scan, labeled_scan in enumerate(labeled_scans):
        if verbose and i_scan % 1000 == 0:
            print(i_scan)
        peaks = labeled_scan.deconvoluted_peak_set
        related = []
        labeled_peaks = labeled_scan.annotations['labeled_peaks']
        for peak in peaks:
            is_on_series = peak in labeled_peaks
            matches = feature_function.find_matches(peak, peaks)
            for match in matches:
                if peak is match:
                    continue
                if track_relations:
                    pr = PeakRelation(
                        peak, match, feature_function, intensity_ratio_function(peak, match))
                    related.append(pr)
                is_match_expected = feature_function.is_valid_match(peak, match, labeled_peaks)
                if is_on_series and is_match_expected:
                    total_on_series_satisfied += 1
                else:
                    total_off_series_satisfied += 1
            if is_on_series:
                total_on_series += 1
            else:
                total_off_series += 1
        if len(related) > 0 and track_relations:
            peak_relations.append((labeled_scan, related))

    total_on_series_satisfied_normalized = total_on_series_satisfied / \
        max(total_on_series, 1)
    total_off_series_satisfied_normalized = total_off_series_satisfied / \
        max(total_off_series, 1)

    return FittedFeature(feature_function, total_on_series_satisfied_normalized,
                         total_off_series_satisfied_normalized, peak_relations,
                         on_count=total_on_series_satisfied,
                         off_count=total_off_series_satisfied)


class FittedFeatureBase(object):
    def find_matches(self, peak, peak_list):
        return self.feature.find_matches(peak, peak_list)

    def is_valid_match(self, from_peak, to_peak, labeled_peaks):
        return self.feature.is_valid_match(from_peak, to_peak, labeled_peaks)

    def _feature_probability(self, p=0.5):
        return (p * self.on_series) / (
            (p * self.on_series) + ((1 - p) * self.off_series))


@total_ordering
class FittedFeature(FittedFeatureBase):
    def __init__(self, feature, on_series, off_series, relations=None, on_count=0, off_count=0):
        if relations is None:
            relations = []
        self.feature = feature

        # forward these attributes directly rather than using a property
        # to avoid overhead
        self.from_charge = feature.from_charge
        self.to_charge = feature.to_charge

        # mu
        self.on_series = on_series
        # v
        self.off_series = off_series

        # tracking purposes only
        self.on_count = on_count
        self.off_count = off_count
        self.relations = relations

    @property
    def _total_on_series(self):
        if self.on_series == 0:
            return 1
        return self.on_count / self.on_series

    @property
    def _total_off_series(self):
        if self.off_series == 0:
            return 1
        return self.off_count / self.off_series

    @property
    def name(self):
        return self.feature.name

    @property
    def observations(self):
        count = (self.on_count + self.off_count)
        if count > 0:
            return count
        count = len(list(self.peak_relations()))
        if count > 0:
            return count
        return 0

    def __hash__(self):
        return hash((self.feature, self.expected))

    def __eq__(self, other):
        v = self.feature == other.feature and self.expected == other.expected
        if not v:
            return v
        v = isclose(self.on_series, other.on_series) and isclose(
            self.off_series, other.off_series)
        if not v:
            return v
        return True

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        v = self.feature < other.feature
        if not v:
            return False
        v = self.on_series < other.on_series
        if not v:
            return False
        v = self.off_series < other.off_series
        if not v:
            return False
        return True

    def __repr__(self):
        temp = ("<FittedFeature {feature.name}, {terms} u:{on_series:0.4g}"
                " v:{off_series:0.4g} @ {series} {count_relations}>")
        return temp.format(
            feature=self.feature,
            terms=', '.join(map(str, self.feature._get_display_fields())),
            on_series=self.on_series, off_series=self.off_series,
            series=self.expected, count_relations=self.on_count + self.off_count)

    def charge_relations(self):
        counter = Counter()
        for rel in self.peak_relations(False):
            counter[rel.from_charge, rel.to_charge] += 1
        return counter

    def intensity_ratio(self):
        counter = Counter()
        for rel in self.peak_relations(False):
            counter[intensity_ratio_function(rel.from_peak, rel.to_peak)] += 1
        return counter

    def charge_intensity_ratio(self):
        counter = Counter()
        for rel in self.peak_relations(False):
            counter[(rel.from_charge, rel.to_charge),
                    intensity_ratio_function(rel.from_peak, rel.to_peak)] += 1
        return counter

    def peak_relations(self, include_noise=True):
        for _, peak_relations in self.relations:
            for pr in peak_relations:
                if not include_noise and not pr.expected:
                    continue
                yield pr

    def partitions(self, minimum_count=10):
        counter = Counter()
        on_counter = Counter()
        off_counter = Counter()
        for rel in self.peak_relations(True):
            fs = FeatureSpecialization(
                from_charge=rel.from_charge,
                to_charge=rel.to_charge,
                intensity_ratio=rel.intensity_ratio,
            )
            counter[fs] += 1
            if rel.expected == self.expected:
                on_counter[fs] += 1
            else:
                off_counter[fs] += 1

        counter = {k: (v, on_counter[k], off_counter[k])
                   for k, v in counter.items() if v >= minimum_count}
        return counter

    def specialize(self, minimum_count=10):
        counter = self.partitions(minimum_count)
        total_off_series = self._total_off_series
        total_on_series = self._total_on_series
        specialized_features = []
        for params, counts in counter.items():
            feature = self.feature.specialize(
                intensity_ratio=params.intensity_ratio,
                from_charge=params.from_charge,
                to_charge=params.to_charge)
            total, on, off = counts
            fit = FittedFeature(
                feature, on / total_on_series, off / total_off_series,
                None, on, off)
            specialized_features.append(fit)
        return specialized_features

    def __call__(self, peak1, peak2, structure=None):
        return self.feature(peak1, peak2, structure)

    def pack(self):
        self.relations = []
        return self

    def to_json(self):
        d = {}
        d['on_series'] = self.on_series
        d['on_count'] = self.on_count
        d['off_series'] = self.off_series
        d['off_count'] = self.off_count
        d['feature'] = self.feature.to_json()
        return d

    @classmethod
    def from_json(cls, d):
        feature = DeltaFeature.from_json(d['feature'])
        inst = cls(
            feature, d["on_series"],
            d["off_series"], relations=None, on_count=d['on_count'],
            off_count=d['off_count']
        )
        return inst


FeatureSpecialization = namedtuple("FeatureSpecialization", [
                                   "from_charge", "to_charge", "intensity_ratio"])


class PeakRelation(object):
    __slots__ = ["from_peak", "to_peak", "feature",
                 "intensity_ratio", "expected", "from_charge",
                 "to_charge"]

    def __init__(self, from_peak, to_peak, feature, intensity_ratio=None, expected=None):
        if intensity_ratio is None:
            intensity_ratio = intensity_ratio_function(from_peak, to_peak)
        self.from_peak = from_peak
        self.to_peak = to_peak
        self.feature = feature
        self.intensity_ratio = intensity_ratio
        self.from_charge = from_peak.charge
        self.to_charge = to_peak.charge
        self.expected = expected

    def __reduce__(self):
        return self.__class__, (self.from_peak, self.to_peak, self.feature, self.intensity_ratio, self.expected)

    def __repr__(self):
        template = "<PeakRelation {s.from_peak.neutral_mass}({s.from_charge}) ->" +\
            " {s.to_peak.neutral_mass}({s.to_charge}) by {s.feature.name} on {s.expected}>"
        return template.format(s=self)

    def peak_key(self):
        if self.from_peak.index.neutral_mass < self.to_peak.index.neutral_mass:
            return self.from_peak, self.to_peak
        else:
            return self.to_peak, self.from_peak
