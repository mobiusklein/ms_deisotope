from collections import namedtuple

import numpy as np

# from brainpy import neutral_mass as calc_neutral_mass
from ms_peak_picker import FittedPeak
from ms_deisotope.averagine import glycan
from ms_deisotope.scoring import g_test_scaled
from .shape_fitter import AdaptiveMultimodalChromatogramShapeFitter
from .lcms_feature import (
    EmptyFeature,
    LCMSFeature,
    LCMSFeatureTreeNode,
    RunningWeightedAverage,
    NodeFeatureSetIterator)


class map_coord(namedtuple("map_coord", ("mz", 'time'))):
    def __repr__(self):
        return "(%0.3f, %0.3f)" % self


class LCMSFeatureSetFit(object):
    def __init__(self, features, theoretical, score, charge,
                 missing_features=0, supporters=None, data=None,
                 neutral_mass=None, n_points=0, scores=None, times=None):
        if supporters is None:
            supporters = []
        if scores is None:
            scores = np.array([])
        if times is None:
            times = np.array([])
        self.features = features
        self.theoretical = theoretical
        self.score = score
        self.charge = charge
        self.data = data
        self.n_points = n_points
        self.missing_features = missing_features
        self.monoisotopic_feature = features[0]
        self.supporters = supporters
        self.mz = theoretical.monoisotopic_mz
        if neutral_mass is None:
            neutral_mass = neutral_mass(self.mz, self.charge)
        self.neutral_mass = neutral_mass
        self.scores = scores
        self.times = times

    def count_null_features(self):
        n_null = 0
        for feature in self.features:
            if feature is None or isinstance(feature, EmptyFeature):
                n_null += 1
        return n_null

    def has_multiple_real_features(self):
        return len(self) - self.count_null_features() > 1

    def clone(self):
        return self.__class__(
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.supporters, self.data,
            self.neutral_mass, self.scores, self.times)

    def __reduce__(self):
        return self.__class__, (
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.supporters, self.data, self.neutral_mass,
            self.scores, self.times)

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

    def __init__(self, time=None, members=None, precursor_information=None):
        if precursor_information is None:
            precursor_information = []
        self._neutral_mass = 0
        self.charge = 0
        super(DeconvolutedLCMSFeatureTreeNode, self).__init__(time, members)
        self.precursor_information = precursor_information

    def _recalculate(self):
        self._calculate_most_abundant_member()
        if self._most_abundant_member is not None:
            self._mz = self._most_abundant_member.mz
            self._neutral_mass = self._most_abundant_member.neutral_mass
            self.charge = self._most_abundant_member.charge

    @property
    def neutral_mass(self):
        if self._neutral_mass == 0:
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

    def __reduce__(self):
        return self.__class__, (self.nodes, self.charge, self.adducts,
                                self.used_as_adduct, self.score, self.feature_id, self.supporters)

    @property
    def precursor_information(self):
        if self._precursor_information is None:
            pinfo = []
            for node in self:
                pinfo.extend(node.precursor_information)
            self._precursor_information = tuple(pinfo)
        return self._precursor_information

    def clone(self, deep=False, cls=None):
        if cls is None:
            cls = self.__class__
        return cls(
            self.nodes.clone(deep=deep), self.charge, self.adducts, self.used_as_adduct, self.score,
            self.n_features, self.feature_id, list(self.supporters))

    def _invalidate(self, reaverage=False):
        self._last_neutral_mass = self._neutral_mass if self._neutral_mass is not None else 0.
        self._neutral_mass = None
        self._precursor_information = None
        super(DeconvolutedLCMSFeature, self)._invalidate(reaverage)

    @property
    def neutral_mass(self):
        if self._neutral_mass is None:
            avger = DeconvolutedRunningWeightedAverage()
            for node in self.nodes:
                avger.update(node.members)
            self._neutral_mass = self._last_neutral_mass = avger.current_mean
        return self._neutral_mass

    def _copy_chunk(self, nodes, *args, **kwargs):
        x = self.__class__(
            nodes, self.charge, list(self.adducts), list(self.used_as_adduct),
            self.score, self.n_features, None, list(self.supporters))
        return x

    def sum(self, other):
        missing = []
        feat_iter = NodeFeatureSetIterator([self, other])
        for nodes in feat_iter:
            base = nodes[0]
            new = nodes[1]
            if base is None:
                missing.append(new)
            elif new is not None:
                base.members[0].intensity += new.members[0].intensity
                base.precursor_information.extend(new.precursor_information)
        if missing:
            for node in missing:
                self.insert_node(DeconvolutedLCMSFeatureTreeNode(
                    node.time, list(node.members), list(node.precursor_information)))
        self.supporters.extend(other.supporters)
        return self

    def __repr__(self):
        return "%s(%0.4f, %d, %0.2f, %0.2f, %0.2f)" % (
            self.__class__.__name__, self.neutral_mass,
            self.charge, self.score,
            self.start_time, self.end_time)


class DeconvolutedRunningWeightedAverage(RunningWeightedAverage):
    def add(self, peak):
        if peak.intensity == 0:
            if self.current_mean == 0 and self.total_weight == 0:
                self.current_mean = peak.neutral_mass
                self.total_weight = 1
            else:
                return
        agg = (self.total_weight * self.current_mean) + \
            (peak.neutral_mass * peak.intensity)
        self.total_weight += peak.intensity
        self.current_mean = agg / self.total_weight
        self.current_count += 1
        return self


class DriftTimeRunningWeightedAverage(RunningWeightedAverage):
    def add(self, peak):
        if peak.intensity == 0:
            if self.current_mean == 0 and self.total_weight == 0:
                self.current_mean = peak.drift_time
                self.total_weight = 1
            else:
                return
        agg = (self.total_weight * self.current_mean) + \
            (peak.drift_time * peak.intensity)
        self.total_weight += peak.intensity
        self.current_mean = agg / self.total_weight
        self.current_count += 1
        return self


class IonMobilityDeconvolutedLCMSFeature(DeconvolutedLCMSFeature):
    def __init__(self, nodes=None, charge=None, adducts=None, used_as_adduct=None, score=0.0,
                 n_features=0, feature_id=None, supporters=None):
        self._drift_time = None
        self._last_drift_time = None
        super(IonMobilityDeconvolutedLCMSFeature, self).__init__(
            nodes=nodes, charge=charge, adducts=adducts, used_as_adduct=used_as_adduct, score=score,
            n_features=n_features, feature_id=feature_id, supporters=supporters)

    def _invalidate(self, reaverage=False):
        self._last_drift_time = self._drift_time if self._drift_time is not None else 0.
        self._drift_time = None
        return super(IonMobilityDeconvolutedLCMSFeature, self)._invalidate(reaverage=reaverage)

    @property
    def drift_time(self):
        if self._drift_time is None:
            avger = DriftTimeRunningWeightedAverage()
            for node in self.nodes:
                avger.update(node.members)
            self._drift_time = self._last_drift_time = avger.current_mean
        return self._drift_time

    def __repr__(self):
        return "%s(%0.4f, %0.4f, %d, %0.2f, %0.2f, %0.2f)" % (
            self.__class__.__name__, self.neutral_mass, self.drift_time,
            self.charge, self.score,
            self.start_time, self.end_time)


def envelope_to_peak_list(envelope):
    return [FittedPeak(e[0], e[1], 0, 0, 0, 0, 0, 0, 0) for e in envelope]


def scale_theoretical_isotopic_pattern(eid, tid):
    total = sum(p.intensity for p in eid)
    for p in tid:
        p.intensity *= total


def isotopic_consistency(eic, averagine=glycan, truncate_after=0.95):
    peak_scores = []
    peak_abundances = []
    for node in eic:
        for peak in node.members:
            eid = envelope_to_peak_list(peak.envelope)
            tid = averagine.isotopic_cluster(peak.mz, peak.charge, truncate_after=truncate_after)
            tid.scale(eid)
            peak_scores.append(abs(g_test_scaled(None, eid, tid.truncated_tid)))
            peak_abundances.append(peak.intensity)
    return max(1 - np.average(peak_scores, weights=peak_abundances), 1e-4)


def spacing_fit(eic):
    times, intensities = eic.as_arrays()
    last_rt = times[0]
    last_int = intensities[0]
    rt_deltas = []
    intensity_deltas = []
    for rt, inten in zip(times[1:], intensities[1:]):
        d_rt = rt - last_rt
        rt_deltas.append(d_rt)
        intensity_deltas.append(abs(last_int - inten))
        last_rt = rt
        last_int = inten
    return max(1 - np.average(rt_deltas, weights=intensity_deltas) * 2, 1e-4)


def shape_fit(eic, smooth=0.15):
    return max(1 - AdaptiveMultimodalChromatogramShapeFitter(eic, smooth=smooth).line_test, 1e-4)


def profile_qc(eic, smooth=0.15, averagine=glycan, truncate_after=0.95):
    v = 1.0
    v *= isotopic_consistency(eic, averagine, truncate_after)
    v *= spacing_fit(eic)
    v *= shape_fit(eic, smooth)
    return v


try:
    has_c = True
    _map_coord = map_coord
    _LCMSFeatureSetFit = LCMSFeatureSetFit
    from ms_deisotope._c.feature_map.feature_fit import (LCMSFeatureSetFit, map_coord)
except ImportError as e:
    print(e)
    has_c = False
