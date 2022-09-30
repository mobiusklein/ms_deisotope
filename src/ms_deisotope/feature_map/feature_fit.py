from collections import namedtuple

import numpy as np

# from brainpy import neutral_mass as calc_neutral_mass
from ms_peak_picker import FittedPeak

from ms_deisotope.averagine import glycan, mass_charge_ratio
from ms_deisotope.scoring import g_test_scaled
from ms_deisotope.peak_set import DeconvolutedPeak, IonMobilityDeconvolutedPeak, IonMobilityProfileDeconvolutedPeakSolution, Envelope

from ms_deisotope.peak_dependency_network.intervals import SimpleInterval

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
    def neutral_mass(self) -> float:
        if self._neutral_mass == 0:
            if self._most_abundant_member is not None:
                self._neutral_mass = self._most_abundant_member.neutral_mass
        return self._neutral_mass


class DeconvolutedLCMSFeature(LCMSFeature):
    __slots__ = ('score', 'n_features', 'supporters',
                 '_neutral_mass', '_last_neutral_mass',
                 '_precursor_information', 'charge', '_mz')

    def __init__(self, nodes=None, charge=None, adducts=None, used_as_adduct=None, score=0.0,
                 n_features=0, feature_id=None, supporters=None):
        if supporters is None:
            supporters = []
        self.charge = charge
        self.score = score
        self._neutral_mass = None
        self._last_neutral_mass = None
        self._precursor_information = None
        self._mz = None
        self.n_features = n_features
        self.supporters = supporters
        super(DeconvolutedLCMSFeature, self).__init__(nodes, adducts, used_as_adduct, feature_id=feature_id)

    def _initialize_averager(self):
        self._peak_averager = DeconvolutedRunningWeightedAverage()

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

    def clone(self, deep=False, cls=None) -> 'DeconvolutedLCMSFeature':
        if cls is None:
            cls = self.__class__
        return cls(
            self.nodes.clone(deep=deep), self.charge, self.adducts, self.used_as_adduct, self.score,
            self.n_features, self.feature_id, list(self.supporters))

    def _invalidate(self, reaverage: bool=False):
        self._last_neutral_mass = self._neutral_mass if self._neutral_mass is not None else 0.
        self._neutral_mass = None
        self._precursor_information = None
        self._mz = None
        super(DeconvolutedLCMSFeature, self)._invalidate(reaverage)

    def _update_from_averager(self):
        self._neutral_mass = self._last_neutral_mass = self._peak_averager.current_mean

    @property
    def neutral_mass(self) -> float:
        if self._neutral_mass is None:
            self._reaverage_and_update()
        return self._neutral_mass

    @property
    def mz(self):
        if self._mz is None:
            self._mz = mass_charge_ratio(self.neutral_mass, self.charge)
        return self._mz

    def _copy_chunk(self, nodes, *args, **kwargs) -> 'DeconvolutedLCMSFeature':
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

    def sum_envelopes(self) -> Envelope:
        mzs = None
        intensities = None
        for node in self:
            peak: DeconvolutedPeak
            for peak in node.members:
                if mzs is None:
                    mzs = [0.0 for i in range(len(peak.envelope))]
                    intensities = [0.0 for i in range(len(peak.envelope))]
                for i, pt in enumerate(peak.envelope):
                    try:
                        mzs[i] += (pt.mz * pt.intensity)
                        intensities[i] += pt.intensity
                    except IndexError:
                        mzs.append(pt.mz * pt.intensity)
                        intensities.append(pt.intensity)
        out = [None for i in range(len(mzs))]
        for i, mz_int in enumerate(mzs):
            inten = intensities[i]
            out[i] = (mz_int / inten, inten)
        return Envelope(out)


class DeconvolutedRunningWeightedAverage(RunningWeightedAverage):
    def add(self, peak: DeconvolutedPeak):
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

try:
    from ms_deisotope._c.feature_map.lcms_feature import RunningWeightedAverageNeutralMass as DeconvolutedRunningWeightedAverage
except ImportError:
    pass


class DeconvolutedDriftTimeRunningWeightedAverage(DeconvolutedRunningWeightedAverage):
    def _initialize(self):
        self.neutral_mass_mean = 0.0
        self.drift_time_mean = 0.0

    def add(self, peak: IonMobilityDeconvolutedPeak):
        if peak.intensity == 0:
            if self.current_mean == 0 and self.total_weight == 0:
                self.neutral_mass_mean = peak.neutral_mass
                self.drift_time_mean = peak.drift_time
                self.total_weight = 1
            else:
                return
        drift_time_agg = (self.total_weight * self.drift_time_mean) + \
            (peak.drift_time * peak.intensity)
        neutral_mass_agg = (self.total_weight * self.neutral_mass_mean) + \
            (peak.neutral_mass * peak.intensity)
        self.total_weight += peak.intensity
        self.drift_time_mean = drift_time_agg / self.total_weight
        self.neutral_mass_mean = neutral_mass_agg / self.total_weight
        return self


class IonMobilityDeconvolutedLCMSFeature(DeconvolutedLCMSFeature):
    def __init__(self, nodes=None, charge=None, adducts=None, used_as_adduct=None, score=0.0,
                 n_features=0, feature_id=None, supporters=None):
        self._drift_time = None
        self._last_drift_time = None
        super(IonMobilityDeconvolutedLCMSFeature, self).__init__(
            nodes=nodes, charge=charge, adducts=adducts, used_as_adduct=used_as_adduct, score=score,
            n_features=n_features, feature_id=feature_id, supporters=supporters)

    def _initialize_averager(self):
        self._peak_averager = DeconvolutedDriftTimeRunningWeightedAverage()

    def _invalidate(self, reaverage: bool=False):
        self._last_drift_time = self._drift_time if self._drift_time is not None else 0.
        self._drift_time = None
        return super(IonMobilityDeconvolutedLCMSFeature, self)._invalidate(reaverage=reaverage)

    def _update_from_averager(self):
        self._neutral_mass = self._last_neutral_mass = self._peak_averager.neutral_mass_mean
        self._drift_time = self._last_drift_time = self._peak_averager.drift_time_mean

    @property
    def drift_time(self) -> float:
        if self._drift_time is None:
            self._reaverage_and_update()
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


class IonMobilityProfileDeconvolutedRunningWeightedAverage(DeconvolutedRunningWeightedAverage):
    def __init__(self, *args, **kwargs):
        self._initialize()

    def _initialize(self):
        super()._initialize()
        self.interval = SimpleInterval(0, 0)

    def add(self, peak: IonMobilityProfileDeconvolutedPeakSolution):
        if peak.intensity == 0:
            if self.current_mean == 0 and self.total_weight == 0:
                self.current_mean = peak.neutral_mass
                self.interval.start = peak.ion_mobility_interval.start
                self.interval.end = peak.ion_mobility_interval.end
                self.total_weight = 1
            else:
                return
        agg = (self.total_weight * self.current_mean) + \
            (peak.neutral_mass * peak.intensity)
        start_agg = (self.interval.start * self.total_weight) + (peak.ion_mobility_interval.start * peak.intensity)
        end_agg = (self.interval.end * self.total_weight) + (peak.ion_mobility_interval.end * peak.intensity)
        self.total_weight += peak.intensity
        self.current_mean = agg / self.total_weight
        self.interval.start = start_agg / self.total_weight
        self.interval.end = end_agg / self.total_weight
        self.current_count += 1
        return self


class IonMobilityProfileDeconvolutedLCMSFeature(DeconvolutedLCMSFeature):
    __slots__ = ('_ion_mobility_interval', '_last_ion_mobility_interval')

    _ion_mobility_interval: SimpleInterval

    def __init__(self, nodes=None, charge=None, adducts=None, used_as_adduct=None, score=0.0,
                 n_features=0, feature_id=None, supporters=None):
        super().__init__(nodes, charge, adducts, used_as_adduct, score,
                         n_features, feature_id, supporters)
        self._ion_mobility_interval = None
        self._last_ion_mobility_interval = None

    def _invalidate(self, reaverage: bool=False):
        self._last_ion_mobility_interval = self._ion_mobility_interval if self._ion_mobility_interval is not None \
            else SimpleInterval(0, 0)
        self._ion_mobility_interval = None
        super()._invalidate(reaverage)

    def _initialize_averager(self):
        self._peak_averager = IonMobilityProfileDeconvolutedRunningWeightedAverage()

    def _update_from_averager(self):
        self._neutral_mass = self._last_neutral_mass = self._peak_averager.current_mean
        self._ion_mobility_interval = self._last_ion_mobility_interval = self._peak_averager.interval

    @property
    def neutral_mass(self) -> float:
        if self._neutral_mass is None:
            self._reaverage_and_update()
        return self._neutral_mass

    @property
    def ion_mobility_interval(self) -> SimpleInterval:
        if self._ion_mobility_interval is None:
            self._reaverage_and_update()
        return self._ion_mobility_interval

    def between_ion_mobilities(self, ion_mobility_start: float, ion_mobility_end: float,
                               time_start: float = 0, time_end: float = float('inf'),
                               include_mz: bool=True, include_envelope: bool=True) -> 'IonMobilityProfileDeconvolutedLCMSFeature':
        nodes = []
        node, i = self.find_time(time_start)
        if i > 0:
            i -= 1
        for node in self:
            if node.time < time_start:
                continue
            if node.time > time_end:
                break
            peaks = []
            for peak in node.members:
                try:
                    peaks.append(
                        peak.from_feature(
                            peak.solution.between_times(ion_mobility_start, ion_mobility_end),
                            include_mz=include_mz, include_envelope=include_envelope))
                except IndexError:
                    continue
            if peaks:
                node = node.__class__(node.time, peaks)
                nodes.append(node)
        return self._copy_chunk(nodes)


try:
    has_c = True
    from ms_deisotope._c.feature_map.feature_fit import (
        _sum_envelopes as _c_sum_envelopes)

    DeconvolutedLCMSFeature.sum_envelopes = _c_sum_envelopes
except ImportError as e:
    print(e)
    has_c = False
