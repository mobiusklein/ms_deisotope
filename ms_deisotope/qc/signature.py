from collections import namedtuple
from numbers import Number

from ms_deisotope.utils import Base
from ms_deisotope.averagine import mass_charge_ratio
from ms_deisotope.data_source.scan.base import BasePeakMethods, PeakSetMethods


class Target(Base):
    def __init__(self, name, neutral_mass, charge_states=None):
        if isinstance(charge_states, Number):
            charge_states = [charge_states]
        elif charge_states is None:
            charge_states = [1]
        self.name = name
        self.neutral_mass = neutral_mass
        self.charge_states = charge_states

    def iter_mz(self):
        for charge in self.charge_states:
            yield mass_charge_ratio(self.neutral_mass, charge)

    def __eq__(self, other):
        return self.name == other.name and self.charge_states == other.charge_states and abs(
            self.neutral_mass - other.neutral_mass) < 1e-3

    def __hash__(self):
        return hash(self.name)


class SignatureIonDetector(object):
    """Detect when a scan contains intense signal for signature ions.

    Attributes
    ----------
    signature_ions: :class:`list` of :class:`Target`
        The target ions to scan for.
    method: str
        The detection method name, as a string. Currently supports "ratio" only.

    """

    METHODS = ('ratio', )

    def __init__(self, signature_ions=None, method='ratio'):
        if signature_ions is None:
            signature_ions = []
        self.signature_ions = signature_ions
        if method not in self.METHODS:
            raise ValueError("Method not recognized: %r" % (method, ))
        self.method = method

    def scan(self, peak_list, error_tolerance=2e-5):
        result = []
        peak_list = PeakSetMethods(peak_list)
        if peak_list.is_deconvoluted:
            for sig in self.signature_ions:
                matches = peak_list.all_peaks_for(sig.mass, error_tolerance)
                for match in matches:
                    if match.charge in sig.charge_states:
                        result.append((sig, match))
        else:
            for sig in self.signature_ions:
                for mz in sig.iter_mz():
                    matches = peak_list.all_peaks_for(mz, error_tolerance)
                    for match in matches:
                        result.append((sig, match))
        return matches

    def ratio(self, peak_list, error_tolerance=2e-5):
        try:
            base_peak = BasePeakMethods(peak_list)()
            maximum = base_peak.intensity
        except (ValueError, TypeError):
            if len(peak_list) == 0:
                return 0
            else:
                raise
        n = len([i for i in self.signature_ions])
        if n == 0:
            return 0
        oxonium = sum(
            p.intensity / maximum for p in self.scan(
                peak_list, error_tolerance))
        return oxonium / n

    def detect(self, peak_list, error_tolerance=2e-5):
        if self.method == 'ratio':
            return self.ratio(peak_list, error_tolerance)
        else:
            raise NotImplementedError(self.method)

    def __call__(self, peak_list, error_tolerance=2e-5):
        return self.detect(peak_list, error_tolerance)

TMTInfo = namedtuple("TMTInfo", ["name", "intact_mass", "reporter_mass_hcd", "reporter_mass_etd"])


TMT0_INFO = [
    TMTInfo('TMT0-126', 224.152478, 125.12044953323, 113.12044853323)
]


TMT2_INFO = [
    TMTInfo('TMT2-126', 225.155833, 125.12044953323, 113.12044853323),
    TMTInfo('TMT2-127', 225.155833, 126.12380453323, 113.12044853323),
]


TMT6_INFO = [
    TMTInfo('TMT6-126', 229.162932, 125.12044953323, 113.12044853323),
    TMTInfo('TMT6-127', 229.162932, 126.11748453323001, 114.11748353323),
    TMTInfo('TMT6-128', 229.162932, 127.12715953323, 115.12715653323),
    TMTInfo('TMT6-129', 229.162932, 128.12419453323, 116.12419153323),
    TMTInfo('TMT6-130', 229.162932, 129.13386853323, 117.13386453323001),
    TMTInfo('TMT6-131', 229.162932, 130.13090353323, 118.13089953323001),
]


TMT10_INFO = [
    TMTInfo('TMT10-126', 229.162932, 125.12044953323, 113.12044853323),
    TMTInfo('TMT10-127N', 229.162932, 126.11748453323001, 114.11748353323),
    TMTInfo('TMT10-127C', 229.162932, 126.12380453323, 113.12044853323),
    TMTInfo('TMT10-128N', 229.162932, 127.12083953323001, 114.12044853323),
    TMTInfo('TMT10-128C', 229.162932, 127.12715953323, 115.12715653323),
    TMTInfo('TMT10-129N', 229.162932, 128.12419453323, 116.12419153323),
    TMTInfo('TMT10-129C', 229.162932, 128.13051353323, 115.12715653323),
    TMTInfo('TMT10-130N', 229.162932, 129.12754853323, 116.12419153323),
    TMTInfo('TMT10-130C', 229.162932, 129.13386853323, 117.13386453323001),
    TMTInfo('TMT10-131', 229.162932, 130.13090353323, 118.13089953323001),
    TMTInfo('TMT11-131C', 229.162932, 130.13722253323, 117.13386453323001),
]


class SignatureIonExtractor(SignatureIonDetector):

    def extract(self, peak_list, error_tolerance=2e-5):
        result = {}
        peak_list = PeakSetMethods(peak_list)
        if peak_list.is_deconvoluted:
            for sig in self.signature_ions:
                result[sig.name] = 0
                matches = peak_list.all_peaks_for(sig.mass, error_tolerance)
                for match in matches:
                    if match.charge in sig.charge_states:
                        result[sig.name] += match.intensity
        else:
            for sig in self.signature_ions:
                result[sig.name] = 0
                for mz in sig.iter_mz():
                    matches = peak_list.all_peaks_for(mz, error_tolerance)
                    for match in matches:
                        result[sig.name] += match.intensity
        return result


    def __call__(self, peak_list, error_tolerance=2e-5):
        return self.extract(peak_list, error_tolerance)


class TMTReporterExtractor(SignatureIonExtractor):
    TMT_REAGENTS = {
        "tmt10": TMT10_INFO,
        "tmt2": TMT2_INFO,
        "tmt6": TMT6_INFO,
        "tmt0": TMT0_INFO
    }

    def _find_reagent(self, reagent):
        label = reagent.lower().replace("plex", "")
        if label in self.TMT_REAGENTS:
            reagents = self.TMT_REAGENTS[label]
        else:
            raise KeyError("Did not recognize reagent %s" % (reagent, ))
        return reagents

    def __init__(self, reagent="TMT10plex"):
        reagents = self._find_reagent(reagent)
        super(TMTReporterExtractor, self).__init__(
            [Target(
                n.name,
                n.reporter_mass_hcd)
                for n in reagents])
