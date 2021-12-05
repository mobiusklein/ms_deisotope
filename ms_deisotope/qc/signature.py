'''A module for signature ion extraction.
'''

from collections import namedtuple
from numbers import Number

from ms_deisotope.utils import Base
from ms_deisotope.averagine import mass_charge_ratio
from ms_deisotope.data_source.scan.base import BasePeakMethods, PeakSetMethods

class Target(Base):
    """An analyte with known mass and expected charge states

    Attributes
    ----------
    name: str
        The name of the analyte
    neutral_mass: float
        The known neutral mass of the analyte
    charge_states: list
        The list of expected charge states, defaults to ``[1]``

    """
    def __init__(self, name, neutral_mass, charge_states=None):
        if isinstance(charge_states, Number):
            charge_states = [charge_states]
        elif charge_states is None:
            charge_states = [1]
        self.name = name
        self.neutral_mass = neutral_mass
        self.charge_states = charge_states

    def iter_mz(self):
        """Iterate over theoretical m/z coordinates for this analyte

        Yields
        ------
        float:
            The theoretical m/z for a charge state
        """
        for charge in self.charge_states:
            yield mass_charge_ratio(self.neutral_mass, charge)

    def __eq__(self, other):
        return self.name == other.name and self.charge_states == other.charge_states and abs(
            self.neutral_mass - other.neutral_mass) < 1e-3

    def __hash__(self):
        return hash(self.name)


class SignatureIonDetector(object):
    """Detect when a scan contains intense signal for signature ions.

    Implements the :class:`~.Callable` interface.

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
        '''Scan the peak list for signature ion masses or m/zs.

        Parameters
        ----------
        peak_list : peak set-like
            The peak set to search. May be any object compatible with
            :class:`~.PeakSetMethods`.
        error_tolerance : float
            The PPM error tolerance for peak mass coordinate matching.

        Returns
        -------
        list
        '''
        result = []
        peak_list = PeakSetMethods(peak_list)
        if peak_list.is_deconvoluted:
            for sig in self.signature_ions:
                matches = peak_list.all_peaks_for(
                    sig.neutral_mass, error_tolerance)
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
        '''Compute an intensity ratio of the signature ions against the
        base peak.

        Parameters
        ----------
        peak_list : peak set-like
            The peak set to search. May be any object compatible with
            :class:`~.PeakSetMethods`.
        error_tolerance : float
            The PPM error tolerance for peak mass coordinate matching.

        Returns
        -------
        float
        '''
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
        '''Test the peak list for the signature ions.

        Parameters
        ----------
        peak_list : peak set-like
            The peak set to search. May be any object compatible with
            :class:`~.PeakSetMethods`.
        error_tolerance : float
            The PPM error tolerance for peak mass coordinate matching.

        Returns
        -------
        float
        '''
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
]


TMT11_INFO = TMT10_INFO[:]
TMT11_INFO.append(
    TMTInfo('TMT11-131C', 229.162932, 130.13722253323, 117.13386453323001))


TMT16_INFO = [
    TMTInfo('TMTpro-126', 304.207100, 126.127726, 0.0),
    TMTInfo('TMTpro-127N', 304.207100, 127.124761, 0.0),
    TMTInfo('TMTpro-127C', 304.207100, 127.131081, 0.0),
    TMTInfo('TMTpro-128N', 304.207100, 128.128116, 0.0),
    TMTInfo('TMTpro-128C', 304.207100, 128.134436, 0.0),
    TMTInfo('TMTpro-129N', 304.207100, 129.131471, 0.0),
    TMTInfo('TMTpro-129C', 304.207100, 129.137790, 0.0),
    TMTInfo('TMTpro-130N', 304.207100, 130.134825, 0.0),
    TMTInfo('TMTpro-130C', 304.207100, 130.141145, 0.0),
    TMTInfo('TMTpro-131N', 304.207100, 131.138180, 0.0),
    TMTInfo('TMTpro-131C', 304.207100, 131.144500, 0.0),
    TMTInfo('TMTpro-132N', 304.207100, 132.141535, 0.0),
    TMTInfo('TMTpro-132C', 304.207100, 132.147855, 0.0),
    TMTInfo('TMTpro-133N', 304.207100, 133.144890, 0.0),
    TMTInfo('TMTpro-133C', 304.207100, 133.151210, 0.0),
    TMTInfo('TMTpro-134N', 304.207100, 134.148245, 0.0)
]


class SignatureIonExtractor(SignatureIonDetector):
    '''Extracts signal for a set of target ions from each scan.
    '''
    def extract(self, peak_list, error_tolerance=2e-5):
        '''Extract the intensity for each signature ion.

        Parameters
        ----------
        peak_list : peak set-like
            The peak set to search. May be any object compatible with
            :class:`~.PeakSetMethods`.
        error_tolerance : float
            The PPM error tolerance for peak mass coordinate matching.

        Returns
        -------
        channels: dict[str, float]
            A mapping from signature ion name to its intensity.
        '''
        result = {}
        peak_list = PeakSetMethods(peak_list)
        if peak_list.is_deconvoluted:
            for sig in self.signature_ions:
                result[sig.name] = 0
                matches = peak_list.all_peaks_for(sig.neutral_mass, error_tolerance)
                for match in matches:
                    if match.charge in sig.charge_states:
                        result[sig.name] += match.intensity
        else:
            for sig in self.signature_ions:
                result[sig.name] = 0
                for mz in sig.iter_mz():
                    match, error = peak_list.get_nearest_peak(mz)
                    if abs(error) / mz < error_tolerance:
                        result[sig.name] += match.intensity
        return result


    def __call__(self, peak_list, error_tolerance=2e-5):
        return self.extract(peak_list, error_tolerance)


class TMTReporterExtractor(SignatureIonExtractor):
    '''A :class:`~.SignatureIonExtractor` for TMT reporter ions.
    '''
    TMT_REAGENTS = {
        "tmt16": TMT16_INFO,
        "tmt11": TMT11_INFO,
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

    def __init__(self, reagent="TMT11plex"):
        reagents = self._find_reagent(reagent)
        super(TMTReporterExtractor, self).__init__(
            [Target(
                n.name,
                n.reporter_mass_hcd)
                for n in reagents])
