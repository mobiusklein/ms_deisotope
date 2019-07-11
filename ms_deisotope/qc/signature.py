from ms_deisotope.data_source.scan.base import BasePeakMethods


class SignatureIonScanner(object):
    def __init__(self, signature_ions=None):
        if signature_ions is None:
            signature_ions = []
        self.signature_ions = signature_ions

    def scan(self, peak_list, error_tolerance=2e-5):
        matches = []
        for ion in self.signature_ions:
            match = peak_list.has_peak(
                ion, error_tolerance)
            if match is not None:
                matches.append(match)
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

    def __call__(self, peak_list, error_tolerance=2e-5):
        return self.ratio(peak_list, error_tolerance)
