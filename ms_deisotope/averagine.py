from brainpy import calculate_mass, neutral_mass, PROTON, isotopic_variants


class Averagine(object):
    def __init__(self, base_composition):
        self.base_composition = base_composition
        self.base_mass = calculate_mass(base_composition)

    def scale(self, mz, charge=1, charge_carrier=PROTON):
        neutral = neutral_mass(mz, charge, charge_carrier)

        scale = neutral / self.base_mass
        scaled = {}
        for elem, count in self.base_composition.items():
            scaled[elem] = count * scale
        return scaled

    def isotopic_cluster(self, mz, charge=1, charge_carrier=PROTON, truncate_after=0.9):
        composition = self.scale(mz, charge, charge_carrier)
        cumsum = 0
        result = []
        for peak in isotopic_variants(composition, charge=charge):
            cumsum += peak.intensity
            result.append(peak)
            if cumsum >= truncate_after:
                break
        return result


peptide = Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})
