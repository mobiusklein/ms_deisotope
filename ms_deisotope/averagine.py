from brainpy import calculate_mass, neutral_mass, PROTON, isotopic_variants, mass_charge_ratio


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

    def __repr__(self):
        return "Averagine(%r)" % self.base_composition


peptide = Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})


class MultipleAveragine(object):
    def __init__(self, compositions):
        self.averagines = list(map(Averagine, compositions))

    def add(self, averagine):
        self.averagines.append(averagine)

    def scale(self, mz, charge=1, charge_carrier=PROTON):
        for case in self.averagines:
            scaled = case.scale(mz, charge=charge, charge_carrier=charge_carrier)
            yield case, scaled

    def isotopic_cluster(self, mz, charge=1, charge_carrier=PROTON, truncate_after=0.9):
        for case, composition in self.scale(mz, charge=charge, charge_carrier=charge_carrier):
            cumsum = 0
            result = []
            for peak in isotopic_variants(composition, charge=charge):
                cumsum += peak.intensity
                result.append(peak)
                if cumsum >= truncate_after:
                    break
            yield case, result

multiple = MultipleAveragine([
    {"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417},
    {"C": 2.9384, "H": 7.7583, "N": 4.3577, "O": 0.4773, "S": 0.0417},
    {"C": 2.9384, "H": 7.7583, "N": 0.3577, "O": 4.4773, "S": 0.0417},
    {"C": 2.9384, "H": 7.7583, "N": 1.3577, "O": 0.4773, "S": 2.0417},
    {"C": 6.9384, "H": 7.7583, "N": 0.3577, "O": 0.4773, "S": 0.0417},
    {"C": 8.}
    ])
