from .cv import Term


class ActivationInformation(object):

    def __init__(self, method, energy, data=None):
        if data is None:
            data = dict()
        self.method = dissociation_methods_map.get(str(method).lower(), method)
        self.energy = energy
        self.data = data

    def __repr__(self):
        return "ActivationInformation(%r, %r%s)" % (
            str(self.method), self.energy,
            "" if not self.data else ", %r" % self.data)

    def __str__(self):
        return str(self.method)

    def is_multiple_dissociation(self):
        return False

    def has_supplemental_dissociation(self):
        return False

    def __eq__(self, other):
        if other is None:
            return False
        if self.is_multiple_dissociation() != other.is_multiple_dissociation():
            return False
        return (self.method == other.method) and abs(self.energy - other.energy) < 1e-3

    def __ne__(self, other):
        return not (self == other)


class MultipleActivationInformation(ActivationInformation):

    def __init__(self, methods, energies, data=None):
        if data is None:
            data = dict()
        self.methods = []
        for method in methods:
            self.methods.append(
                dissociation_methods_map.get(
                    str(method).lower(), method))
        self.energies = list(energies)
        self.data = data

    @property
    def method(self):
        return self.methods[0]

    @property
    def energy(self):
        return self.energies[0]

    def __repr__(self):
        return "MultipleActivationInformation(methods=%r, energies=%r%s)" % (
            list(map(str, self.methods)), self.energies,
            "" if not self.data else ", %r" % self.data)

    def __str__(self):
        return ', '.join(map(str, self.methods))

    def is_multiple_dissociation(self):
        return True

    def has_supplemental_dissociation(self):
        return any([d.is_supplemental() for d in self.methods])

    def __eq__(self, other):
        if other is None:
            return False
        if self.is_multiple_dissociation() != other.is_multiple_dissociation():
            return False
        return (self.methods == other.methods) and all(
            abs(self_energy - other_energy) < 1e-3 for self_energy, other_energy in zip(
                self.energies, other.energies))


class DissociationMethod(Term):

    def is_supplemental(self):
        return "supplemental" in self.name


dissociation_methods = [
    DissociationMethod(u'sustained off-resonance irradiation', u'MS:1000282',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'post-source decay', u'MS:1000135',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'plasma desorption', u'MS:1000134',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'surface-induced dissociation', u'MS:1000136',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'collision-induced dissociation', u'MS:1000133',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'pulsed q dissociation', u'MS:1000599',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'electron transfer dissociation', u'MS:1000598',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'in-source collision-induced dissociation', u'MS:1001880',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'infrared multiphoton dissociation', u'MS:1000262',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'blackbody infrared radiative dissociation', u'MS:1000242',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'low-energy collision-induced dissociation', u'MS:1000433',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'photodissociation', u'MS:1000435',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'LIFT', u'MS:1002000', 'dissociation method',
                       [u'dissociation method']),
    DissociationMethod(u'Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)',
                       u'MS:1002631', 'dissociation method', [u'dissociation method']),
    DissociationMethod(u'electron capture dissociation', u'MS:1000250',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'trap-type collision-induced dissociation', u'MS:1002472',
                       'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'beam-type collision-induced dissociation', u'MS:1000422',
                       'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'supplemental collision-induced dissociation', u'MS:1002679',
                       'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'higher energy beam-type collision-induced dissociation', u'MS:1002481', 'dissociation method',
                       [u'beam-type collision-induced dissociation',
                        u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'supplemental beam-type collision-induced dissociation', u'MS:1002678', 'dissociation method',
                       [u'beam-type collision-induced dissociation',
                        u'collision-induced dissociation', u'dissociation method']),
]


UnknownDissociation = DissociationMethod(
    "unknown dissociation", None, 'dissociation method',
    [u'dissociation method'])


dissociation_methods_map = {
    "": UnknownDissociation,
    None: UnknownDissociation,
    UnknownDissociation.name: UnknownDissociation
}


for method in dissociation_methods:
    dissociation_methods_map[method.name] = method

CID = dissociation_methods_map.get("collision-induced dissociation")
HCD = dissociation_methods_map.get("beam-type collision-induced dissociation")
ETD = dissociation_methods_map.get("electron transfer dissociation")
ECD = dissociation_methods_map.get("electron capture dissociation")

supplemental_term_map = {
    dissociation_methods_map["beam-type collision-induced dissociation"]: dissociation_methods_map[
        "supplemental beam-type collision-induced dissociation"],
    dissociation_methods_map["collision-induced dissociation"]: dissociation_methods_map[
        "supplemental collision-induced dissociation"]
}


supplemental_energy = "supplemental collision energy"


energy_terms = set([
    "collision energy",
    supplemental_energy,
    "activation energy",
    "collision energy ramp start",
    "collision energy ramp end",
    "percent collision energy ramp start",
    "percent collision energy ramp end",
])


dissociation_methods_map.update({
    'cad': CID,
    'cid': CID,
    "hcd": HCD,
    "etd": ETD,
    "ecd": ECD
})


ActivationInformation.dissociation_methods = dissociation_methods_map
