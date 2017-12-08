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


dissociation_methods = [
    Term(u'sustained off-resonance irradiation', u'MS:1000282',
         'dissociation method', [u'dissociation method']),
    Term(u'post-source decay', u'MS:1000135',
         'dissociation method', [u'dissociation method']),
    Term(u'plasma desorption', u'MS:1000134',
         'dissociation method', [u'dissociation method']),
    Term(u'surface-induced dissociation', u'MS:1000136',
         'dissociation method', [u'dissociation method']),
    Term(u'collision-induced dissociation', u'MS:1000133',
         'dissociation method', [u'dissociation method']),
    Term(u'pulsed q dissociation', u'MS:1000599',
         'dissociation method', [u'dissociation method']),
    Term(u'electron transfer dissociation', u'MS:1000598',
         'dissociation method', [u'dissociation method']),
    Term(u'in-source collision-induced dissociation', u'MS:1001880',
         'dissociation method', [u'dissociation method']),
    Term(u'infrared multiphoton dissociation', u'MS:1000262',
         'dissociation method', [u'dissociation method']),
    Term(u'blackbody infrared radiative dissociation', u'MS:1000242',
         'dissociation method', [u'dissociation method']),
    Term(u'low-energy collision-induced dissociation', u'MS:1000433',
         'dissociation method', [u'dissociation method']),
    Term(u'photodissociation', u'MS:1000435',
         'dissociation method', [u'dissociation method']),
    Term(u'LIFT', u'MS:1002000', 'dissociation method',
         [u'dissociation method']),
    Term(u'Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)',
         u'MS:1002631', 'dissociation method', [u'dissociation method']),
    Term(u'electron capture dissociation', u'MS:1000250',
         'dissociation method', [u'dissociation method']),
    Term(u'trap-type collision-induced dissociation', u'MS:1002472',
         'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    Term(u'beam-type collision-induced dissociation', u'MS:1000422',
         'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    Term(u'supplemental collision-induced dissociation', u'MS:1002679',
         'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    Term(u'higher energy beam-type collision-induced dissociation', u'MS:1002481', 'dissociation method',
         [u'beam-type collision-induced dissociation', u'collision-induced dissociation', u'dissociation method']),
    Term(u'supplemental beam-type collision-induced dissociation', u'MS:1002678', 'dissociation method',
         [u'beam-type collision-induced dissociation', u'collision-induced dissociation', u'dissociation method']),
]


UnknownDissociation = Term(
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


dissociation_methods_map.update({
    'cad': CID,
    'cid': CID,
    "hcd": HCD,
    "etd": ETD,
    "ecd": ECD
})


ActivationInformation.dissociation_methods = dissociation_methods_map
