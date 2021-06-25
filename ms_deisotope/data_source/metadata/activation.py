from ms_deisotope.utils import _MappingOverAttributeProxy

from .cv import Term, TermSet


class ActivationInformation(object):
    """Describes the dissociation process used to produce an MSn scan.

    Attributes
    ----------
    data : :class:`dict`
        Additional metadata describing the dissociation
    energy : float
        The quantity of energy used
    method : :class:`DissociationMethod`
        The dissociation method used
    """
    __slots__ = ('energy', 'method', 'data')

    def __init__(self, method, energy, data=None):  # pylint: disable=redefined-outer-name
        if data is None:
            data = dict()
        self.method = dissociation_methods_map.get(str(method).lower(), None)
        if self.method is None:
            try:
                self.method = dissociation_methods_map.get(
                    method.accession, None)
            except AttributeError:
                pass
            if self.method is None:
                self.method = self._make_unknown_method(method)
        self.energy = energy
        self.data = data

    def to_dict(self):
        d = {
            "method": str(self.method),
            "energy": self.energy,
            "data": self.data
        }
        return d

    @classmethod
    def from_dict(cls, d):
        if 'methods' in d and 'energies' in d:
            return MultipleActivationInformation(**d)
        else:
            return ActivationInformation(**d)

    @staticmethod
    def _make_unknown_method(method):  # pylint: disable=redefined-outer-name
        return DissociationMethod(UnknownDissociation.name, method, method, method, [method])

    def __repr__(self):
        return "ActivationInformation(%r, %r%s)" % (
            str(self.method), self.energy,
            "" if not self.data else ", %r" % self.data)

    def __str__(self):
        return str(self.method)

    def __hash__(self):
        return hash(str(self))

    def is_multiple_dissociation(self):
        """Determine if multiple dissociation methods were
        used.

        Returns
        -------
        bool
        """
        return False

    def has_supplemental_dissociation(self):
        """Determine if supplemental activation was used

        Returns
        -------
        bool
        """
        return False

    def __eq__(self, other):
        if other is None:
            return False
        if self.is_multiple_dissociation() != other.is_multiple_dissociation():
            return False
        if self.energy is not None:
            if other.energy is not None:
                energy_match = abs(self.energy - other.energy) < 1e-3
            else:
                energy_match = False
        else:
            if other.energy is not None:
                energy_match = False
            else:
                energy_match = True
        return (self.method == other.method) and energy_match

    def __ne__(self, other):
        return not self == other

    def has_dissociation_type(self, dissociation_type):
        return self.method.is_a(dissociation_type)

    @property
    def __dict__(self):
        return _MappingOverAttributeProxy(self)

    def __reduce__(self):
        return self.__class__, (self.method, self.energy, self.data)


class MultipleActivationInformation(ActivationInformation):
    """Describes the dissociation processes used to produce an MSn
    scan with multiple types of activation.

    Attributes
    ----------
    data : :class:`dict`
        Additional metadata describing the dissociation
    energies : :class:`list` of float
        A list of dissociation energies used. Parallel to :attr:`methods`
    methods : :class:`list` of :class:`DissociationMethod`
        A list of dissociation methods used. Parallel to :attr:`energies`
    """

    def __init__(self, methods, energies, data=None):  # pylint: disable=super-init-not-called
        if data is None:
            data = dict()
        self.methods = []
        for method in methods:  # pylint: disable=redefined-outer-name
            self.methods.append(
                dissociation_methods_map.get(
                    str(method).lower(), method))
        self.energies = list(energies)
        self.data = data

    def to_dict(self):
        d = {
            "methods": [str(m) for m in self.methods],
            "energies": self.energies,
            "data": self.data
        }
        return d

    @property
    def method(self):
        return self.methods[0]

    @property
    def energy(self):
        return self.energies[0]

    def __iter__(self):
        return iter(self.methods)

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
        if len(self.energies) != len(other.energies):
            return False
        energies_equal = []
        for se, oe in zip(self.energies, other.energies):
            if se is None and oe is None:
                energies_equal.append(True)
            elif se is None or oe is None:
                energies_equal.append(False)
            else:
                energies_equal.append(abs(se - oe) < 1e-3)
        return (self.methods == other.methods) and all(energies_equal)

    def has_dissociation_type(self, dissociation_type):
        return any(method.is_a(dissociation_type) for method in self.methods)

    def __reduce__(self):
        return self.__class__, (self.methods, self.energies, self.data)


class DissociationMethod(Term):
    """Controlled vocabulary describing a dissociation process
    """

    def is_supplemental(self):
        """Returns whether or not the dissociation process is
        supplemental

        Returns
        -------
        bool
        """
        return "supplemental" in self.name


dissociation_methods = []


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('dissociation method', term_cls_name="DissociationMethod", writer=cog.out)
# ]]]
# CV Version: 4.1.55
dissociation_methods = TermSet([
    DissociationMethod('collision-induced dissociation', 'MS:1000133',
                       ('The dissociation of an ion after collisional excitation. The '
                        'term collisional-activated dissociation is not recommended.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('plasma desorption', 'MS:1000134',
                       ('The ionization of material in a solid sample by bombarding '
                        'it with ionic or neutral atoms formed as a result of the '
                        'fission of a suitable nuclide, typically 252Cf. Synonymous '
                        'with fission fragment ionization.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('post-source decay', 'MS:1000135',
                       ('A technique specific to reflectron time-of-flight mass '
                        'spectrometers where product ions of metastable transitions '
                        'or collision-induced dissociations generated in the drift '
                        'tube prior to entering the reflectron are m/z separated to '
                        'yield product ion spectra.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('surface-induced dissociation', 'MS:1000136',
                       ('Fragmentation that results from the collision of an ion with '
                        'a surface.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('blackbody infrared radiative dissociation', 'MS:1000242',
                       ('A special case of infrared multiphoton dissociation wherein '
                        'excitation of the reactant ion is caused by absorption of '
                        'infrared photons radiating from heated blackbody '
                        'surroundings, which are usually the walls of a vacuum '
                        'chamber. See also infrared multiphoton dissociation.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('electron capture dissociation', 'MS:1000250',
                       ('A process in which a multiply protonated molecules interacts '
                        'with a low energy electrons. Capture of the electron leads '
                        'the liberation of energy and a reduction in charge state of '
                        'the ion with the production of the (M + nH) (n-1)+ odd '
                        'electron ion, which readily fragments.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('infrared multiphoton dissociation', 'MS:1000262',
                       ('Multiphoton ionization where the reactant ion dissociates as '
                        'a result of the absorption of multiple infrared photons.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('sustained off-resonance irradiation', 'MS:1000282',
                       ('A technique associated with Fourier transform ion cyclotron '
                        'resonance (FT-ICR) mass spectrometry to carry out '
                        'ion/neutral reactions such as low-energy collision-induced '
                        'dissociation. A radio-frequency electric field of slightly '
                        'off-resonance to the cyclotron frequency of the reactant ion '
                        'cyclically accelerates and decelerates the reactant ion that '
                        "is confined in the Penning ion trap. The ion's orbit does "
                        'not exceed the dimensions of ion trap while the ion '
                        'undergoes an ion/neutral species process that produces a '
                        'high average translational energy for an extended time.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('low-energy collision-induced dissociation', 'MS:1000433',
                       ('A collision-induced dissociation process wherein the '
                        'precursor ion has the translational energy lower than '
                        'approximately 1000 eV. This process typically requires '
                        'multiple collisions and the collisional excitation is '
                        'cumulative.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('photodissociation', 'MS:1000435',
                       ('A process wherein the reactant ion is dissociated as a '
                        'result of absorption of one or more photons.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('electron transfer dissociation', 'MS:1000598',
                       ('A process to fragment ions in a mass spectrometer by '
                        'inducing fragmentation of cations (e.g. peptides or '
                        'proteins) by transferring electrons to them.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('pulsed q dissociation', 'MS:1000599',
                       ('A process that involves precursor ion activation at high Q, '
                        'a time delay to allow the precursor to fragment, then a '
                        'rapid pulse to low Q where all fragment ions are trapped. '
                        'The product ions can then be scanned out of the ion trap and '
                        'detected.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('in-source collision-induced dissociation', 'MS:1001880',
                       ('The dissociation of an ion as a result of collisional '
                        'excitation during ion transfer from an atmospheric pressure '
                        'ion source and the mass spectrometer vacuum.'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('LIFT', 'MS:1002000',
                       ("A Bruker's proprietary technique where molecular ions are "
                        'initially accelerated at lower energy, then collide with '
                        "inert gas in a collision cell that is then 'lifted' to high "
                        'potential. The use of inert gas is optional, as it could '
                        'lift also fragments provided by LID." '
                        '[DOI:10.1007/s00216-003-2057-0'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)', 'MS:1002631',
                       ('A dissociation process combining electron-transfer and '
                        'higher-energy collision dissociation (EThcD). It combines '
                        'ETD (reaction time) followed by HCD (activation energy).'),
                       'dissociation method',
                       ['dissociation method']),
    DissociationMethod('beam-type collision-induced dissociation', 'MS:1000422',
                       ('A collision-induced dissociation process that occurs in a '
                        'beam-type collision cell.'),
                       'dissociation method',
                       ['collision-induced dissociation', 'dissociation method']),
    DissociationMethod('trap-type collision-induced dissociation', 'MS:1002472',
                       ('A collision-induced dissociation process that occurs in a '
                        'trap-type collision cell.'),
                       'dissociation method',
                       ['collision-induced dissociation', 'dissociation method']),
    DissociationMethod('supplemental collision-induced dissociation', 'MS:1002679',
                       ('The dissociation of an ion after supplemental collisional '
                        'excitation.'),
                       'dissociation method',
                       ['collision-induced dissociation', 'dissociation method']),
    DissociationMethod('higher energy beam-type collision-induced dissociation', 'MS:1002481',
                       ('A collision-induced dissociation process wherein the '
                        'projectile ion has the translational energy higher than '
                        'approximately 1000 eV.'),
                       'dissociation method',
                       ['beam-type collision-induced dissociation', 'collision-induced dissociation', 'dissociation method']),
    DissociationMethod('supplemental beam-type collision-induced dissociation', 'MS:1002678',
                       ('A supplemental collision-induced dissociation process that '
                        'occurs in a beam-type collision cell in addition to another '
                        'primary type of dissociation.'),
                       'dissociation method',
                       ['beam-type collision-induced dissociation', 'collision-induced dissociation', 'dissociation method']),
])
# [[[end]]]


UnknownDissociation = DissociationMethod(
    "unknown dissociation", None, None, 'dissociation method',
    [u'dissociation method'])


dissociation_methods_map = {
    "": UnknownDissociation,
    None: UnknownDissociation,
    UnknownDissociation.name: UnknownDissociation
}


method = None
for method in dissociation_methods:
    dissociation_methods_map[method.name] = method
    dissociation_methods_map[method.id] = method
del method


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


__all__ = [
    "ActivationInformation", "MultipleActivationInformation", "DissociationMethod",
    "dissociation_methods", "UnknownDissociation", "dissociation_methods_map",
    "CID", "HCD", "ETD", "ECD", "supplemental_energy", "supplemental_term_map",
    "energy_terms"
]
