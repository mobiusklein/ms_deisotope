from .cv import Term, render_list


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
        return (self.method == other.method) and abs(self.energy - other.energy) < 1e-3

    def __ne__(self, other):
        return not (self == other)


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


def __generate_list_code():
    '''Prints the code to generate these static lists
    '''
    render_list('dissociation method', term_cls_name="DissociationMethod")


dissociation_methods = []


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('dissociation method', term_cls_name="DissociationMethod", writer=cog.out)
# ]]]
dissociation_methods = [
    DissociationMethod(u'sustained off-resonance irradiation', u'MS:1000282',
                       u"A technique associated with Fourier transform ion cyclotron resonance (FT-ICR) mass spectrometry to carry out ion/neutral reactions such as low-energy collision-induced dissociation. A radio-frequency electric field of slightly off-resonance to the cyclotron frequency of the reactant ion cyclically accelerates and decelerates the reactant ion that is confined in the Penning ion trap. The ion's orbit does not exceed the dimensions of ion trap while the ion undergoes an ion/neutral species process that produces a high average translational energy for an extended time.",       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'post-source decay', u'MS:1000135', u'A technique specific to reflectron time-of-flight mass spectrometers where product ions of metastable transitions or collision-induced dissociations generated in the drift tube prior to entering the reflectron are m/z separated to yield product ion spectra.',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'plasma desorption', u'MS:1000134', u'The ionization of material in a solid sample by bombarding it with ionic or neutral atoms formed as a result of the fission of a suitable nuclide, typically 252Cf. Synonymous with fission fragment ionization.',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'surface-induced dissociation', u'MS:1000136',
                       u'Fragmentation that results from the collision of an ion with a surface.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'collision-induced dissociation', u'MS:1000133',
                       u'The dissociation of an ion after collisional excitation. The term collisional-activated dissociation is not recommended.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'pulsed q dissociation', u'MS:1000599', u'A process that involves precursor ion activation at high Q, a time delay to allow the precursor to fragment, then a rapid pulse to low Q where all fragment ions are trapped. The product ions can then be scanned out of the ion trap and detected.',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'electron transfer dissociation', u'MS:1000598',
                       u'A process to fragment ions in a mass spectrometer by inducing fragmentation of cations (e.g. peptides or proteins) by transferring electrons to them.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'in-source collision-induced dissociation', u'MS:1001880',
                       u'The dissociation of an ion as a result of collisional excitation during ion transfer from an atmospheric pressure ion source and the mass spectrometer vacuum.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'infrared multiphoton dissociation', u'MS:1000262',
                       u'Multiphoton ionization where the reactant ion dissociates as a result of the absorption of multiple infrared photons.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'blackbody infrared radiative dissociation', u'MS:1000242',
                       u'A special case of infrared multiphoton dissociation wherein excitation of the reactant ion is caused by absorption of infrared photons radiating from heated blackbody surroundings, which are usually the walls of a vacuum chamber. See also infrared multiphoton dissociation.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'low-energy collision-induced dissociation', u'MS:1000433',
                       u'A collision-induced dissociation process wherein the precursor ion has the translational energy lower than approximately 1000 eV. This process typically requires multiple collisions and the collisional excitation is cumulative.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'photodissociation', u'MS:1000435', u'A process wherein the reactant ion is dissociated as a result of absorption of one or more photons.',
                       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'LIFT', u'MS:1002000',
                       u'A Bruker\'s proprietary technique where molecular ions are initially accelerated at lower energy, then collide with inert gas in a collision cell that is then \'lifted\' to high potential. The use of inert gas is optional, as it could lift also fragments provided by LID." [DOI:10.1007/s00216-003-2057-0 ',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)', u'MS:1002631',
                       u'A dissociation process combining electron-transfer and higher-energy collision dissociation (EThcD). It combines ETD (reaction time) followed by HCD (activation energy).',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'electron capture dissociation', u'MS:1000250',
                       u'A process in which a multiply protonated molecules interacts with a low energy electrons. Capture of the electron leads the liberation of energy and a reduction in charge state of the ion with the production of the (M + nH) (n-1)+ odd electron ion, which readily fragments.',       'dissociation method', [u'dissociation method']),
    DissociationMethod(u'trap-type collision-induced dissociation', u'MS:1002472', u'A collision-induced dissociation process that occurs in a trap-type collision cell.',
                       'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'beam-type collision-induced dissociation', u'MS:1000422', u'A collision-induced dissociation process that occurs in a beam-type collision cell.',
                       'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'supplemental collision-induced dissociation', u'MS:1002679', u'The dissociation of an ion after supplemental collisional excitation.',
                       'dissociation method', [u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'higher energy beam-type collision-induced dissociation', u'MS:1002481', u'A collision-induced dissociation process wherein the projectile ion has the translational energy higher than approximately 1000 eV.',
                       'dissociation method', [u'beam-type collision-induced dissociation', u'collision-induced dissociation', u'dissociation method']),
    DissociationMethod(u'supplemental beam-type collision-induced dissociation', u'MS:1002678', u'A supplemental collision-induced dissociation process that occurs in a beam-type collision cell in addition to another primary type of dissociation.',
                       'dissociation method', [u'beam-type collision-induced dissociation', u'collision-induced dissociation', u'dissociation method']),
]
# [[[end]]]


UnknownDissociation = DissociationMethod(
    "unknown dissociation", None, None, 'dissociation method',
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


if __name__ == '__main__':
    __generate_list_code()
