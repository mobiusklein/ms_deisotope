from __future__ import print_function


from .cv import Term, render_list


class Component(Term):
    """Describes a named component of a mass spectrometer, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`Component` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


def __generate_list_code():
    '''Prints the code to generate these static lists
    '''
    import sys
    render_list('ionization type', term_cls_name="Component")
    sys.stdout.write("\n\n")
    render_list('detector type', term_cls_name="Component")
    sys.stdout.write("\n\n")
    render_list('mass analyzer type', 'analyzer_types',
                term_cls_name="Component")
    sys.stdout.write("\n\n")
    render_list('inlet type', term_cls_name="Component")


ionization_types = [
    Component(u'multiphoton ionization', u'MS:1000227',
              u'Photoionization of an atom or molecule in which in two or more photons are absorbed.', 'ionization type', [u'ionization type']),
    Component(u'fast ion bombardment', u'MS:1000446',
              u'The ionization of any species by the interaction of a focused beam of ions having a translational energy of several thousand eV with a solid sample.', 'ionization type', [u'ionization type']),
    Component(u'resonance enhanced multiphoton ionization', u'MS:1000276',
              u'Multiphoton ionization in which the ionization cross section is significantly enhanced because the energy of the incident photons is resonant with an intermediate excited state of the neutral species.', 'ionization type', [u'ionization type']),
    Component(u'pyrolysis mass spectrometry', u'MS:1000274',
              u'A mass spectrometry technique in which the sample is heated to the point of decomposition and the gaseous decomposition products are introduced into the ion source.', 'ionization type', [u'ionization type']),
    Component(u'neutralization reionization mass spectrometry', u'MS:1000272',
              u'With this technique, m/z selected ions form neutrals by charge transfer to a collision gas or by dissociation. The neutrals are separated from the remaining ions and ionized in collisions with a second gas. This method is used to investigate reaction intermediates and other unstable species.', 'ionization type', [u'ionization type']),
    Component(u'photoionization', u'MS:1000273', u'The ionization of an atom or molecule by a photon, written M + h? ? M^+ + e. The term photon impact is not recommended.',
              'ionization type', [u'ionization type']),
    Component(u'Negative Ion chemical ionization', u'MS:1000271',
              u'Chemical ionization that results in the formation of negative ions.', 'ionization type', [u'ionization type']),
    Component(u'chemical ionization', u'MS:1000071',
              u'The formation of a new ion by the reaction of a neutral species with an ion. The process may involve transfer of an electron, a proton or other charged species between the reactants. When a positive ion results from chemical ionization the term may be used without qualification. When a negative ion results the term negative ion chemical ionization should be used. Note that this term is not synonymous with chemi-ionization.', 'ionization type', [u'ionization type']),
    Component(u'electrospray ionization', u'MS:1000073',
              u'A process in which ionized species in the gas phase are produced from an analyte-containing solution via highly charged fine droplets, by means of spraying the solution from a narrow-bore needle tip at atmospheric pressure in the presence of a high electric field. When a pressurized gas is used to aid in the formation of a stable spray, the term pneumatically assisted electrospray ionization is used. The term ion spray is not recommended.', 'ionization type', [u'ionization type']),
    Component(u'fast atom bombardment ionization', u'MS:1000074',
              u'The ionization of any species by the interaction of a focused beam of neutral atoms having a translational energy of several thousand eV with a sample that is typically dissolved in a solvent matrix. See also secondary ionization.', 'ionization type', [u'ionization type']),
    Component(u'flowing afterglow', u'MS:1000255',
              u'An ion source immersed in a flow of helium or other inert buffer gas that carries the ions through a meter-long reactor at pressures around 100 Pa.', 'ionization type', [u'ionization type']),
    Component(u'desorption ionization', u'MS:1000247',
              u'The formation of ions from a solid or liquid material after the rapid vaporization of that sample.', 'ionization type', [u'ionization type']),
    Component(u'atmospheric pressure ionization', u'MS:1000240',
              u'Any ionization process in which ions are formed in the gas phase at atmospheric pressure.', 'ionization type', [u'ionization type']),
    Component(u'spark ionization', u'MS:1000404',
              u'The formation of ions from a solid material by an intermittent electrical discharge.', 'ionization type', [u'ionization type']),
    Component(u'thermal ionization', u'MS:1000407',
              u'The ionization of a neutral species through contact with a high temperature surface.', 'ionization type', [u'ionization type']),
    Component(u'surface ionization', u'MS:1000406',
              u'The ionization of a neutral species when it interacts with a solid surface with an appropriate work function and temperature.', 'ionization type', [u'ionization type']),
    Component(u'plasma desorption ionization', u'MS:1000400',
              u'The ionization of material in a solid sample by bombarding it with ionic or neutral atoms formed as a result of the fission of a suitable nuclide, typically 252Cf. Synonymous with fission fragment ionization.', 'ionization type', [u'ionization type']),
    Component(u'soft ionization', u'MS:1000403', u'The formation of gas-phase ions without extensive fragmentation.',
              'ionization type', [u'ionization type']),
    Component(u'secondary ionization', u'MS:1000402',
              u'The process in which ions are ejected from a sample surface as a result of bombardment by a primary beam of atoms or ions.', 'ionization type', [u'ionization type']),
    Component(u'vertical ionization', u'MS:1000408',
              u'A process in which an electron is removed from or added to a molecule without a change in the positions of the atoms. The resulting ion is typically in an excited vibrational state.', 'ionization type', [u'ionization type']),
    Component(u'autodetachment', u'MS:1000383', u'The formation of a neutral when a negative ion in a discrete state with an energy greater than the detachment threshold loses an electron spontaneously without further interaction with an energy source.',
              'ionization type', [u'ionization type']),
    Component(u'adiabatic ionization', u'MS:1000380',
              u'A process whereby an electron is removed from an atom, ion, or molecule to produce an ion in its lowest energy state.', 'ionization type', [u'ionization type']),
    Component(u'associative ionization', u'MS:1000381',
              u'An ionization process in which two excited atoms or molecules react to form a single positive ion and an electron.', 'ionization type', [u'ionization type']),
    Component(u'chemi-ionization', u'MS:1000386',
              u'The reaction of a neutral molecule with an internally excited molecule to form an ion. Note that this term is not synonymous with chemical ionization.', 'ionization type', [u'ionization type']),
    Component(u'autoionization', u'MS:1000384', u'The formation of an ion when an atom or molecule in a discrete state with an energy greater than the ionization threshold loses an electron spontaneously without further interaction with an energy source.',
              'ionization type', [u'ionization type']),
    Component(u'charge exchange ionization', u'MS:1000385',
              u'The interaction of an ion with an atom or molecule in which the charge on the ion is transferred to the neutral without the dissociation of either. Synonymous with charge transfer ionization.', 'ionization type', [u'ionization type']),
    Component(u'dissociative ionization', u'MS:1000388',
              u'The reaction of a gas-phase molecule that results in its decomposition to form products, one of which is an ion.', 'ionization type', [u'ionization type']),
    Component(u'electron ionization', u'MS:1000389',
              u"The ionization of an atom or molecule by electrons that are typically accelerated to energies between 50 and 150 eV. Usually 70 eV electrons are used to produce positive ions. The term 'electron impact' is not recommended.", 'ionization type', [u'ionization type']),
    Component(u'field ionization', u'MS:1000258',
              u'The removal of electrons from any species by interaction with a high electric field.', 'ionization type', [u'ionization type']),
    Component(u'glow discharge ionization', u'MS:1000259',
              u'The formation of ions in the gas phase and from solid samples at the cathode by application of a voltage to a low pressure gas.', 'ionization type', [u'ionization type']),
    Component(u'liquid secondary ionization', u'MS:1000395',
              u'The ionization of any species by the interaction of a focused beam of ions with a sample that is dissolved in a solvent matrix. See also fast atom bombardment and secondary ionization.', 'ionization type', [u'ionization type']),
    Component(u'penning ionization', u'MS:1000399',
              u'Ionization that occurs through the interaction of two or more neutral gaseous species, at least one of which is internally excited.', 'ionization type', [u'ionization type']),
    Component(u'microelectrospray', u'MS:1000397', u'Electrospray ionization at a solvent flow rate of 300-800 nL/min where the flow is a result of a mechanical pump. See nanoelectrospray.',
              'ionization type', [u'electrospray ionization', u'ionization type']),
    Component(u'nanoelectrospray', u'MS:1000398', u'Electrospray ionization at a flow rate less than ~25 nL/min. Nanoelectrospray is synonymous with nanospray. The flow is dependent on the potenial on the tip of the electrospray needle and/or a gas presure to push the sample through the needle. See also electrospray ionization and microelectrospray.',
              'ionization type', [u'electrospray ionization', u'ionization type']),
    Component(u'matrix-assisted laser desorption ionization', u'MS:1000075',
              u'The formation of gas-phase ions from molecules that are present in a solid or solvent matrix that is irradiated with a pulsed laser. See also laser desorption/ionization.', 'ionization type', [u'desorption ionization', u'ionization type']),
    Component(u'field desorption', u'MS:1000257', u'The formation of gas-phase ions from a material deposited on a solid surface in the presence of a high electric field. Because this process may encompass ionization by field ionization or other mechanisms, it is not recommended as a synonym for field desorption ionization.',
              'ionization type', [u'desorption ionization', u'ionization type']),
    Component(u'surface-assisted laser desorption ionization', u'MS:1000405',
              u'The formation of gas-phase ions from molecules that are deposited on a particular surface substrate that is irradiated with a pulsed laser. See also matrix-assisted laser desorption ionization.', 'ionization type', [u'desorption ionization', u'ionization type']),
    Component(u'desorption/ionization on silicon', u'MS:1000387', u'The formation of ions by laser desorption ionization of a sample deposited on a porous silicon surface.',
              'ionization type', [u'desorption ionization', u'ionization type']),
    Component(u'laser desorption ionization', u'MS:1000393', u'The formation of gas-phase ions by the interaction of a pulsed laser with a solid or liquid material.',
              'ionization type', [u'desorption ionization', u'ionization type']),
    Component(u'atmospheric pressure matrix-assisted laser desorption ionization', u'MS:1000239',
              u'Matrix-assisted laser desorption ionization in which the sample target is at atmospheric pressure and the ions formed by the pulsed laser are sampled through a small aperture into the mass spectrometer.', 'ionization type', [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'atmospheric pressure chemical ionization', u'MS:1000070', u'Chemical ionization that takes place at atmospheric pressure as opposed to the reduced pressure is normally used for chemical ionization.',
              'ionization type', [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'desorption electrospray ionization', u'MS:1002011',
              u'Combination of electrospray and desorption ionization method that ionizes gases, liquids and solids in open air under atmospheric pressure." [DOI:10.1126/science.1104404', 'ionization type', [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'atmospheric pressure photoionization', u'MS:1000382', u'Atmospheric pressure chemical ionization in which the reactant ions are generated by photo-ionization.',
              'ionization type', [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'surface enhanced laser desorption ionization', u'MS:1000278',
              u'The formation of ionized species in the gas phase from analytes deposited on a particular surface substrate which is irradiated with a laser beam of which wavelength is absorbed by the surface. See also desorption/ionization on silicon and laser desorption/ionization.', 'ionization type', [u'surface ionization', u'ionization type']),
    Component(u'surface enhanced neat desorption', u'MS:1000279', u'Matrix-assisted laser desorption ionization in which the matrix is covalently linked to the target surface.',
              'ionization type', [u'surface ionization', u'ionization type']),
]


detector_types = [
    Component(u'channeltron', u'MS:1000107', u'A horn-shaped (or cone-shaped) continuous dynode particle multiplier. The ion strikes the inner surface of the device and induces the production of secondary electrons that in turn impinge on the inner surfaces to produce more secondary electrons. This avalanche effect produces an increase in signal in the final measured current pulse.',
              'detector type', [u'detector type']),
    Component(u'photomultiplier', u'MS:1000116',
              u'A detector for conversion of the ion/electron signal into photon(s) which are then amplified and detected.', 'detector type', [u'detector type']),
    Component(u'multi-collector', u'MS:1000115',
              u'A detector system commonly used in inductively coupled plasma mass spectrometers.', 'detector type', [u'detector type']),
    Component(u'faraday cup', u'MS:1000112', u'A conducting cup or chamber that intercepts a charged particle beam and is electrically connected to a current measuring device.',
              'detector type', [u'detector type']),
    Component(u'daly detector', u'MS:1000110',
              u'Detector consisting of a conversion dynode, scintillator and photomultiplier. The metal knob at high potential emits secondary electrons when ions impinge on the surface. The secondary electrons are accelerated onto the scintillator that produces light that is then detected by the photomultiplier detector.', 'detector type', [u'detector type']),
    Component(u'electron multiplier', u'MS:1000253',
              u'A device to amplify the current of a beam or packet of charged particles or photons by incidence upon the surface of an electrode to produce secondary electrons. The secondary electrons are then accelerated to other electrodes or parts of a continuous electrode to produce further secondary electrons.', 'detector type', [u'detector type']),
    Component(u'fluorescence detector', u'MS:1002308',
              u'A detector using a fluorescent signal after excitation with light.', 'detector type', [u'detector type']),
    Component(u'conversion dynode', u'MS:1000346',
              u'A surface that is held at high potential such that ions striking the surface produce electrons that are subsequently detected.', 'detector type', [u'detector type']),
    Component(u'dynode', u'MS:1000347', u'One of a series of electrodes in a photomultiplier tube. Such an arrangement is able to amplify the current emitted by the photocathode.',
              'detector type', [u'detector type']),
    Component(u'array detector', u'MS:1000345',
              u'Detector comprising several ion collection elements, arranged in a line or grid where each element is an individual detector.', 'detector type', [u'detector type']),
    Component(u'focal plane collector', u'MS:1000348',
              u'A detector for spatially disperse ion beams in which all ions simultaneously impinge on the detector plane.', 'detector type', [u'detector type']),
    Component(u'ion-to-photon detector', u'MS:1000349',
              u'A detector in which ions strike a conversion dynode to produce electrons that in turn strike a phosphor and the resulting photons are detected by a photomultiplier.', 'detector type', [u'detector type']),
    Component(u'postacceleration detector', u'MS:1000351',
              u'A detector in which the charged particles are accelerated to a high velocity and impinge on a conversion dynode, emitting secondary electrons. The electrons are accelerated onto a phosphor screen, which emits photons that are in turn detected using a photomultiplier or other photon detector.', 'detector type', [u'detector type']),
    Component(u'point collector', u'MS:1000350',
              u'A detector in which the ion beam is focused onto a point and the individual ions arrive sequentially.', 'detector type', [u'detector type']),
    Component(u'inductive detector', u'MS:1000624',
              u'Inductive detector.', 'detector type', [u'detector type']),
    Component(u'electron multiplier tube', u'MS:1000111', u'A device to amplify the current of a beam or packet of charged particles or photons by incidence upon the surface of an electrode to produce secondary electrons.',
              'detector type', [u'electron multiplier', u'detector type']),
    Component(u'Acquity UPLC FLR', u'MS:1000819', u'Acquity UPLC Fluorescence Detector.', 'detector type', [
              u'Waters instrument model', u'fluorescence detector', u'instrument model', u'detector type']),
    Component(u'conversion dynode electron multiplier', u'MS:1000108',
              u'A surface that is held at high potential so that ions striking the surface produce electrons that are subsequently detected.', 'detector type', [u'conversion dynode', u'detector type']),
    Component(u'conversion dynode photomultiplier', u'MS:1000109', u'A detector in which ions strike a conversion dynode to produce electrons that in turn generate photons through a phosphorescent screen that are detected by a photomultiplier.',
              'detector type', [u'conversion dynode', u'detector type']),
    Component(u'microchannel plate detector', u'MS:1000114', u'A thin plate that contains a closely spaced array of channels that each act as a continuous dynode particle multiplier. A charged particle, fast neutral particle, or photon striking the plate causes a cascade of secondary electrons that ultimately exits the opposite side of the plate.',
              'detector type', [u'array detector', u'detector type']),
    Component(u'photodiode array detector', u'MS:1000621', u'An array detector used to record spectra in the ultraviolet and visible region of light.',
              'detector type', [u'array detector', u'detector type']),
    Component(u'focal plane array', u'MS:1000113', u'An array of detectors for spatially disperse ion beams in which all ions simultaneously impinge on the detector plane.',
              'detector type', [u'focal plane collector', u'detector type']),
    Component(u'Acquity UPLC PDA', u'MS:1000818', u'Acquity UPLC Photodiode Array Detector.', 'detector type', [
              u'Waters instrument model', u'photodiode array detector', u'instrument model', u'array detector', u'detector type']),
]


analyzer_types = [
    Component(u'cyclotron', u'MS:1000288', u'A device that uses an oscillating electric field and magnetic field to accelerate charged particles.',
              'mass analyzer type', [u'mass analyzer type']),
    Component(u'orbitrap', u'MS:1000484', u'An ion trapping device that consists of an outer barrel-like electrode and a coaxial inner spindle-like electrode that form an electrostatic field with quadro-logarithmic potential distribution. The frequency of harmonic oscillations of the orbitally trapped ions along the axis of the electrostatic field is independent of the ion velocity and is inversely proportional to the square root of m/z so that the trap can be used as a mass analyzer.',
              'mass analyzer type', [u'mass analyzer type']),
    Component(u'ion trap', u'MS:1000264', u'A device for spatially confining ions using electric and magnetic fields alone or in combination.',
              'mass analyzer type', [u'mass analyzer type']),
    Component(u'fourier transform ion cyclotron resonance mass spectrometer', u'MS:1000079',
              u'A mass spectrometer based on the principle of ion cyclotron resonance in which an ion in a magnetic field moves in a circular orbit at a frequency characteristic of its m/z value. Ions are coherently excited to a larger radius orbit using a pulse of radio frequency energy and their image charge is detected on receiver plates as a time domain signal. Fourier transformation of the time domain signal results in a frequency domain signal which is converted to a mass spectrum based in the inverse relationship between frequency and m/z.', 'mass analyzer type', [u'mass analyzer type']),
    Component(u'electrostatic energy analyzer', u'MS:1000254',
              u'A device consisting of conducting parallel plates, concentric cylinders or concentric spheres that separates charged particles according to their kinetic energy by means of an electric field that is constant in time.', 'mass analyzer type', [u'mass analyzer type']),
    Component(u'quadrupole', u'MS:1000081', u'A mass spectrometer that consists of four parallel rods whose centers form the corners of a square and whose opposing poles are connected. The voltage applied to the rods is a superposition of a static potential and a sinusoidal radio frequency potential. The motion of an ion in the x and y dimensions is described by the Matthieu equation whose solutions show that ions in a particular m/z range can be transmitted along the z axis.',
              'mass analyzer type', [u'mass analyzer type']),
    Component(u'magnetic sector', u'MS:1000080', u'A device that produces a magnetic field perpendicular to a charged particle beam that deflects the beam to an extent that is proportional to the particle momentum per unit charge. For a monoenergetic beam, the deflection is proportional to m/z.',
              'mass analyzer type', [u'mass analyzer type']),
    Component(u'time-of-flight', u'MS:1000084', u'Instrument that separates ions by m/z in a field-free region after acceleration to a fixed acceleration energy.',
              'mass analyzer type', [u'mass analyzer type']),
    Component(u'stored waveform inverse fourier transform', u'MS:1000284',
              u'A technique to create excitation waveforms for ions in FT-ICR mass spectrometer or Paul ion trap. An excitation waveform in the time-domain is generated by taking the inverse Fourier transform of an appropriate frequency-domain programmed excitation spectrum, in which the resonance frequencies of ions to be excited are included. This technique may be used for selection of precursor ions in MS2 experiments.', 'mass analyzer type', [u'mass analyzer type']),
    Component(u'quadrupole ion trap', u'MS:1000082', u'Quadrupole Ion Trap mass analyzer captures the ions in a three dimensional ion trap and then selectively ejects them by varying the RF and DC potentials.',
              'mass analyzer type', [u'ion trap', u'mass analyzer type']),
    Component(u'linear ion trap', u'MS:1000291', u'A two dimensional Paul ion trap in which ions are confined in the axial dimension by means of an electric field at the ends of the trap.',
              'mass analyzer type', [u'ion trap', u'mass analyzer type']),
    Component(u'axial ejection linear ion trap', u'MS:1000078', u'A linear ion trap mass spectrometer where ions are ejected along the axis of the analyzer.',
              'mass analyzer type', [u'linear ion trap', u'ion trap', u'mass analyzer type']),
    Component(u'radial ejection linear ion trap', u'MS:1000083', u'A linear ion trap mass spectrometer where ions are ejected along the radius of the analyzer.',
              'mass analyzer type', [u'linear ion trap', u'ion trap', u'mass analyzer type']),
]


inlet_types = [
    Component(u'flow injection analysis', u'MS:1000058',
              u'Sample is directly injected or infused into the ionization source.', 'inlet type', [u'inlet type']),
    Component(u'inductively coupled plasma', u'MS:1000059',
              u'A gas discharge ion source in which the energy to the plasma is supplied by electromagnetic induction.', 'inlet type', [u'inlet type']),
    Component(u'direct inlet', u'MS:1000056',
              u'The sample is directly inserted into the ion source, usually on the end of a heatable probe.', 'inlet type', [u'inlet type']),
    Component(u'electrospray inlet', u'MS:1000057',
              u'Inlet used for introducing the liquid sample into an electrospray ionization source.', 'inlet type', [u'inlet type']),
    Component(u'continuous flow fast atom bombardment', u'MS:1000055',
              u'Fast atom bombardment ionization in which the analyte in solution is entrained in a flowing liquid matrix.', 'inlet type', [u'inlet type']),
    Component(u'reservoir', u'MS:1000067',
              u'A sample inlet method involving a reservoir.', 'inlet type', [u'inlet type']),
    Component(u'particle beam', u'MS:1000066',
              u'Method for generating ions from a solution of an analyte.', 'inlet type', [u'inlet type']),
    Component(u'open split', u'MS:1000065',
              u'A division of flowing stream of liquid into two streams.', 'inlet type', [u'inlet type']),
    Component(u'moving wire', u'MS:1000064',
              u'Continuous moving surface in the form of a wire which passes through an ion source carrying analyte molecules.', 'inlet type', [u'inlet type']),
    Component(u'moving belt', u'MS:1000063',
              u'Continuous moving surface in the form of a belt which passes through an ion source carrying analyte molecules.', 'inlet type', [u'inlet type']),
    Component(u'membrane separator', u'MS:1000062',
              u'A device to separate carrier molecules from analyte molecules on the basis of ease of diffusion across a semipermeable membrane.', 'inlet type', [u'inlet type']),
    Component(u'jet separator', u'MS:1000061',
              u'A device that separates carrier gas from gaseous analyte molecules on the basis of diffusivity.', 'inlet type', [u'inlet type']),
    Component(u'infusion', u'MS:1000060',
              u'The continuous flow of solution of a sample into the ionization source.', 'inlet type', [u'inlet type']),
    Component(u'thermospray inlet', u'MS:1000069',
              u'A method for generating gas phase ions from a solution of an analyte by rapid heating of the sample.', 'inlet type', [u'inlet type']),
    Component(u'septum', u'MS:1000068',
              u'A disc composed of a flexible material that seals the entrance to the reservoir. Can also be entrance to the vacuum chamber.', 'inlet type', [u'inlet type']),
    Component(u'direct liquid introduction', u'MS:1000249',
              u'The delivery of a liquid sample into a mass spectrometer for spray or desorption ionization.', 'inlet type', [u'inlet type']),
    Component(u'direct insertion probe', u'MS:1000248',
              u'A device for introducing a solid or liquid sample into a mass spectrometer ion source for desorption ionization.', 'inlet type', [u'inlet type']),
    Component(u'membrane inlet', u'MS:1000396',
              u'A semi-permeable membrane separator that permits the passage of gas sample directly to the mass spectrometer ion source.', 'inlet type', [u'inlet type']),
    Component(u'nanospray inlet', u'MS:1000485', u'Nanospray Inlet.',
              'inlet type', [u'electrospray inlet', u'inlet type']),
]


all_components = ionization_types + detector_types + analyzer_types + inlet_types

all_components_by_name = {c.name: c for c in all_components}


def component(name):
    try:
        return all_components_by_name[name]
    except KeyError:
        return Component(name, name, name, name, [name])


class ComponentGroup(object):

    def __init__(self, type, members, order):
        self.type = type
        self.members = list(members)
        self.order = int(order)

    def __repr__(self):
        t = "{s.__class__.__name__}({s.type!r}, {s.members}, order={s.order})"
        return t.format(s=self)

    def __getitem__(self, i):
        return self.members[i]

    def __setitem__(self, i, v):
        self.members[i] = v

    def add(self, v):
        self.members.append(v)

    def __len__(self):
        return len(self.members)


class InstrumentInformation(object):

    def __init__(self, id, groups):
        self.id = id
        self.groups = sorted(groups, key=lambda x: x.order)
        self.analyzers = []

        for group in self.groups:
            if group.type == 'analyzer':
                self.analyzers.extend(group)

    def __getitem__(self, i):
        return self.groups[i]

    def __len__(self):
        return len(self.groups)

    def __repr__(self):
        return "{self.__class__.__name__}({self.id!r}, {self.groups})".format(
            self=self)


if __name__ == '__main__':
    __generate_list_code()
