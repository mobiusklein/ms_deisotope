from __future__ import print_function

from ms_deisotope.utils import uid


from .cv import Term, TermSet, render_list


class Component(Term):
    """Describes a named component of a mass spectrometer, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`Component` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('ionization type', term_cls_name="Component", writer=cog.out)
# ]]]
# CV Version: 4.1.55
ionization_types = TermSet([
    Component('chemical ionization', 'MS:1000071',
              ('The formation of a new ion by the reaction of a neutral '
               'species with an ion. The process may involve transfer of an '
               'electron, a proton or other charged species between the '
               'reactants. When a positive ion results from chemical '
               'ionization the term may be used without qualification. When '
               'a negative ion results the term negative ion chemical '
               'ionization should be used. Note that this term is not '
               'synonymous with chemi-ionization.'),
              'ionization type',
              ['ionization type']),
    Component('electrospray ionization', 'MS:1000073',
              ('A process in which ionized species in the gas phase are '
               'produced from an analyte-containing solution via highly '
               'charged fine droplets, by means of spraying the solution '
               'from a narrow-bore needle tip at atmospheric pressure in the '
               'presence of a high electric field. When a pressurized gas is '
               'used to aid in the formation of a stable spray, the term '
               'pneumatically assisted electrospray ionization is used. The '
               'term ion spray is not recommended.'),
              'ionization type',
              ['ionization type']),
    Component('fast atom bombardment ionization', 'MS:1000074',
              ('The ionization of any species by the interaction of a '
               'focused beam of neutral atoms having a translational energy '
               'of several thousand eV with a sample that is typically '
               'dissolved in a solvent matrix. See also secondary '
               'ionization.'),
              'ionization type',
              ['ionization type']),
    Component('multiphoton ionization', 'MS:1000227',
              ('Photoionization of an atom or molecule in which in two or '
               'more photons are absorbed.'),
              'ionization type',
              ['ionization type']),
    Component('atmospheric pressure ionization', 'MS:1000240',
              ('Any ionization process in which ions are formed in the gas '
               'phase at atmospheric pressure.'),
              'ionization type',
              ['ionization type']),
    Component('desorption ionization', 'MS:1000247',
              ('The formation of ions from a solid or liquid material after '
               'the rapid vaporization of that sample.'),
              'ionization type',
              ['ionization type']),
    Component('flowing afterglow', 'MS:1000255',
              ('An ion source immersed in a flow of helium or other inert '
               'buffer gas that carries the ions through a meter-long '
               'reactor at pressures around 100 Pa.'),
              'ionization type',
              ['ionization type']),
    Component('field ionization', 'MS:1000258',
              ('The removal of electrons from any species by interaction '
               'with a high electric field.'),
              'ionization type',
              ['ionization type']),
    Component('glow discharge ionization', 'MS:1000259',
              ('The formation of ions in the gas phase and from solid '
               'samples at the cathode by application of a voltage to a low '
               'pressure gas.'),
              'ionization type',
              ['ionization type']),
    Component('Negative Ion chemical ionization', 'MS:1000271',
              ('Chemical ionization that results in the formation of '
               'negative ions.'),
              'ionization type',
              ['ionization type']),
    Component('neutralization reionization mass spectrometry', 'MS:1000272',
              ('With this technique, m/z selected ions form neutrals by '
               'charge transfer to a collision gas or by dissociation. The '
               'neutrals are separated from the remaining ions and ionized '
               'in collisions with a second gas. This method is used to '
               'investigate reaction intermediates and other unstable '
               'species.'),
              'ionization type',
              ['ionization type']),
    Component('photoionization', 'MS:1000273',
              ('The ionization of an atom or molecule by a photon, written M '
               '+ h? ? M^+ + e. The term photon impact is not recommended.'),
              'ionization type',
              ['ionization type']),
    Component('pyrolysis mass spectrometry', 'MS:1000274',
              ('A mass spectrometry technique in which the sample is heated '
               'to the point of decomposition and the gaseous decomposition '
               'products are introduced into the ion source.'),
              'ionization type',
              ['ionization type']),
    Component('resonance enhanced multiphoton ionization', 'MS:1000276',
              ('Multiphoton ionization in which the ionization cross section '
               'is significantly enhanced because the energy of the incident '
               'photons is resonant with an intermediate excited state of '
               'the neutral species.'),
              'ionization type',
              ['ionization type']),
    Component('adiabatic ionization', 'MS:1000380',
              ('A process whereby an electron is removed from an atom, ion, '
               'or molecule to produce an ion in its lowest energy state.'),
              'ionization type',
              ['ionization type']),
    Component('associative ionization', 'MS:1000381',
              ('An ionization process in which two excited atoms or '
               'molecules react to form a single positive ion and an '
               'electron.'),
              'ionization type',
              ['ionization type']),
    Component('autodetachment', 'MS:1000383',
              ('The formation of a neutral when a negative ion in a discrete '
               'state with an energy greater than the detachment threshold '
               'loses an electron spontaneously without further interaction '
               'with an energy source.'),
              'ionization type',
              ['ionization type']),
    Component('autoionization', 'MS:1000384',
              ('The formation of an ion when an atom or molecule in a '
               'discrete state with an energy greater than the ionization '
               'threshold loses an electron spontaneously without further '
               'interaction with an energy source.'),
              'ionization type',
              ['ionization type']),
    Component('charge exchange ionization', 'MS:1000385',
              ('The interaction of an ion with an atom or molecule in which '
               'the charge on the ion is transferred to the neutral without '
               'the dissociation of either. Synonymous with charge transfer '
               'ionization.'),
              'ionization type',
              ['ionization type']),
    Component('chemi-ionization', 'MS:1000386',
              ('The reaction of a neutral molecule with an internally '
               'excited molecule to form an ion. Note that this term is not '
               'synonymous with chemical ionization.'),
              'ionization type',
              ['ionization type']),
    Component('dissociative ionization', 'MS:1000388',
              ('The reaction of a gas-phase molecule that results in its '
               'decomposition to form products, one of which is an ion.'),
              'ionization type',
              ['ionization type']),
    Component('electron ionization', 'MS:1000389',
              ('The ionization of an atom or molecule by electrons that are '
               'typically accelerated to energies between 50 and 150 eV. '
               'Usually 70 eV electrons are used to produce positive ions. '
               "The term 'electron impact' is not recommended."),
              'ionization type',
              ['ionization type']),
    Component('liquid secondary ionization', 'MS:1000395',
              ('The ionization of any species by the interaction of a '
               'focused beam of ions with a sample that is dissolved in a '
               'solvent matrix. See also fast atom bombardment and secondary '
               'ionization.'),
              'ionization type',
              ['ionization type']),
    Component('penning ionization', 'MS:1000399',
              ('Ionization that occurs through the interaction of two or '
               'more neutral gaseous species, at least one of which is '
               'internally excited.'),
              'ionization type',
              ['ionization type']),
    Component('plasma desorption ionization', 'MS:1000400',
              ('The ionization of material in a solid sample by bombarding '
               'it with ionic or neutral atoms formed as a result of the '
               'fission of a suitable nuclide, typically 252Cf. Synonymous '
               'with fission fragment ionization.'),
              'ionization type',
              ['ionization type']),
    Component('secondary ionization', 'MS:1000402',
              ('The process in which ions are ejected from a sample surface '
               'as a result of bombardment by a primary beam of atoms or '
               'ions.'),
              'ionization type',
              ['ionization type']),
    Component('soft ionization', 'MS:1000403',
              ('The formation of gas-phase ions without extensive '
               'fragmentation.'),
              'ionization type',
              ['ionization type']),
    Component('spark ionization', 'MS:1000404',
              ('The formation of ions from a solid material by an '
               'intermittent electrical discharge.'),
              'ionization type',
              ['ionization type']),
    Component('surface ionization', 'MS:1000406',
              ('The ionization of a neutral species when it interacts with a '
               'solid surface with an appropriate work function and '
               'temperature.'),
              'ionization type',
              ['ionization type']),
    Component('thermal ionization', 'MS:1000407',
              ('The ionization of a neutral species through contact with a '
               'high temperature surface.'),
              'ionization type',
              ['ionization type']),
    Component('vertical ionization', 'MS:1000408',
              ('A process in which an electron is removed from or added to a '
               'molecule without a change in the positions of the atoms. The '
               'resulting ion is typically in an excited vibrational state.'),
              'ionization type',
              ['ionization type']),
    Component('fast ion bombardment', 'MS:1000446',
              ('The ionization of any species by the interaction of a '
               'focused beam of ions having a translational energy of '
               'several thousand eV with a solid sample.'),
              'ionization type',
              ['ionization type']),
    Component('microelectrospray', 'MS:1000397',
              ('Electrospray ionization at a solvent flow rate of 300-800 '
               'nL/min where the flow is a result of a mechanical pump. See '
               'nanoelectrospray.'),
              'ionization type',
              ['electrospray ionization', 'ionization type']),
    Component('nanoelectrospray', 'MS:1000398',
              ('Electrospray ionization at a flow rate less than ~25 nL/min. '
               'Nanoelectrospray is synonymous with nanospray. The flow is '
               'dependent on the potenial on the tip of the electrospray '
               'needle and/or a gas presure to push the sample through the '
               'needle. See also electrospray ionization and '
               'microelectrospray.'),
              'ionization type',
              ['electrospray ionization', 'ionization type']),
    Component('atmospheric pressure chemical ionization', 'MS:1000070',
              ('Chemical ionization that takes place at atmospheric pressure '
               'as opposed to the reduced pressure is normally used for '
               'chemical ionization.'),
              'ionization type',
              ['atmospheric pressure ionization', 'ionization type']),
    Component('atmospheric pressure matrix-assisted laser desorption ionization', 'MS:1000239',
              ('Matrix-assisted laser desorption ionization in which the '
               'sample target is at atmospheric pressure and the ions formed '
               'by the pulsed laser are sampled through a small aperture '
               'into the mass spectrometer.'),
              'ionization type',
              ['atmospheric pressure ionization', 'ionization type']),
    Component('atmospheric pressure photoionization', 'MS:1000382',
              ('Atmospheric pressure chemical ionization in which the '
               'reactant ions are generated by photo-ionization.'),
              'ionization type',
              ['atmospheric pressure ionization', 'ionization type']),
    Component('desorption electrospray ionization', 'MS:1002011',
              ('Combination of electrospray and desorption ionization method '
               'that ionizes gases, liquids and solids in open air under '
               'atmospheric pressure." [DOI:10.1126/science.1104404'),
              'ionization type',
              ['atmospheric pressure ionization', 'ionization type']),
    Component('matrix-assisted laser desorption ionization', 'MS:1000075',
              ('The formation of gas-phase ions from molecules that are '
               'present in a solid or solvent matrix that is irradiated with '
               'a pulsed laser. See also laser desorption/ionization.'),
              'ionization type',
              ['desorption ionization', 'ionization type']),
    Component('field desorption', 'MS:1000257',
              ('The formation of gas-phase ions from a material deposited on '
               'a solid surface in the presence of a high electric field. '
               'Because this process may encompass ionization by field '
               'ionization or other mechanisms, it is not recommended as a '
               'synonym for field desorption ionization.'),
              'ionization type',
              ['desorption ionization', 'ionization type']),
    Component('desorption/ionization on silicon', 'MS:1000387',
              ('The formation of ions by laser desorption ionization of a '
               'sample deposited on a porous silicon surface.'),
              'ionization type',
              ['desorption ionization', 'ionization type']),
    Component('laser desorption ionization', 'MS:1000393',
              ('The formation of gas-phase ions by the interaction of a '
               'pulsed laser with a solid or liquid material.'),
              'ionization type',
              ['desorption ionization', 'ionization type']),
    Component('surface-assisted laser desorption ionization', 'MS:1000405',
              ('The formation of gas-phase ions from molecules that are '
               'deposited on a particular surface substrate that is '
               'irradiated with a pulsed laser. See also matrix-assisted '
               'laser desorption ionization.'),
              'ionization type',
              ['desorption ionization', 'ionization type']),
    Component('surface enhanced laser desorption ionization', 'MS:1000278',
              ('The formation of ionized species in the gas phase from '
               'analytes deposited on a particular surface substrate which '
               'is irradiated with a laser beam of which wavelength is '
               'absorbed by the surface. See also desorption/ionization on '
               'silicon and laser desorption/ionization.'),
              'ionization type',
              ['surface ionization', 'ionization type']),
    Component('surface enhanced neat desorption', 'MS:1000279',
              ('Matrix-assisted laser desorption ionization in which the '
               'matrix is covalently linked to the target surface.'),
              'ionization type',
              ['surface ionization', 'ionization type']),
])
# [[[end]]]


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('detector type', term_cls_name="Component", writer=cog.out)
# ]]]
# CV Version: 4.1.55
detector_types = TermSet([
    Component('channeltron', 'MS:1000107',
              ('A horn-shaped (or cone-shaped) continuous dynode particle '
               'multiplier. The ion strikes the inner surface of the device '
               'and induces the production of secondary electrons that in '
               'turn impinge on the inner surfaces to produce more secondary '
               'electrons. This avalanche effect produces an increase in '
               'signal in the final measured current pulse.'),
              'detector type',
              ['detector type']),
    Component('daly detector', 'MS:1000110',
              ('Detector consisting of a conversion dynode, scintillator and '
               'photomultiplier. The metal knob at high potential emits '
               'secondary electrons when ions impinge on the surface. The '
               'secondary electrons are accelerated onto the scintillator '
               'that produces light that is then detected by the '
               'photomultiplier detector.'),
              'detector type',
              ['detector type']),
    Component('faraday cup', 'MS:1000112',
              ('A conducting cup or chamber that intercepts a charged '
               'particle beam and is electrically connected to a current '
               'measuring device.'),
              'detector type',
              ['detector type']),
    Component('multi-collector', 'MS:1000115',
              ('A detector system commonly used in inductively coupled '
               'plasma mass spectrometers.'),
              'detector type',
              ['detector type']),
    Component('photomultiplier', 'MS:1000116',
              ('A detector for conversion of the ion/electron signal into '
               'photon(s) which are then amplified and detected.'),
              'detector type',
              ['detector type']),
    Component('electron multiplier', 'MS:1000253',
              ('A device to amplify the current of a beam or packet of '
               'charged particles or photons by incidence upon the surface '
               'of an electrode to produce secondary electrons. The '
               'secondary electrons are then accelerated to other electrodes '
               'or parts of a continuous electrode to produce further '
               'secondary electrons.'),
              'detector type',
              ['detector type']),
    Component('array detector', 'MS:1000345',
              ('Detector comprising several ion collection elements, '
               'arranged in a line or grid where each element is an '
               'individual detector.'),
              'detector type',
              ['detector type']),
    Component('conversion dynode', 'MS:1000346',
              ('A surface that is held at high potential such that ions '
               'striking the surface produce electrons that are subsequently '
               'detected.'),
              'detector type',
              ['detector type']),
    Component('dynode', 'MS:1000347',
              ('One of a series of electrodes in a photomultiplier tube. '
               'Such an arrangement is able to amplify the current emitted '
               'by the photocathode.'),
              'detector type',
              ['detector type']),
    Component('focal plane collector', 'MS:1000348',
              ('A detector for spatially disperse ion beams in which all '
               'ions simultaneously impinge on the detector plane.'),
              'detector type',
              ['detector type']),
    Component('ion-to-photon detector', 'MS:1000349',
              ('A detector in which ions strike a conversion dynode to '
               'produce electrons that in turn strike a phosphor and the '
               'resulting photons are detected by a photomultiplier.'),
              'detector type',
              ['detector type']),
    Component('point collector', 'MS:1000350',
              ('A detector in which the ion beam is focused onto a point and '
               'the individual ions arrive sequentially.'),
              'detector type',
              ['detector type']),
    Component('postacceleration detector', 'MS:1000351',
              ('A detector in which the charged particles are accelerated to '
               'a high velocity and impinge on a conversion dynode, emitting '
               'secondary electrons. The electrons are accelerated onto a '
               'phosphor screen, which emits photons that are in turn '
               'detected using a photomultiplier or other photon detector.'),
              'detector type',
              ['detector type']),
    Component('inductive detector', 'MS:1000624',
              ('Inductive detector.'),
              'detector type',
              ['detector type']),
    Component('fluorescence detector', 'MS:1002308',
              ('A detector using a fluorescent signal after excitation with '
               'light.'),
              'detector type',
              ['detector type']),
    Component('electron multiplier tube', 'MS:1000111',
              ('A device to amplify the current of a beam or packet of '
               'charged particles or photons by incidence upon the surface '
               'of an electrode to produce secondary electrons.'),
              'detector type',
              ['electron multiplier', 'detector type']),
    Component('microchannel plate detector', 'MS:1000114',
              ('A thin plate that contains a closely spaced array of '
               'channels that each act as a continuous dynode particle '
               'multiplier. A charged particle, fast neutral particle, or '
               'photon striking the plate causes a cascade of secondary '
               'electrons that ultimately exits the opposite side of the '
               'plate.'),
              'detector type',
              ['array detector', 'detector type']),
    Component('photodiode array detector', 'MS:1000621',
              ('An array detector used to record spectra in the ultraviolet '
               'and visible region of light.'),
              'detector type',
              ['array detector', 'detector type']),
    Component('conversion dynode electron multiplier', 'MS:1000108',
              ('A surface that is held at high potential so that ions '
               'striking the surface produce electrons that are subsequently '
               'detected.'),
              'detector type',
              ['conversion dynode', 'detector type']),
    Component('conversion dynode photomultiplier', 'MS:1000109',
              ('A detector in which ions strike a conversion dynode to '
               'produce electrons that in turn generate photons through a '
               'phosphorescent screen that are detected by a '
               'photomultiplier.'),
              'detector type',
              ['conversion dynode', 'detector type']),
    Component('focal plane array', 'MS:1000113',
              ('An array of detectors for spatially disperse ion beams in '
               'which all ions simultaneously impinge on the detector plane.'),
              'detector type',
              ['focal plane collector', 'detector type']),
    Component('Acquity UPLC FLR', 'MS:1000819',
              ('Acquity UPLC Fluorescence Detector.'),
              'detector type',
              ['Waters instrument model', 'fluorescence detector', 'instrument model', 'detector type']),
    Component('Acquity UPLC PDA', 'MS:1000818',
              ('Acquity UPLC Photodiode Array Detector.'),
              'detector type',
              ['Waters instrument model', 'photodiode array detector', 'instrument model', 'array detector', 'detector type']),
])
# [[[end]]]


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('mass analyzer type', 'analyzer_types', term_cls_name="Component", writer=cog.out)
# ]]]
# CV Version: 4.1.55
analyzer_types = TermSet([
    Component('fourier transform ion cyclotron resonance mass spectrometer', 'MS:1000079',
              ('A mass spectrometer based on the principle of ion cyclotron '
               'resonance in which an ion in a magnetic field moves in a '
               'circular orbit at a frequency characteristic of its m/z '
               'value. Ions are coherently excited to a larger radius orbit '
               'using a pulse of radio frequency energy and their image '
               'charge is detected on receiver plates as a time domain '
               'signal. Fourier transformation of the time domain signal '
               'results in a frequency domain signal which is converted to a '
               'mass spectrum based in the inverse relationship between '
               'frequency and m/z.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('magnetic sector', 'MS:1000080',
              ('A device that produces a magnetic field perpendicular to a '
               'charged particle beam that deflects the beam to an extent '
               'that is proportional to the particle momentum per unit '
               'charge. For a monoenergetic beam, the deflection is '
               'proportional to m/z.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('quadrupole', 'MS:1000081',
              ('A mass spectrometer that consists of four parallel rods '
               'whose centers form the corners of a square and whose '
               'opposing poles are connected. The voltage applied to the '
               'rods is a superposition of a static potential and a '
               'sinusoidal radio frequency potential. The motion of an ion '
               'in the x and y dimensions is described by the Matthieu '
               'equation whose solutions show that ions in a particular m/z '
               'range can be transmitted along the z axis.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('time-of-flight', 'MS:1000084',
              ('Instrument that separates ions by m/z in a field-free region '
               'after acceleration to a fixed acceleration energy.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('electrostatic energy analyzer', 'MS:1000254',
              ('A device consisting of conducting parallel plates, '
               'concentric cylinders or concentric spheres that separates '
               'charged particles according to their kinetic energy by means '
               'of an electric field that is constant in time.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('ion trap', 'MS:1000264',
              ('A device for spatially confining ions using electric and '
               'magnetic fields alone or in combination.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('stored waveform inverse fourier transform', 'MS:1000284',
              ('A technique to create excitation waveforms for ions in FT- '
               'ICR mass spectrometer or Paul ion trap. An excitation '
               'waveform in the time-domain is generated by taking the '
               'inverse Fourier transform of an appropriate frequency-domain '
               'programmed excitation spectrum, in which the resonance '
               'frequencies of ions to be excited are included. This '
               'technique may be used for selection of precursor ions in MS2 '
               'experiments.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('cyclotron', 'MS:1000288',
              ('A device that uses an oscillating electric field and '
               'magnetic field to accelerate charged particles.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('orbitrap', 'MS:1000484',
              ('An ion trapping device that consists of an outer barrel-like '
               'electrode and a coaxial inner spindle-like electrode that '
               'form an electrostatic field with quadro-logarithmic '
               'potential distribution. The frequency of harmonic '
               'oscillations of the orbitally trapped ions along the axis of '
               'the electrostatic field is independent of the ion velocity '
               'and is inversely proportional to the square root of m/z so '
               'that the trap can be used as a mass analyzer.'),
              'mass analyzer type',
              ['mass analyzer type']),
    Component('quadrupole ion trap', 'MS:1000082',
              ('Quadrupole Ion Trap mass analyzer captures the ions in a '
               'three dimensional ion trap and then selectively ejects them '
               'by varying the RF and DC potentials.'),
              'mass analyzer type',
              ['ion trap', 'mass analyzer type']),
    Component('linear ion trap', 'MS:1000291',
              ('A two dimensional Paul ion trap in which ions are confined '
               'in the axial dimension by means of an electric field at the '
               'ends of the trap.'),
              'mass analyzer type',
              ['ion trap', 'mass analyzer type']),
    Component('axial ejection linear ion trap', 'MS:1000078',
              ('A linear ion trap mass spectrometer where ions are ejected '
               'along the axis of the analyzer.'),
              'mass analyzer type',
              ['linear ion trap', 'ion trap', 'mass analyzer type']),
    Component('radial ejection linear ion trap', 'MS:1000083',
              ('A linear ion trap mass spectrometer where ions are ejected '
               'along the radius of the analyzer.'),
              'mass analyzer type',
              ['linear ion trap', 'ion trap', 'mass analyzer type']),
])
# [[[end]]]


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('inlet type', term_cls_name="Component", writer=cog.out)
# ]]]
# CV Version: 4.1.55
inlet_types = TermSet([
    Component('continuous flow fast atom bombardment', 'MS:1000055',
              ('Fast atom bombardment ionization in which the analyte in '
               'solution is entrained in a flowing liquid matrix.'),
              'inlet type',
              ['inlet type']),
    Component('direct inlet', 'MS:1000056',
              ('The sample is directly inserted into the ion source, usually '
               'on the end of a heatable probe.'),
              'inlet type',
              ['inlet type']),
    Component('electrospray inlet', 'MS:1000057',
              ('Inlet used for introducing the liquid sample into an '
               'electrospray ionization source.'),
              'inlet type',
              ['inlet type']),
    Component('flow injection analysis', 'MS:1000058',
              ('Sample is directly injected or infused into the ionization '
               'source.'),
              'inlet type',
              ['inlet type']),
    Component('inductively coupled plasma', 'MS:1000059',
              ('A gas discharge ion source in which the energy to the plasma '
               'is supplied by electromagnetic induction.'),
              'inlet type',
              ['inlet type']),
    Component('infusion', 'MS:1000060',
              ('The continuous flow of solution of a sample into the '
               'ionization source.'),
              'inlet type',
              ['inlet type']),
    Component('jet separator', 'MS:1000061',
              ('A device that separates carrier gas from gaseous analyte '
               'molecules on the basis of diffusivity.'),
              'inlet type',
              ['inlet type']),
    Component('membrane separator', 'MS:1000062',
              ('A device to separate carrier molecules from analyte '
               'molecules on the basis of ease of diffusion across a '
               'semipermeable membrane.'),
              'inlet type',
              ['inlet type']),
    Component('moving belt', 'MS:1000063',
              ('Continuous moving surface in the form of a belt which passes '
               'through an ion source carrying analyte molecules.'),
              'inlet type',
              ['inlet type']),
    Component('moving wire', 'MS:1000064',
              ('Continuous moving surface in the form of a wire which passes '
               'through an ion source carrying analyte molecules.'),
              'inlet type',
              ['inlet type']),
    Component('open split', 'MS:1000065',
              ('A division of flowing stream of liquid into two streams.'),
              'inlet type',
              ['inlet type']),
    Component('particle beam', 'MS:1000066',
              ('Method for generating ions from a solution of an analyte.'),
              'inlet type',
              ['inlet type']),
    Component('reservoir', 'MS:1000067',
              ('A sample inlet method involving a reservoir.'),
              'inlet type',
              ['inlet type']),
    Component('septum', 'MS:1000068',
              ('A disc composed of a flexible material that seals the '
               'entrance to the reservoir. Can also be entrance to the '
               'vacuum chamber.'),
              'inlet type',
              ['inlet type']),
    Component('thermospray inlet', 'MS:1000069',
              ('A method for generating gas phase ions from a solution of an '
               'analyte by rapid heating of the sample.'),
              'inlet type',
              ['inlet type']),
    Component('direct insertion probe', 'MS:1000248',
              ('A device for introducing a solid or liquid sample into a '
               'mass spectrometer ion source for desorption ionization.'),
              'inlet type',
              ['inlet type']),
    Component('direct liquid introduction', 'MS:1000249',
              ('The delivery of a liquid sample into a mass spectrometer for '
               'spray or desorption ionization.'),
              'inlet type',
              ['inlet type']),
    Component('membrane inlet', 'MS:1000396',
              ('A semi-permeable membrane separator that permits the passage '
               'of gas sample directly to the mass spectrometer ion source.'),
              'inlet type',
              ['inlet type']),
    Component('nanospray inlet', 'MS:1000485',
              ('Nanospray Inlet.'),
              'inlet type',
              ['electrospray inlet', 'inlet type']),
])
# [[[end]]]


all_components = ionization_types + detector_types + analyzer_types + inlet_types

# more consistent alias
components = all_components

all_components_by_name = {c.name: c for c in all_components}


def component(name):
    try:
        return all_components_by_name[name]
    except KeyError:
        return Component(name, name, name, name, [name])


class InstrumentModel(Term):
    """Describes a named model of a mass spectrometer, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`InstrumentModel` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('instrument model', term_cls_name="InstrumentModel", writer=cog.out)
# ]]]
# CV Version: 4.1.55
instrument_models = TermSet([
    InstrumentModel('SCIEX instrument model', 'MS:1000121',
                    ('The brand of instruments from the joint venture between '
                     'Applied Biosystems and MDS Analytical Technologies (formerly '
                     'MDS SCIEX). Previously branded as \\"Applied Biosystems|MDS '
                     'SCIEX\\".'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Bruker Daltonics instrument model', 'MS:1000122',
                    ("Bruker Daltonics' instrument model."),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Shimadzu instrument model', 'MS:1000124',
                    ('Shimadzu corporation instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Waters instrument model', 'MS:1000126',
                    ('Waters Corporation instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Thermo Fisher Scientific instrument model', 'MS:1000483',
                    ('Thermo Fisher Scientific instrument model. The company has '
                     'gone through several names including Thermo Finnigan, Thermo '
                     'Scientific.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Hitachi instrument model', 'MS:1000488',
                    ('Hitachi instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Varian instrument model', 'MS:1000489',
                    ('Varian instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Agilent instrument model', 'MS:1000490',
                    ('Agilent instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Dionex instrument model', 'MS:1000491',
                    ('Dionex instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('Applied Biosystems instrument model', 'MS:1000495',
                    ('Applied Biosystems instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('LECO instrument model', 'MS:1001800',
                    ('LECO instrument model.'),
                    'instrument model',
                    ['instrument model']),
    InstrumentModel('4000 QTRAP', 'MS:1000139',
                    ('Applied Biosystems/MDS SCIEX Q 4000 TRAP MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 150EX', 'MS:1000143',
                    ('Applied Biosystems/MDS SCIEX API 150EX MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 150EX Prep', 'MS:1000144',
                    ('Applied Biosystems/MDS SCIEX API 150EX Prep MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 2000', 'MS:1000145',
                    ('Applied Biosystems/MDS SCIEX API 2000 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 3000', 'MS:1000146',
                    ('Applied Biosystems/MDS SCIEX API 3000 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 4000', 'MS:1000147',
                    ('Applied Biosystems/MDS SCIEX API 4000 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('proteomics solution 1', 'MS:1000186',
                    ('Applied Biosystems/MDS SCIEX Proteomics Solution 1 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Q TRAP', 'MS:1000187',
                    ('Applied Biosystems/MDS SCIEX Q TRAP MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QSTAR', 'MS:1000190',
                    ('Applied Biosystems/MDS SCIEX QSTAR MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('SymBiot I', 'MS:1000194',
                    ('Applied Biosystems/MDS SCIEX SymBiot I MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('SymBiot XVI', 'MS:1000195',
                    ('Applied Biosystems/MDS SCIEX SymBiot XVI MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('3200 QTRAP', 'MS:1000651',
                    ('SCIEX or Applied Biosystems|MDS SCIEX QTRAP 3200.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('4800 Plus MALDI TOF/TOF', 'MS:1000652',
                    ('SCIEX or Applied Biosystems|MDS SCIEX 4800 Plus MALDI TOF- '
                     'TOF Analyzer.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 3200', 'MS:1000653',
                    ('SCIEX or Applied Biosystems|MDS SCIEX API 3200 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 5000', 'MS:1000654',
                    ('SCIEX or Applied Biosystems|MDS SCIEX API 5000 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QSTAR Elite', 'MS:1000655',
                    ('SCIEX or Applied Biosystems|MDS SCIEX QSTAR Elite.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QSTAR Pulsar', 'MS:1000656',
                    ('Applied Biosystems|MDS SCIEX QSTAR Pulsar.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QSTAR XL', 'MS:1000657',
                    ('Applied Biosystems|MDS SCIEX QSTAR XL.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QTRAP 5500', 'MS:1000931',
                    ('Applied Biosystems|MDS SCIEX QTRAP 5500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('TripleTOF 5600', 'MS:1000932',
                    ('SCIEX TripleTOF 5600, a quadrupole - quadrupole - time-of- '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('5800 TOF/TOF', 'MS:1001482',
                    ('SCIEX 5800 TOF-TOF Analyzer.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('TripleTOF 6600', 'MS:1002533',
                    ('SCIEX TripleTOF 6600, a quadrupole - quadrupole - time-of- '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('2000 QTRAP', 'MS:1002577',
                    ('SCIEX 2000 QTRAP.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('2500 QTRAP', 'MS:1002578',
                    ('SCIEX 2500 QTRAP.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('3500 QTRAP', 'MS:1002579',
                    ('SCIEX 3500 QTRAP.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QTRAP 4500', 'MS:1002580',
                    ('SCIEX QTRAP 4500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QTRAP 6500', 'MS:1002581',
                    ('SCIEX QTRAP 6500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('QTRAP 6500+', 'MS:1002582',
                    ('SCIEX QTRAP 6500+.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('TripleTOF 4600', 'MS:1002583',
                    ('SCIEX TripleTOF 4600 time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('TripleTOF 5600+', 'MS:1002584',
                    ('SCIEX TripleTOF 5600+ time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 100', 'MS:1002585',
                    ('Applied Biosystems/MDS SCIEX API 100 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 100LC', 'MS:1002586',
                    ('Applied Biosystems/MDS SCIEX API 100LC MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 165', 'MS:1002587',
                    ('Applied Biosystems/MDS SCIEX API 165 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 300', 'MS:1002588',
                    ('Applied Biosystems/MDS SCIEX API 300 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 350', 'MS:1002589',
                    ('Applied Biosystems/MDS SCIEX API 350 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('API 365', 'MS:1002590',
                    ('Applied Biosystems/MDS SCIEX API 365 MS.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Triple Quad 3500', 'MS:1002591',
                    ('SCIEX Triple Quad 3500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Triple Quad 4500', 'MS:1002592',
                    ('SCIEX Triple Quad 4500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Triple Quad 5500', 'MS:1002593',
                    ('SCIEX Triple Quad 5500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Triple Quad 6500', 'MS:1002594',
                    ('SCIEX Triple Quad 6500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Triple Quad 6500+', 'MS:1002595',
                    ('SCIEX Triple Quad 6500+.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('X500R QTOF', 'MS:1002674',
                    ('SCIEX X500R QTOF, a quadrupole - quadrupole - time-of-flight '
                     'mass spectrometer.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Triple Quad 7500', 'MS:1003144',
                    ('SCIEX Triple Quad 7500.'),
                    'instrument model',
                    ['SCIEX instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics HCT Series', 'MS:1000697',
                    ("Bruker Daltonics' HCT Series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics esquire series', 'MS:1001533',
                    ("Bruker Daltonics' esquire series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics flex series', 'MS:1001534',
                    ("Bruker Daltonics' flex series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics BioTOF series', 'MS:1001535',
                    ("Bruker Daltonics' BioTOF series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics micrOTOF series', 'MS:1001536',
                    ("Bruker Daltonics' micrOTOF series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics amaZon series', 'MS:1001545',
                    ("Bruker Daltonics' amaZon series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics maXis series', 'MS:1001547',
                    ("Bruker Daltonics' maXis series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics solarix series', 'MS:1001548',
                    ("Bruker Daltonics' solarix: ESI quadrupole ion trap, APCI, "
                     'APPI, ETD, PTR.'),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics apex series', 'MS:1001556',
                    ("Bruker Daltonics' apex series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics SCION series', 'MS:1002293',
                    ("Bruker Daltonics' SCION series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics EVOQ series', 'MS:1002294',
                    ("Bruker Daltonics' EVOQ series."),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Bruker Daltonics timsTOF series', 'MS:1003123',
                    ('Bruker Daltonics timsTOF series'),
                    'instrument model',
                    ['Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('Shimadzu MALDI-TOF instrument model', 'MS:1000602',
                    ('Shimadzu MALDI-TOF instrument model.'),
                    'instrument model',
                    ['Shimadzu instrument model', 'instrument model']),
    InstrumentModel('Shimadzu Scientific Instruments instrument model', 'MS:1000603',
                    ('Shimadzu Scientific Instruments instrument model.'),
                    'instrument model',
                    ['Shimadzu instrument model', 'instrument model']),
    InstrumentModel('Auto Spec Ultima NT', 'MS:1000150',
                    ('Waters magnetic sector based AutoSpec Ultima NT MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('GCT', 'MS:1000159',
                    ('Waters oa-ToF based GCT.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('IsoPrime', 'MS:1000164',
                    ('Waters IsoPrime MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('IsoProbe', 'MS:1000165',
                    ('Waters IsoProbe MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('IsoProbe T', 'MS:1000166',
                    ('Waters IsoProbe T MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('M@LDI L', 'MS:1000170',
                    ('Waters oa-ToF based MALDI L.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('M@LDI LR', 'MS:1000171',
                    ('Waters oa-ToF based MALDI LR.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('NG-5400', 'MS:1000180',
                    ('Waters NG-5400 MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Platform ICP', 'MS:1000184',
                    ('Waters Platform ICP MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Q-Tof micro', 'MS:1000188',
                    ('Waters oa-ToF based Q-Tof micro.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Q-Tof Ultima', 'MS:1000189',
                    ('Waters oa-ToF based Q-Tof Ultima.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('quattro micro', 'MS:1000191',
                    ('Waters (triple) quadrupole based micro.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Quattro Ultima', 'MS:1000192',
                    ('Waters (triple) quadrupole based Ultima.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Q-Tof Premier', 'MS:1000632',
                    ('Waters oa-ToF based Q-Tof Premier.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Acquity UPLC PDA', 'MS:1000818',
                    ('Acquity UPLC Photodiode Array Detector.'),
                    'instrument model',
                    ['Waters instrument model', 'photodiode array detector', 'instrument model', 'array detector', 'detector type']),
    InstrumentModel('Acquity UPLC FLR', 'MS:1000819',
                    ('Acquity UPLC Fluorescence Detector.'),
                    'instrument model',
                    ['Waters instrument model', 'fluorescence detector', 'instrument model', 'detector type']),
    InstrumentModel('ACQUITY UPLC', 'MS:1001761',
                    ('Waters LC-system ACQUITY UPLC.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('ACQUITY UPLC H-Class', 'MS:1001762',
                    ('Waters LC-system ACQUITY UPLC H-Class.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('ACQUITY UPLC H-Class Bio', 'MS:1001763',
                    ('Waters LC-system ACQUITY UPLC H-Class Bio.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('ACQUITY UPLC I-Class', 'MS:1001764',
                    ('Waters LC-system ACQUITY UPLC I-Class.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('ACQUITY UPLC Systems with 2D Technology', 'MS:1001765',
                    ('Waters LC-system ACQUITY UPLC Systems with 2D Technology.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('nanoACQUITY UPLC', 'MS:1001766',
                    ('Waters LC-system nanoACQUITY UPLC.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('nanoACQUITY UPLC System with 1D Technology', 'MS:1001767',
                    ('Waters LC-system nanoACQUITY UPLC System with 1D Technology.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('nanoACQUITY UPLC with HDX Technology', 'MS:1001768',
                    ('Waters LC-system nanoACQUITY UPLC with HDX Technology.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('TRIZAIC UPLC nanoTile', 'MS:1001769',
                    ('Waters LC-system TRIZAIC UPLC nanoTile.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('GCT Premier', 'MS:1001770',
                    ('Waters oa-ToF based GCT Premier.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('MALDI Synapt G2 HDMS', 'MS:1001771',
                    ('Waters oa-ToF based MALDI Synapt G2 HDMS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('MALDI Synapt G2 MS', 'MS:1001772',
                    ('Waters oa-ToF based MALDI Synapt G2 MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('MALDI Synapt G2-S HDMS', 'MS:1001773',
                    ('Waters oa-ToF based MALDI Synapt G2 MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('MALDI Synapt G2-S MS', 'MS:1001774',
                    ('Waters oa-ToF based MALDI Synapt G2-S MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('MALDI Synapt HDMS', 'MS:1001775',
                    ('Waters oa-ToF based MALDI Synapt HDMS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('MALDI Synapt MS', 'MS:1001776',
                    ('Waters oa-ToF based MALDI Synapt MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Synapt G2 HDMS', 'MS:1001777',
                    ('Waters oa-ToF based Synapt G2 HDMS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Synapt G2 MS', 'MS:1001778',
                    ('Waters oa-ToF based Synapt G2 MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Synapt G2-S HDMS', 'MS:1001779',
                    ('Waters oa-ToF based Synapt G2-S HDMS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Synapt G2-S MS', 'MS:1001780',
                    ('Waters oa-ToF based Synapt G2-S MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Synapt HDMS', 'MS:1001781',
                    ('Waters oa-ToF based Synapt HDMS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Synapt MS', 'MS:1001782',
                    ('Waters oa-ToF based Synapt MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo G2 Q-Tof', 'MS:1001783',
                    ('Waters oa-ToF based Xevo G2 Q-Tof.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo G2 Tof', 'MS:1001784',
                    ('Waters oa-ToF based Xevo G2 Tof.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo Q-Tof', 'MS:1001785',
                    ('Waters oa-ToF based Xevo Q-Tof.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('3100', 'MS:1001786',
                    ('Waters quadrupole based 3100.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Acquity SQD', 'MS:1001787',
                    ('Waters quadrupole based Acquity SQD.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Acquity TQD', 'MS:1001788',
                    ('Waters quadrupole based Acquity TQD.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Quattro micro GC', 'MS:1001789',
                    ('Waters (triple) quadrupole based Quattro micro GC.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo TQ MS', 'MS:1001790',
                    ('Waters quadrupole based Xevo TQ MS.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo TQD', 'MS:1001791',
                    ('Waters quadrupole based Xevo TQD.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo TQ-S', 'MS:1001792',
                    ('Waters quadrupole based Xevo TQ-S.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('SQ Detector 2', 'MS:1002274',
                    ('Waters quadrupole based SQ Detector 2.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo G2-S Tof', 'MS:1002275',
                    ('Waters oa-ToF based Xevo G2-S Tof.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo G2-S QTof', 'MS:1002276',
                    ('Waters oa-ToF based Xevo G2-S QTof.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('AutoSpec Premier', 'MS:1002277',
                    ('Waters AutoSpec Premier magnetic sector instrument.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('SYNAPT G2-Si', 'MS:1002726',
                    ('Waters Corporation SYNAPT G2-Si orthogonal acceleration '
                     'time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('MALDI SYNAPT G2-Si', 'MS:1002727',
                    ('Waters Corporation MALDI SYNAPT G2-Si orthogonal '
                     'acceleration time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Vion IMS QTof', 'MS:1002728',
                    ('Waters Corporation Vion IMS QTof orthogonal acceleration '
                     'time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo G2 XS Tof', 'MS:1002729',
                    ('Waters Corporation Xevo G2 XS Tof orthogonal acceleration '
                     'time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo TQ-XS', 'MS:1002730',
                    ('Waters Corporation Xevo TQ-XS triple quadrupole mass '
                     'spectrometer.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Xevo TQ-S micro', 'MS:1002731',
                    ('Waters Corporation Xevo TQ-S micro triple quadrupole mass '
                     'spectrometer.'),
                    'instrument model',
                    ['Waters instrument model', 'instrument model']),
    InstrumentModel('Thermo Finnigan instrument model', 'MS:1000125',
                    ('ThermoFinnigan from Thermo Electron Corporation instrument '
                     'model.'),
                    'instrument model',
                    ['Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Thermo Electron instrument model', 'MS:1000492',
                    ('Thermo Electron Corporation instrument model.'),
                    'instrument model',
                    ['Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Finnigan MAT instrument model', 'MS:1000493',
                    ('Finnigan MAT instrument model.'),
                    'instrument model',
                    ['Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Thermo Scientific instrument model', 'MS:1000494',
                    ('Thermo Scientific instrument model.'),
                    'instrument model',
                    ['Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('IonSpec instrument model', 'MS:1000123',
                    ('IonSpec corporation instrument model.'),
                    'instrument model',
                    ['Varian instrument model', 'instrument model']),
    InstrumentModel('1200 series LC/MSD SL', 'MS:1000467',
                    ('The 1200 Series LC/MSD SL ion trap belongs to the Agilent '
                     'LC/MSD ion trap family. It provides fast polarity switching '
                     'and multisignal data acquisition capabilities in a single '
                     'run while also providing 5 stages of automated data '
                     'dependent MS2 and 11 stages of manual MS2.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6110 Quadrupole LC/MS', 'MS:1000468',
                    ('The 6110 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with an entry level '
                     'single quadrupole mass spectrometer from the 6100 Series of '
                     'Agilent quadrupole mass spectrometers. 6110 Quadrupole mass '
                     'spectrometer has m/z range of 10-1500 and 2500 u/s scan '
                     'speed. It proves useful for wide range of SIM quantitative '
                     'applications.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6120A Quadrupole LC/MS', 'MS:1000469',
                    ('The 6120A Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a single quadrupole '
                     'mass spectrometer from the 6100 Series of Agilent mass '
                     'spectrometers. 6120 quadrupole mass spectrometer has m/z '
                     'range of 10-1500, 2500 u/s scan speed and utilizes multiple '
                     'signal acquisition.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6130 Quadrupole LC/MS', 'MS:1000470',
                    ('The 6130 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a single quadrupole '
                     'mass spectrometer from the 6100 series of Agilent mass '
                     'spectrometers. The 6130 quadrupole mass spectrometer has m/z '
                     'range of 2-3000, 2500 u/s scan speed in standard mode and '
                     '5250 u/s speed in fast-scan mode. It also uses multiple '
                     'signal acquisition.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6140 Quadrupole LC/MS', 'MS:1000471',
                    ('The 6140 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a single quadrupole '
                     'mass spectrometer from the 6100 Series of Agilent quadrupole '
                     'mass spectrometers. 6140 Quadrupole mass spectrometer has '
                     'm/z range of 10-1350, 2500 u/s scan speed in standard mode '
                     'and 10000 u/s speed in fast-scan mode. It also uses multiple '
                     'signal acquisition.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6210 Time-of-Flight LC/MS', 'MS:1000472',
                    ('The 6210 Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer. This time of flight mass '
                     'spectrometer has a m/z range of 50-12000, mass accuracy of '
                     'less than 2 ppm and resolution greater than 13,000 at m/z '
                     '2722. It has multiple ion sources and can be used with '
                     'multimode ion sources.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6310 Ion Trap LC/MS', 'MS:1000473',
                    ('The 6310 Ion Trap LC/MS is a Agilent liquid chromatography '
                     'instrument combined with a 6300 series Agilent ion trap. It '
                     'has a mass range of 50-2200 between 0.6 to 0.35 resolution '
                     'and mass range of 200-4000 with resolution of 3-4. The scan '
                     'speed varies from 1650-27000 for the respective mass ranges.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6320 Ion Trap LC/MS', 'MS:1000474',
                    ('The 6320 Ion Trap LC/MS is a Agilent liquid chromatography '
                     'instrument combined with a 6300 series Agilent ion trap. It '
                     'has a mass range of 50-2200 between 0.6 to 0.25 resolution '
                     'and mass range of 200-4000 with resolution of less than 3. '
                     'The scan speed varies from 1650-27000 for the respective '
                     'mass ranges.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6330 Ion Trap LC/MS', 'MS:1000475',
                    ('The 6330 Ion Trap LC/MS is a Agilent liquid chromatography '
                     'instrument combined with a 6300 series Agilent ion trap. It '
                     'has a mass range of 50-2200 between 0.6 to 0.25 resolution '
                     'and mass range of 200-4000 with resolution of less than 3. '
                     'The scan speed varies from 1650-27000 for the respective '
                     'mass ranges.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6340 Ion Trap LC/MS', 'MS:1000476',
                    ('The 6340 Ion Trap LC/MS is a Agilent liquid chromatography '
                     'instrument combined with a 6300 series Agilent ion trap. It '
                     'has a mass range of 50-2200 between 0.6 to 0.25 resolution '
                     'and mass range of 200-4000 with resolution of less than 3. '
                     'The scan speed varies from 1650-27000 for the respective '
                     'mass ranges.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6410 Triple Quadrupole LC/MS', 'MS:1000477',
                    ('The 6410 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer. Mass range of the mass '
                     'spectrometer is 15-1650 m/z, resolution is at three settings '
                     'of 0.7 u (unit), 1.2 u (wide) and 2.5 u (widest). The mass '
                     'accuracy for 6410 mass spectrometer is 0.1 across the mass '
                     'range. The collision cell is a hexapole with linear '
                     'acceleration.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('1200 series LC/MSD VL', 'MS:1000478',
                    ('The LC/MSD VL ion trap is part of the family of Agilent ion '
                     'trap mass spectrometers. It has ESI, APCI and APPI ion '
                     'sources and is a useful ion trap when the amount of sample '
                     'is not the limiting factor.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6220 Time-of-Flight LC/MS', 'MS:1000675',
                    ('The 6220 Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer. This time of flight mass '
                     'spectrometer has a m/z range of 50-12000, mass accuracy of '
                     'less than 2 ppm and resolution greater than 13,000 at m/z '
                     '2722. It has multiple ion sources and can be used with '
                     'multimode ion sources.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6510 Quadrupole Time-of-Flight LC/MS', 'MS:1000676',
                    ('The 6510 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer. This time of flight mass '
                     'spectrometer has a m/z range of 50-12000, mass accuracy of '
                     'less than 2 ppm and resolution greater than 13,000 at m/z '
                     '2722. It has multiple ion sources and can be used with '
                     'multimode ion sources.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6520A Quadrupole Time-of-Flight LC/MS', 'MS:1000677',
                    ('The 6520A Quadrupole Time-of-Flight LC/MS is a Agilent '
                     'liquid chromatography instrument combined with a Agilent '
                     'time of flight mass spectrometer. This time of flight mass '
                     'spectrometer has a m/z range of 50-12000, mass accuracy of '
                     'less than 2 ppm and resolution greater than 26,000 at m/z '
                     '2722. It has multiple ion sources and can be used with '
                     'multimode ion sources.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6420 Triple Quadrupole LC/MS', 'MS:1002444',
                    ('The 6420 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6460 Triple Quadrupole LC/MS', 'MS:1002445',
                    ('The 6460 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer. It is similar to the 6420 but '
                     'adds Agilent Jet Stream (AJS) technology to increase '
                     'sensitivity.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6490 Triple Quadrupole LC/MS', 'MS:1002446',
                    ('The 6490 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer. It is similar to the 6420 but '
                     'adds the Agilent iFunnel technology to increase sensitivity.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6550 iFunnel Q-TOF LC/MS', 'MS:1002783',
                    ('The 6550 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6550A iFunnel Q-TOF LC/MS', 'MS:1002784',
                    ('The 6550A Quadrupole Time-of-Flight LC/MS is a Agilent '
                     'liquid chromatography instrument combined with a Agilent '
                     'time of flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6520B Q-TOF LC/MS', 'MS:1002785',
                    ('The 6520B Quadrupole Time-of-Flight LC/MS is a Agilent '
                     'liquid chromatography instrument combined with a Agilent '
                     'time of flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6530A Q-TOF LC/MS', 'MS:1002786',
                    ('The 6530A Quadrupole Time-of-Flight LC/MS is a Agilent '
                     'liquid chromatography instrument combined with a Agilent '
                     'time of flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6530B Q-TOF LC/MS', 'MS:1002787',
                    ('The 6530B Quadrupole Time-of-Flight LC/MS is a Agilent '
                     'liquid chromatography instrument combined with a Agilent '
                     'time of flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6538 Q-TOF LC/MS', 'MS:1002788',
                    ('The 6538 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6540 Q-TOF LC/MS', 'MS:1002789',
                    ('The 6540 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6542 Q-TOF LC/MS', 'MS:1002790',
                    ('The 6542 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6545 Q-TOF LC/MS', 'MS:1002791',
                    ('The 6545 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6560 Q-TOF LC/MS', 'MS:1002792',
                    ('The 6560 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6570 Q-TOF LC/MS', 'MS:1002793',
                    ('The 6570 Quadrupole Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6120B Quadrupole LC/MS', 'MS:1002794',
                    ('The 6120B Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a single quadrupole '
                     'mass spectrometer from the 6100 Series of Agilent mass '
                     'spectrometers.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6150 Quadrupole LC/MS', 'MS:1002795',
                    ('The 6150 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a single quadrupole '
                     'mass spectrometer from the 6100 Series of Agilent mass '
                     'spectrometers.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6224 Time-of-Flight LC/MS', 'MS:1002796',
                    ('The 6224 Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6230A Time-of-Flight LC/MS', 'MS:1002797',
                    ('The 6230A Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6230B Time-of-Flight LC/MS', 'MS:1002798',
                    ('The 6230B Time-of-Flight LC/MS is a Agilent liquid '
                     'chromatography instrument combined with a Agilent time of '
                     'flight mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6430 Triple Quadrupole LC/MS', 'MS:1002799',
                    ('The 6430 Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6495A Triple Quadrupole LC/MS', 'MS:1002800',
                    ('The 6495A Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('6495B Triple Quadrupole LC/MS', 'MS:1002801',
                    ('The 6495B Quadrupole LC/MS system is a Agilent liquid '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('7000A Triple Quadrupole GC/MS', 'MS:1002802',
                    ('The 7000A Quadrupole GC/MS system is a Agilent gas '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('7000B Triple Quadrupole GC/MS', 'MS:1002803',
                    ('The 7000B Quadrupole GC/MS system is a Agilent gas '
                     'chromatography instrument combined with a Agilent triple '
                     'quadrupole mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('7800 Quadrupole ICP-MS', 'MS:1002804',
                    ('The 7800 Quadrupole ICP-MS system is a Agilent inductively '
                     'couple plasma instrument combined with a Agilent quadrupole '
                     'mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('8800 Triple Quadrupole ICP-MS', 'MS:1002805',
                    ('The 8800 Quadrupole ICP-MS system is a Agilent inductively '
                     'couple plasma instrument combined with a Agilent quadrupole '
                     'mass spectrometer.'),
                    'instrument model',
                    ['Agilent instrument model', 'instrument model']),
    InstrumentModel('4700 Proteomics Analyzer', 'MS:1000140',
                    ('Applied Biosystems/MDS SCIEX 4700 Proteomics Analyzer MS.'),
                    'instrument model',
                    ['Applied Biosystems instrument model', 'instrument model']),
    InstrumentModel('Voyager-DE PRO', 'MS:1000203',
                    ('Applied Biosystems/MDS SCIEX Voyager-DE PRO MS.'),
                    'instrument model',
                    ['Applied Biosystems instrument model', 'instrument model']),
    InstrumentModel('Voyager-DE STR', 'MS:1000204',
                    ('Applied Biosystems/MDS SCIEX Voyager-DE STR MS.'),
                    'instrument model',
                    ['Applied Biosystems instrument model', 'instrument model']),
    InstrumentModel('4800 Proteomics Analyzer', 'MS:1000658',
                    ('Applied Biosystems|MDS SCIEX 4800 Proteomics Analyzer.'),
                    'instrument model',
                    ['Applied Biosystems instrument model', 'instrument model']),
    InstrumentModel('Pegasus HRT', 'MS:1001801',
                    ('LECO high resolution time-of-flight GC mass spectrometer.'),
                    'instrument model',
                    ['LECO instrument model', 'instrument model']),
    InstrumentModel('Citius HRT', 'MS:1001802',
                    ('LECO high resolution time-of-flight LC mass spectrometer.'),
                    'instrument model',
                    ['LECO instrument model', 'instrument model']),
    InstrumentModel('Pegasus', 'MS:1001803',
                    ('LECO GC time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['LECO instrument model', 'instrument model']),
    InstrumentModel('TruTOF', 'MS:1001804',
                    ('LECO bench-top GC time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['LECO instrument model', 'instrument model']),
    InstrumentModel('Pegasus 4D', 'MS:1001945',
                    ('LECO nominal mass resolution time-of-flight GCxGC mass '
                     'spectrometer.'),
                    'instrument model',
                    ['LECO instrument model', 'instrument model']),
    InstrumentModel('Pegasus III', 'MS:1002278',
                    ('LECO nominal mass resolution time-of-flight GC mass '
                     'spectrometer.'),
                    'instrument model',
                    ['LECO instrument model', 'instrument model']),
    InstrumentModel('Pegasus BT', 'MS:1002719',
                    ('LECO bench-top GC time-of-flight mass spectrometer.'),
                    'instrument model',
                    ['LECO instrument model', 'instrument model']),
    InstrumentModel('HCT', 'MS:1000160',
                    ("Bruker Daltonics' HCT: ESI Q-TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    ['Bruker Daltonics HCT Series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('HCTplus', 'MS:1000161',
                    ("Bruker Daltonics' HCTplus: ESI Q-TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    ['Bruker Daltonics HCT Series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('HCTultra', 'MS:1000698',
                    ("Bruker Daltonics' HCTultra: ESI TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    ['Bruker Daltonics HCT Series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('HCTultra PTM', 'MS:1000699',
                    ("Bruker Daltonics' HCTultra PTM: ESI TOF, Nanospray, APCI, "
                     'APPI, PTR.'),
                    'instrument model',
                    ['Bruker Daltonics HCT Series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('HCTultra ETD II', 'MS:1000700',
                    ("Bruker Daltonics' HCTultra ETD II: ESI Q-TOF, Nanospray, "
                     'APCI, APPI, ETD.'),
                    'instrument model',
                    ['Bruker Daltonics HCT Series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('esquire 4000', 'MS:1000156',
                    ("Bruker Daltonics' esquire 4000: linear ion trap, ESI, MALDI, "
                     'Nanospray, APCI, APPI.'),
                    'instrument model',
                    ['Bruker Daltonics esquire series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('esquire 6000', 'MS:1000157',
                    ("Bruker Daltonics' esquire 6000: linear ion trap, ESI, MALDI, "
                     'Nanospray, APCI, APPI.'),
                    'instrument model',
                    ['Bruker Daltonics esquire series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('autoflex II', 'MS:1000148',
                    ("Bruker Daltonics' autoflex II: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('autoflex TOF/TOF', 'MS:1000149',
                    ("Bruker Daltonics' autoflex TOF/TOF MS: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('microflex', 'MS:1000177',
                    ("Bruker Daltonics' microflex: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('OmniFlex', 'MS:1000183',
                    ("Bruker Daltonics' OmniFlex: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('ultraflex', 'MS:1000201',
                    ("Bruker Daltonics' ultraflex: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('ultraflex TOF/TOF', 'MS:1000202',
                    ("Bruker Daltonics' ultraflex TOF/TOF: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('autoflex III smartbeam', 'MS:1000696',
                    ("Bruker Daltonics' autoflex III smartbeam: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('microflex LT', 'MS:1000701',
                    ("Bruker Daltonics' microflex LT: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('ultraflex III TOF/TOF', 'MS:1000705',
                    ("Bruker Daltonics' ultraflex III TOF/TOF: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('microflex LRF', 'MS:1001543',
                    ("Bruker Daltonics' microflex LRF: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('ultrafleXtreme', 'MS:1001544',
                    ("Bruker Daltonics' ultrafleXtreme: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('microflex II', 'MS:1001550',
                    ("Bruker Daltonics' microflex II: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('autoflex II TOF/TOF', 'MS:1001553',
                    ("Bruker Daltonics' autoflex II TOF/TOF: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('autoflex III TOF/TOF smartbeam', 'MS:1001554',
                    ("Bruker Daltonics' autoflex III TOF/TOF smartbeam: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('autoflex', 'MS:1001555',
                    ("Bruker Daltonics' autoflex: MALDI TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('rapifleX', 'MS:1003122',
                    ("Bruker Daltonics' rapiflex: MALDI TOF/TOF."),
                    'instrument model',
                    ['Bruker Daltonics flex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('BioTOF II', 'MS:1000151',
                    ("Bruker Daltonics' BioTOF II: ESI TOF."),
                    'instrument model',
                    ['Bruker Daltonics BioTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('BioTOF-Q', 'MS:1000152',
                    ("Bruker Daltonics' BioTOF-Q: ESI Q-TOF."),
                    'instrument model',
                    ['Bruker Daltonics BioTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('BioTOF', 'MS:1001537',
                    ("Bruker Daltonics' BioTOF: ESI TOF."),
                    'instrument model',
                    ['Bruker Daltonics BioTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('BioTOF III', 'MS:1001538',
                    ("Bruker Daltonics' BioTOF III: ESI TOF."),
                    'instrument model',
                    ['Bruker Daltonics BioTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('UltroTOF-Q', 'MS:1001539',
                    ("Bruker Daltonics' UltroTOF-Q: ESI Q-TOF (MALDI optional)."),
                    'instrument model',
                    ['Bruker Daltonics BioTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('microTOF LC', 'MS:1000178',
                    ("Bruker Daltonics' microTOF LC: ESI TOF, Nanospray, APCI, "
                     'APPI.'),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('micrOTOF', 'MS:1000702',
                    ("Bruker Daltonics' micrOTOF: ESI TOF, APCI, APPI."),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('micrOTOF-Q', 'MS:1000703',
                    ("Bruker Daltonics' micrOTOF-Q: ESI Q-TOF, Nanospray, APCI, "
                     'APPI.'),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('micrOTOF-Q II', 'MS:1000704',
                    ("Bruker Daltonics' micrOTOF-Q II: ESI Q-TOF, Nanospray, APCI, "
                     'APPI.'),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('micrOTOF II', 'MS:1001540',
                    ("Bruker Daltonics' micrOTOF II: ESI TOF, Nanospray, APCI, "
                     'APPI.'),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('impact', 'MS:1002077',
                    ("Bruker Daltonics' impact: ESI Q-TOF, Nanospray, APCI, APPI, "
                     'GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('compact', 'MS:1002280',
                    ("Bruker Daltonics' compact: ESI Q-TOF, Nanospray, APCI, APPI, "
                     'GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('micrOTOF-Q III', 'MS:1002299',
                    ("Bruker Daltonics' micrOTOF-Q III: ESI Q-TOF, Nanospray, "
                     'APCI, APPI, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('impact II', 'MS:1002666',
                    ("Bruker Daltonics' impact II."),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('impact HD', 'MS:1002667',
                    ("Bruker Daltonics' impact HD."),
                    'instrument model',
                    ['Bruker Daltonics micrOTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('amaZon ETD', 'MS:1001542',
                    ("Bruker Daltonics' amaZon ETD: ESI quadrupole ion trap, "
                     'Nanospray, APCI, APPI, ETD, PTR.'),
                    'instrument model',
                    ['Bruker Daltonics amaZon series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('amaZon X', 'MS:1001546',
                    ("Bruker Daltonics' amaZon X: ESI quadrupole ion trap, APCI, "
                     'APPI, ETD, PTR.'),
                    'instrument model',
                    ['Bruker Daltonics amaZon series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('amaZon Speed ETD', 'MS:1002300',
                    ("Bruker Daltonics' amaZon Speed ETD: ESI quadrupole ion trap, "
                     'Nanospray, APCI, APPI, ETD, PTR, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    ['Bruker Daltonics amaZon series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('amaZon Speed', 'MS:1002301',
                    ("Bruker Daltonics' amaZon ETD: ESI quadrupole ion trap, "
                     'Nanospray, APCI, APPI, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    ['Bruker Daltonics amaZon series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('maXis', 'MS:1001541',
                    ("Bruker Daltonics' maXis: ESI Q-TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    ['Bruker Daltonics maXis series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('maXis 4G', 'MS:1002279',
                    ("Bruker Daltonics' maXis 4G: ESI Q-TOF, Nanospray, APCI, "
                     'APPI, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    ['Bruker Daltonics maXis series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('maXis II', 'MS:1003004',
                    ("Bruker Daltonics' maXis II."),
                    'instrument model',
                    ['Bruker Daltonics maXis series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('solariX', 'MS:1001549',
                    ("Bruker Daltonics' solariX: ESI, MALDI, Qh-FT_ICR."),
                    'instrument model',
                    ['Bruker Daltonics solarix series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('apex IV', 'MS:1000141',
                    ("Bruker Daltonics' apex IV: ESI, MALDI, Nanospray, APCI, "
                     'APPI, Qh-FT_ICR.'),
                    'instrument model',
                    ['Bruker Daltonics apex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('apex Q', 'MS:1000142',
                    ("Bruker Daltonics' apex Q: ESI, MALDI, Nanospray, APCI, APPI, "
                     'Qh-FT_ICR.'),
                    'instrument model',
                    ['Bruker Daltonics apex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('apex ultra', 'MS:1000695',
                    ("Bruker Daltonics' apex ultra: ESI, MALDI, Nanospray, APCI, "
                     'APPI, Qh-FT_ICR.'),
                    'instrument model',
                    ['Bruker Daltonics apex series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('SCION SQ', 'MS:1002295',
                    ("Bruker Daltonics' SCION SQ: GC-single quadrupole."),
                    'instrument model',
                    ['Bruker Daltonics SCION series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('SCION TQ', 'MS:1002296',
                    ("Bruker Daltonics' SCION TQ: GC-triple quadrupole."),
                    'instrument model',
                    ['Bruker Daltonics SCION series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('EVOQ Elite', 'MS:1002297',
                    ("Bruker Daltonics' EVOQ Elite: LC-triple quadrupole."),
                    'instrument model',
                    ['Bruker Daltonics EVOQ series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('EVOQ Qube', 'MS:1002298',
                    ("Bruker Daltonics' EVOQ Qube: LC-triple quadrupole."),
                    'instrument model',
                    ['Bruker Daltonics EVOQ series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('timsTOF Pro', 'MS:1003005',
                    ("Bruker Daltonics' timsTOF Pro."),
                    'instrument model',
                    ['Bruker Daltonics timsTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('timsTOF fleX', 'MS:1003124',
                    ("Bruker Daltonics' timsTOF fleX"),
                    'instrument model',
                    ['Bruker Daltonics timsTOF series', 'Bruker Daltonics instrument model', 'instrument model']),
    InstrumentModel('AXIMA CFR MALDI-TOF', 'MS:1000607',
                    ('Shimadzu Biotech AXIMA CFR MALDI-TOF MS.'),
                    'instrument model',
                    ['Shimadzu MALDI-TOF instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('AXIMA-QIT', 'MS:1000608',
                    ('Shimadzu Biotech AXIMA-QIT MS.'),
                    'instrument model',
                    ['Shimadzu MALDI-TOF instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('AXIMA-CFR plus', 'MS:1000609',
                    ('Shimadzu Biotech AXIMA-CFR plus MS.'),
                    'instrument model',
                    ['Shimadzu MALDI-TOF instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('AXIMA Performance MALDI-TOF/TOF', 'MS:1000610',
                    ('Shimadzu Biotech AXIMA Performance MALDI-TOF/TOF MS.'),
                    'instrument model',
                    ['Shimadzu MALDI-TOF instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('AXIMA Confidence MALDI-TOF', 'MS:1000611',
                    ('Shimadzu Biotech AXIMA Confidence MALDI-TOF (curved field '
                     'reflectron) MS.'),
                    'instrument model',
                    ['Shimadzu MALDI-TOF instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('AXIMA Assurance Linear MALDI-TOF', 'MS:1000612',
                    ('Shimadzu Biotech AXIMA Assurance Linear MALDI-TOF MS.'),
                    'instrument model',
                    ['Shimadzu MALDI-TOF instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('Shimadzu MALDI-7090', 'MS:1002382',
                    ('Shimadzu MALDI-7090: MALDI-TOF-TOF.'),
                    'instrument model',
                    ['Shimadzu MALDI-TOF instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-IT-TOF', 'MS:1000604',
                    ('Shimadzu Scientific Instruments LCMS-IT-TOF MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-2010EV', 'MS:1000605',
                    ('Shimadzu Scientific Instruments LCMS-2010EV MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-2010A', 'MS:1000606',
                    ('Shimadzu Scientific Instruments LCMS-2010A MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-9030', 'MS:1002998',
                    ('Shimadzu Scientific Instruments LCMS-9030 Q-TOF MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-8060', 'MS:1002999',
                    ('Shimadzu Scientific Instruments LCMS-8060 MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-8050', 'MS:1003000',
                    ('Shimadzu Scientific Instruments LCMS-8050 MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-8045', 'MS:1003001',
                    ('Shimadzu Scientific Instruments LCMS-8045 MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-8040', 'MS:1003002',
                    ('Shimadzu Scientific Instruments LCMS-8040 MS.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('LCMS-2020', 'MS:1003003',
                    ('Shimadzu Scientific Instruments LCMS-2020.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('GCMS-QP2010SE', 'MS:1003152',
                    ('Shimadzu Scientific Instruments GCMS-QP2010SE.'),
                    'instrument model',
                    ['Shimadzu Scientific Instruments instrument model', 'Shimadzu instrument model', 'instrument model']),
    InstrumentModel('DELTA plusAdvantage', 'MS:1000153',
                    ('ThermoFinnigan DELTA plusAdvantage MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('DELTAplusXP', 'MS:1000154',
                    ('ThermoFinnigan DELTAplusXP MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LCQ Advantage', 'MS:1000167',
                    ('ThermoFinnigan LCQ Advantage MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LCQ Classic', 'MS:1000168',
                    ('ThermoFinnigan LCQ Classic MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LCQ Deca XP Plus', 'MS:1000169',
                    ('ThermoFinnigan LCQ Deca XP Plus MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('neptune', 'MS:1000179',
                    ('ThermoFinnigan NEPTUNE MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('PolarisQ', 'MS:1000185',
                    ('ThermoFinnigan PolarisQ MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Surveyor MSQ', 'MS:1000193',
                    ('ThermoFinnigan Surveyor MSQ MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TEMPUS TOF', 'MS:1000196',
                    ('ThermoFinnigan TEMPUS TOF MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TRACE DSQ', 'MS:1000197',
                    ('ThermoFinnigan TRACE DSQ MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TRITON', 'MS:1000198',
                    ('ThermoFinnigan TRITON MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Quantum', 'MS:1000199',
                    ('ThermoFinnigan TSQ Quantum MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LCQ Deca', 'MS:1000554',
                    ('ThermoFinnigan LCQ Deca.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('GC Quantum', 'MS:1000558',
                    ('GC Quantum.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LCQ Fleet', 'MS:1000578',
                    ('LCQ Fleet.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('DSQ', 'MS:1000634',
                    ('ThermoFinnigan DSQ GC-MS.'),
                    'instrument model',
                    ['Thermo Finnigan instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('MAT253', 'MS:1000172',
                    ('ThermoFinnigan MAT253 MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('MAT900XP', 'MS:1000173',
                    ('ThermoFinnigan MAT900XP MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('MAT900XP Trap', 'MS:1000174',
                    ('ThermoFinnigan MAT900XP Trap MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('MAT95XP', 'MS:1000175',
                    ('ThermoFinnigan MAT95XP MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('MAT95XP Trap', 'MS:1000176',
                    ('ThermoFinnigan MAT95XP Trap MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('SSQ 7000', 'MS:1000748',
                    ('ThermoFinnigan SSQ 7000 MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ 7000', 'MS:1000749',
                    ('ThermoFinnigan TSQ 7000 MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ', 'MS:1000750',
                    ('ThermoFinnigan TSQ MS.'),
                    'instrument model',
                    ['Finnigan MAT instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ', 'MS:1000447',
                    ('Finnigan LTQ MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ FT', 'MS:1000448',
                    ('Finnigan LTQ FT MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap', 'MS:1000449',
                    ('Finnigan LTQ Orbitrap MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LXQ', 'MS:1000450',
                    ('Finnigan LXQ MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap Discovery', 'MS:1000555',
                    ('LTQ Orbitrap Discovery.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap XL', 'MS:1000556',
                    ('LTQ Orbitrap XL.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ FT Ultra', 'MS:1000557',
                    ('LTQ FT Ultra.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Surveyor PDA', 'MS:1000622',
                    ('Surveyor PDA.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Accela PDA', 'MS:1000623',
                    ('Accela PDA.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('ITQ 700', 'MS:1000635',
                    ('Thermo Scientific ITQ 700 GC-MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('ITQ 900', 'MS:1000636',
                    ('Thermo Scientific ITQ 900 GC-MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('ITQ 1100', 'MS:1000637',
                    ('Thermo Scientific ITQ 1100 GC-MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ XL ETD', 'MS:1000638',
                    ('Thermo Scientific LTQ XL MS with ETD.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap XL ETD', 'MS:1000639',
                    ('Thermo Scientific LTQ Orbitrap XL MS with ETD.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('DFS', 'MS:1000640',
                    ('Thermo Scientific DFS HR GC-MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('DSQ II', 'MS:1000641',
                    ('Thermo Scientific DSQ II GC-MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('MALDI LTQ XL', 'MS:1000642',
                    ('Thermo Scientific MALDI LTQ XL MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('MALDI LTQ Orbitrap', 'MS:1000643',
                    ('Thermo Scientific MALDI LTQ Orbitrap MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Quantum Access', 'MS:1000644',
                    ('Thermo Scientific TSQ Quantum Access MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Element XR', 'MS:1000645',
                    ('Thermo Scientific Element XR HR-ICP-MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Element 2', 'MS:1000646',
                    ('Thermo Scientific Element 2 HR-ICP-MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Element GD', 'MS:1000647',
                    ('Thermo Scientific Element GD Glow Discharge MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('GC IsoLink', 'MS:1000648',
                    ('Thermo Scientific GC IsoLink Isotope Ratio MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Exactive', 'MS:1000649',
                    ('Thermo Scientific Exactive MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Quantum Ultra AM', 'MS:1000743',
                    ('Thermo Scientific TSQ Quantum Ultra AM.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Quantum Ultra', 'MS:1000751',
                    ('Thermo Scientific TSQ Quantum Ultra.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ XL', 'MS:1000854',
                    ('Thermo Scientific LTQ XL MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Velos', 'MS:1000855',
                    ('Thermo Scientific LTQ Velos MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Velos ETD', 'MS:1000856',
                    ('Thermo Scientific LTQ Velos MS with ETD.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Vantage', 'MS:1001510',
                    ('TSQ Vantage.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap Velos', 'MS:1001742',
                    ('Finnigan LTQ Orbitrap Velos MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('ISQ', 'MS:1001908',
                    ('Thermo Scientific ISQ single quadrupole MS with the '
                     'ExtractraBrite source.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Velos Plus', 'MS:1001909',
                    ('Thermo Scientific second generation Velos.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap Elite', 'MS:1001910',
                    ('Thermo Scientific LTQ Orbitrap Elite, often just referred to '
                     'as the Orbitrap Elite.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Q Exactive', 'MS:1001911',
                    ('Thermo Scientific Q Exactive.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap Fusion', 'MS:1002416',
                    ('Thermo Scientific Orbitrap Fusion.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap Fusion ETD', 'MS:1002417',
                    ('Thermo Scientific Orbitrap Fusion with ETD.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Quantiva', 'MS:1002418',
                    ('Thermo Scientific TSQ Quantiva MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Endura', 'MS:1002419',
                    ('Thermo Scientific TSQ Endura MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Q Exactive HF', 'MS:1002523',
                    ('Thermo Scientific Q Exactive.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ 8000 Evo', 'MS:1002525',
                    ('Thermo Scientific TSQ 8000 Evo MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Exactive Plus', 'MS:1002526',
                    ('Thermo Scientific Exactive Plus MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Q Exactive Plus', 'MS:1002634',
                    ('Thermo Scientific Q Exactive Plus.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap Fusion Lumos', 'MS:1002732',
                    ('Thermo Scientific Orbitrap Fusion Lumos mass spectrometer '
                     'with Tribrid architecture consisting of quadrupole mass '
                     'filter, linear ion trap and Orbitrap mass analyzers.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap Classic', 'MS:1002835',
                    ('Thermo Fisher Scientific LTQ Orbitrap Classic.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Altis', 'MS:1002874',
                    ('Thermo Scientific TSQ Altis Triple Quadrupole MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ Quantis', 'MS:1002875',
                    ('Thermo Scientific TSQ Quantis Triple Quadrupole MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('TSQ 9000', 'MS:1002876',
                    ('Thermo Scientific TSQ 9000 Triple Quadrupole MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Q Exactive HF-X', 'MS:1002877',
                    ('Thermo Scientific Q Exactive HF-X Hybrid Quadrupole Orbitrap '
                     'MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap Exploris 480', 'MS:1003028',
                    ('Thermo Scientific Orbitrap Exploris 480 Quadrupole Orbitrap '
                     'MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap Eclipse', 'MS:1003029',
                    ('Thermo Scientific Orbitrap Eclipse mass spectrometer with '
                     'Tribrid architecture consisting of quadrupole mass filter, '
                     'linear ion trap and Orbitrap mass analyzers.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap Exploris 240', 'MS:1003094',
                    ('Thermo Scientific Orbitrap Exploris 240 Quadrupole Orbitrap '
                     'MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap Exploris 120', 'MS:1003095',
                    ('Thermo Scientific Orbitrap Exploris 120 Quadrupole Orbitrap '
                     'MS.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('LTQ Orbitrap Velos Pro', 'MS:1003096',
                    ('Thermo Scientific LTQ Orbitrap Velos Pro, often just '
                     'referred to as the Orbitrap Velos Pro.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('Orbitrap ID-X', 'MS:1003112',
                    ('Thermo Scientific Orbitrap ID-X mass spectrometer with '
                     'Tribrid architecture consisting of quadrupole mass filter, '
                     'linear ion trap and Orbitrap mass analyzers.'),
                    'instrument model',
                    ['Thermo Scientific instrument model', 'Thermo Fisher Scientific instrument model', 'instrument model']),
    InstrumentModel('explorer', 'MS:1000158',
                    ('IonSpec Explorer MS.'),
                    'instrument model',
                    ['IonSpec instrument model', 'Varian instrument model', 'instrument model']),
    InstrumentModel('HiRes ESI', 'MS:1000162',
                    ('IonSpec HiResESI MS.'),
                    'instrument model',
                    ['IonSpec instrument model', 'Varian instrument model', 'instrument model']),
    InstrumentModel('HiRes MALDI', 'MS:1000163',
                    ('IonSpec HiResMALDI MS.'),
                    'instrument model',
                    ['IonSpec instrument model', 'Varian instrument model', 'instrument model']),
    InstrumentModel('OMEGA', 'MS:1000181',
                    ('IonSpec OMEGA MS.'),
                    'instrument model',
                    ['IonSpec instrument model', 'Varian instrument model', 'instrument model']),
    InstrumentModel('OMEGA-2001', 'MS:1000182',
                    ('IonSpec OMEGA-2001 MS.'),
                    'instrument model',
                    ['IonSpec instrument model', 'Varian instrument model', 'instrument model']),
    InstrumentModel('ultima', 'MS:1000200',
                    ('IonSpec Ultima MS.'),
                    'instrument model',
                    ['IonSpec instrument model', 'Varian instrument model', 'instrument model']),
])
# [[[end]]]


class ComponentGroup(object):
    """Represent a collection of :class:`Component` objects which are
    part of a single instrument subsystem.

    Attributes
    ----------
    type: str
        The name of the instrument subsystem
    members: list
        The components contained
    order: int
        The position of this group in the sequence of steps
        that analytes take when passing through the instrument.
    """

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
    '''Describes a single configuration of an instrument for acquiring spectra.

    Parameters
    ----------
    id: object
        The within-run unique identifier of this configuration. May be an integer or string
    groups: :class:`list` of :class:`ComponentGroup`
        The component groups, sorted by :attr:`ComponentGroup.order`, in this configuraton
    analyzers: list
        A convenience list for storing the :class:`Component` objects which are ``analyzer``s
    '''

    def __init__(self, id=None, groups=None, model=None, serial_number=None, name=None):
        self.id = id or uid()
        self.groups = sorted(groups or [], key=lambda x: x.order)
        self.analyzers = []

        self.model = model
        self.serial_number = serial_number
        self.name = name

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


__all__ = [
    "InstrumentInformation", "ComponentGroup", "Component",
    "all_components_by_name", "ionization_types", "detector_types",
    "analyzer_types", "inlet_types", "component"
]
