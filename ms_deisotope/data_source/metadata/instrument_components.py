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
ionization_types = TermSet([
    Component(u'chemical ionization', u'MS:1000071',
              (u'The formation of a new ion by the reaction of a neutral'
               u'species with an ion. The process may involve transfer of an'
               u'electron, a proton or other charged species between the'
               u'reactants. When a positive ion results from chemical'
               u'ionization the term may be used without qualification. When'
               u'a negative ion results the term negative ion chemical'
               u'ionization should be used. Note that this term is not'
               u'synonymous with chemi-ionization.'),
              'ionization type',
              [u'ionization type']),
    Component(u'electrospray ionization', u'MS:1000073',
              (u'A process in which ionized species in the gas phase are'
               u'produced from an analyte-containing solution via highly'
               u'charged fine droplets, by means of spraying the solution'
               u'from a narrow-bore needle tip at atmospheric pressure in the'
               u'presence of a high electric field. When a pressurized gas is'
               u'used to aid in the formation of a stable spray, the term'
               u'pneumatically assisted electrospray ionization is used. The'
               u'term ion spray is not recommended.'),
              'ionization type',
              [u'ionization type']),
    Component(u'fast atom bombardment ionization', u'MS:1000074',
              (u'The ionization of any species by the interaction of a'
               u'focused beam of neutral atoms having a translational energy'
               u'of several thousand eV with a sample that is typically'
               u'dissolved in a solvent matrix. See also secondary'
               u'ionization.'),
              'ionization type',
              [u'ionization type']),
    Component(u'multiphoton ionization', u'MS:1000227',
              (u'Photoionization of an atom or molecule in which in two or'
               u'more photons are absorbed.'),
              'ionization type',
              [u'ionization type']),
    Component(u'atmospheric pressure ionization', u'MS:1000240',
              (u'Any ionization process in which ions are formed in the gas'
               u'phase at atmospheric pressure.'),
              'ionization type',
              [u'ionization type']),
    Component(u'desorption ionization', u'MS:1000247',
              (u'The formation of ions from a solid or liquid material after'
               u'the rapid vaporization of that sample.'),
              'ionization type',
              [u'ionization type']),
    Component(u'flowing afterglow', u'MS:1000255',
              (u'An ion source immersed in a flow of helium or other inert'
               u'buffer gas that carries the ions through a meter-long'
               u'reactor at pressures around 100 Pa.'),
              'ionization type',
              [u'ionization type']),
    Component(u'field ionization', u'MS:1000258',
              (u'The removal of electrons from any species by interaction'
               u'with a high electric field.'),
              'ionization type',
              [u'ionization type']),
    Component(u'glow discharge ionization', u'MS:1000259',
              (u'The formation of ions in the gas phase and from solid'
               u'samples at the cathode by application of a voltage to a low'
               u'pressure gas.'),
              'ionization type',
              [u'ionization type']),
    Component(u'Negative Ion chemical ionization', u'MS:1000271',
              (u'Chemical ionization that results in the formation of'
               u'negative ions.'),
              'ionization type',
              [u'ionization type']),
    Component(u'neutralization reionization mass spectrometry', u'MS:1000272',
              (u'With this technique, m/z selected ions form neutrals by'
               u'charge transfer to a collision gas or by dissociation. The'
               u'neutrals are separated from the remaining ions and ionized'
               u'in collisions with a second gas. This method is used to'
               u'investigate reaction intermediates and other unstable'
               u'species.'),
              'ionization type',
              [u'ionization type']),
    Component(u'photoionization', u'MS:1000273',
              (u'The ionization of an atom or molecule by a photon, written M'
               u'+ h? ? M^+ + e. The term photon impact is not recommended.'),
              'ionization type',
              [u'ionization type']),
    Component(u'pyrolysis mass spectrometry', u'MS:1000274',
              (u'A mass spectrometry technique in which the sample is heated'
               u'to the point of decomposition and the gaseous decomposition'
               u'products are introduced into the ion source.'),
              'ionization type',
              [u'ionization type']),
    Component(u'resonance enhanced multiphoton ionization', u'MS:1000276',
              (u'Multiphoton ionization in which the ionization cross section'
               u'is significantly enhanced because the energy of the incident'
               u'photons is resonant with an intermediate excited state of'
               u'the neutral species.'),
              'ionization type',
              [u'ionization type']),
    Component(u'adiabatic ionization', u'MS:1000380',
              (u'A process whereby an electron is removed from an atom, ion,'
               u'or molecule to produce an ion in its lowest energy state.'),
              'ionization type',
              [u'ionization type']),
    Component(u'associative ionization', u'MS:1000381',
              (u'An ionization process in which two excited atoms or'
               u'molecules react to form a single positive ion and an'
               u'electron.'),
              'ionization type',
              [u'ionization type']),
    Component(u'autodetachment', u'MS:1000383',
              (u'The formation of a neutral when a negative ion in a discrete'
               u'state with an energy greater than the detachment threshold'
               u'loses an electron spontaneously without further interaction'
               u'with an energy source.'),
              'ionization type',
              [u'ionization type']),
    Component(u'autoionization', u'MS:1000384',
              (u'The formation of an ion when an atom or molecule in a'
               u'discrete state with an energy greater than the ionization'
               u'threshold loses an electron spontaneously without further'
               u'interaction with an energy source.'),
              'ionization type',
              [u'ionization type']),
    Component(u'charge exchange ionization', u'MS:1000385',
              (u'The interaction of an ion with an atom or molecule in which'
               u'the charge on the ion is transferred to the neutral without'
               u'the dissociation of either. Synonymous with charge transfer'
               u'ionization.'),
              'ionization type',
              [u'ionization type']),
    Component(u'chemi-ionization', u'MS:1000386',
              (u'The reaction of a neutral molecule with an internally'
               u'excited molecule to form an ion. Note that this term is not'
               u'synonymous with chemical ionization.'),
              'ionization type',
              [u'ionization type']),
    Component(u'dissociative ionization', u'MS:1000388',
              (u'The reaction of a gas-phase molecule that results in its'
               u'decomposition to form products, one of which is an ion.'),
              'ionization type',
              [u'ionization type']),
    Component(u'electron ionization', u'MS:1000389',
              (u'The ionization of an atom or molecule by electrons that are'
               u'typically accelerated to energies between 50 and 150 eV.'
               u'Usually 70 eV electrons are used to produce positive ions.'
               u"The term 'electron impact' is not recommended."),
              'ionization type',
              [u'ionization type']),
    Component(u'liquid secondary ionization', u'MS:1000395',
              (u'The ionization of any species by the interaction of a'
               u'focused beam of ions with a sample that is dissolved in a'
               u'solvent matrix. See also fast atom bombardment and secondary'
               u'ionization.'),
              'ionization type',
              [u'ionization type']),
    Component(u'penning ionization', u'MS:1000399',
              (u'Ionization that occurs through the interaction of two or'
               u'more neutral gaseous species, at least one of which is'
               u'internally excited.'),
              'ionization type',
              [u'ionization type']),
    Component(u'plasma desorption ionization', u'MS:1000400',
              (u'The ionization of material in a solid sample by bombarding'
               u'it with ionic or neutral atoms formed as a result of the'
               u'fission of a suitable nuclide, typically 252Cf. Synonymous'
               u'with fission fragment ionization.'),
              'ionization type',
              [u'ionization type']),
    Component(u'secondary ionization', u'MS:1000402',
              (u'The process in which ions are ejected from a sample surface'
               u'as a result of bombardment by a primary beam of atoms or'
               u'ions.'),
              'ionization type',
              [u'ionization type']),
    Component(u'soft ionization', u'MS:1000403',
              (u'The formation of gas-phase ions without extensive'
               u'fragmentation.'),
              'ionization type',
              [u'ionization type']),
    Component(u'spark ionization', u'MS:1000404',
              (u'The formation of ions from a solid material by an'
               u'intermittent electrical discharge.'),
              'ionization type',
              [u'ionization type']),
    Component(u'surface ionization', u'MS:1000406',
              (u'The ionization of a neutral species when it interacts with a'
               u'solid surface with an appropriate work function and'
               u'temperature.'),
              'ionization type',
              [u'ionization type']),
    Component(u'thermal ionization', u'MS:1000407',
              (u'The ionization of a neutral species through contact with a'
               u'high temperature surface.'),
              'ionization type',
              [u'ionization type']),
    Component(u'vertical ionization', u'MS:1000408',
              (u'A process in which an electron is removed from or added to a'
               u'molecule without a change in the positions of the atoms. The'
               u'resulting ion is typically in an excited vibrational state.'),
              'ionization type',
              [u'ionization type']),
    Component(u'fast ion bombardment', u'MS:1000446',
              (u'The ionization of any species by the interaction of a'
               u'focused beam of ions having a translational energy of'
               u'several thousand eV with a solid sample.'),
              'ionization type',
              [u'ionization type']),
    Component(u'microelectrospray', u'MS:1000397',
              (u'Electrospray ionization at a solvent flow rate of 300-800'
               u'nL/min where the flow is a result of a mechanical pump. See'
               u'nanoelectrospray.'),
              'ionization type',
              [u'electrospray ionization', u'ionization type']),
    Component(u'nanoelectrospray', u'MS:1000398',
              (u'Electrospray ionization at a flow rate less than ~25 nL/min.'
               u'Nanoelectrospray is synonymous with nanospray. The flow is'
               u'dependent on the potenial on the tip of the electrospray'
               u'needle and/or a gas presure to push the sample through the'
               u'needle. See also electrospray ionization and'
               u'microelectrospray.'),
              'ionization type',
              [u'electrospray ionization', u'ionization type']),
    Component(u'atmospheric pressure chemical ionization', u'MS:1000070',
              (u'Chemical ionization that takes place at atmospheric pressure'
               u'as opposed to the reduced pressure is normally used for'
               u'chemical ionization.'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'atmospheric pressure matrix-assisted laser desorption ionization', u'MS:1000239',
              (u'Matrix-assisted laser desorption ionization in which the'
               u'sample target is at atmospheric pressure and the ions formed'
               u'by the pulsed laser are sampled through a small aperture'
               u'into the mass spectrometer.'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'atmospheric pressure photoionization', u'MS:1000382',
              (u'Atmospheric pressure chemical ionization in which the'
               u'reactant ions are generated by photo-ionization.'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'desorption electrospray ionization', u'MS:1002011',
              (u'Combination of electrospray and desorption ionization method'
               u'that ionizes gases, liquids and solids in open air under'
               u'atmospheric pressure." [DOI:10.1126/science.1104404'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'matrix-assisted laser desorption ionization', u'MS:1000075',
              (u'The formation of gas-phase ions from molecules that are'
               u'present in a solid or solvent matrix that is irradiated with'
               u'a pulsed laser. See also laser desorption/ionization.'),
              'ionization type',
              [u'desorption ionization', u'ionization type']),
    Component(u'field desorption', u'MS:1000257',
              (u'The formation of gas-phase ions from a material deposited on'
               u'a solid surface in the presence of a high electric field.'
               u'Because this process may encompass ionization by field'
               u'ionization or other mechanisms, it is not recommended as a'
               u'synonym for field desorption ionization.'),
              'ionization type',
              [u'desorption ionization', u'ionization type']),
    Component(u'desorption/ionization on silicon', u'MS:1000387',
              (u'The formation of ions by laser desorption ionization of a'
               u'sample deposited on a porous silicon surface.'),
              'ionization type',
              [u'desorption ionization', u'ionization type']),
    Component(u'laser desorption ionization', u'MS:1000393',
              (u'The formation of gas-phase ions by the interaction of a'
               u'pulsed laser with a solid or liquid material.'),
              'ionization type',
              [u'desorption ionization', u'ionization type']),
    Component(u'surface-assisted laser desorption ionization', u'MS:1000405',
              (u'The formation of gas-phase ions from molecules that are'
               u'deposited on a particular surface substrate that is'
               u'irradiated with a pulsed laser. See also matrix-assisted'
               u'laser desorption ionization.'),
              'ionization type',
              [u'desorption ionization', u'ionization type']),
    Component(u'surface enhanced laser desorption ionization', u'MS:1000278',
              (u'The formation of ionized species in the gas phase from'
               u'analytes deposited on a particular surface substrate which'
               u'is irradiated with a laser beam of which wavelength is'
               u'absorbed by the surface. See also desorption/ionization on'
               u'silicon and laser desorption/ionization.'),
              'ionization type',
              [u'surface ionization', u'ionization type']),
    Component(u'surface enhanced neat desorption', u'MS:1000279',
              (u'Matrix-assisted laser desorption ionization in which the'
               u'matrix is covalently linked to the target surface.'),
              'ionization type',
              [u'surface ionization', u'ionization type']),
])
# [[[end]]]


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('detector type', term_cls_name="Component", writer=cog.out)
# ]]]
detector_types = TermSet([
    Component(u'channeltron', u'MS:1000107',
              (u'A horn-shaped (or cone-shaped) continuous dynode particle'
               u'multiplier. The ion strikes the inner surface of the device'
               u'and induces the production of secondary electrons that in'
               u'turn impinge on the inner surfaces to produce more secondary'
               u'electrons. This avalanche effect produces an increase in'
               u'signal in the final measured current pulse.'),
              'detector type',
              [u'detector type']),
    Component(u'daly detector', u'MS:1000110',
              (u'Detector consisting of a conversion dynode, scintillator and'
               u'photomultiplier. The metal knob at high potential emits'
               u'secondary electrons when ions impinge on the surface. The'
               u'secondary electrons are accelerated onto the scintillator'
               u'that produces light that is then detected by the'
               u'photomultiplier detector.'),
              'detector type',
              [u'detector type']),
    Component(u'faraday cup', u'MS:1000112',
              (u'A conducting cup or chamber that intercepts a charged'
               u'particle beam and is electrically connected to a current'
               u'measuring device.'),
              'detector type',
              [u'detector type']),
    Component(u'multi-collector', u'MS:1000115',
              (u'A detector system commonly used in inductively coupled'
               u'plasma mass spectrometers.'),
              'detector type',
              [u'detector type']),
    Component(u'photomultiplier', u'MS:1000116',
              (u'A detector for conversion of the ion/electron signal into'
               u'photon(s) which are then amplified and detected.'),
              'detector type',
              [u'detector type']),
    Component(u'electron multiplier', u'MS:1000253',
              (u'A device to amplify the current of a beam or packet of'
               u'charged particles or photons by incidence upon the surface'
               u'of an electrode to produce secondary electrons. The'
               u'secondary electrons are then accelerated to other electrodes'
               u'or parts of a continuous electrode to produce further'
               u'secondary electrons.'),
              'detector type',
              [u'detector type']),
    Component(u'array detector', u'MS:1000345',
              (u'Detector comprising several ion collection elements,'
               u'arranged in a line or grid where each element is an'
               u'individual detector.'),
              'detector type',
              [u'detector type']),
    Component(u'conversion dynode', u'MS:1000346',
              (u'A surface that is held at high potential such that ions'
               u'striking the surface produce electrons that are subsequently'
               u'detected.'),
              'detector type',
              [u'detector type']),
    Component(u'dynode', u'MS:1000347',
              (u'One of a series of electrodes in a photomultiplier tube.'
               u'Such an arrangement is able to amplify the current emitted'
               u'by the photocathode.'),
              'detector type',
              [u'detector type']),
    Component(u'focal plane collector', u'MS:1000348',
              (u'A detector for spatially disperse ion beams in which all'
               u'ions simultaneously impinge on the detector plane.'),
              'detector type',
              [u'detector type']),
    Component(u'ion-to-photon detector', u'MS:1000349',
              (u'A detector in which ions strike a conversion dynode to'
               u'produce electrons that in turn strike a phosphor and the'
               u'resulting photons are detected by a photomultiplier.'),
              'detector type',
              [u'detector type']),
    Component(u'point collector', u'MS:1000350',
              (u'A detector in which the ion beam is focused onto a point and'
               u'the individual ions arrive sequentially.'),
              'detector type',
              [u'detector type']),
    Component(u'postacceleration detector', u'MS:1000351',
              (u'A detector in which the charged particles are accelerated to'
               u'a high velocity and impinge on a conversion dynode, emitting'
               u'secondary electrons. The electrons are accelerated onto a'
               u'phosphor screen, which emits photons that are in turn'
               u'detected using a photomultiplier or other photon detector.'),
              'detector type',
              [u'detector type']),
    Component(u'inductive detector', u'MS:1000624',
              (u'Inductive detector.'),
              'detector type',
              [u'detector type']),
    Component(u'fluorescence detector', u'MS:1002308',
              (u'A detector using a fluorescent signal after excitation with'
               u'light.'),
              'detector type',
              [u'detector type']),
    Component(u'electron multiplier tube', u'MS:1000111',
              (u'A device to amplify the current of a beam or packet of'
               u'charged particles or photons by incidence upon the surface'
               u'of an electrode to produce secondary electrons.'),
              'detector type',
              [u'electron multiplier', u'detector type']),
    Component(u'microchannel plate detector', u'MS:1000114',
              (u'A thin plate that contains a closely spaced array of'
               u'channels that each act as a continuous dynode particle'
               u'multiplier. A charged particle, fast neutral particle, or'
               u'photon striking the plate causes a cascade of secondary'
               u'electrons that ultimately exits the opposite side of the'
               u'plate.'),
              'detector type',
              [u'array detector', u'detector type']),
    Component(u'photodiode array detector', u'MS:1000621',
              (u'An array detector used to record spectra in the ultraviolet'
               u'and visible region of light.'),
              'detector type',
              [u'array detector', u'detector type']),
    Component(u'conversion dynode electron multiplier', u'MS:1000108',
              (u'A surface that is held at high potential so that ions'
               u'striking the surface produce electrons that are subsequently'
               u'detected.'),
              'detector type',
              [u'conversion dynode', u'detector type']),
    Component(u'conversion dynode photomultiplier', u'MS:1000109',
              (u'A detector in which ions strike a conversion dynode to'
               u'produce electrons that in turn generate photons through a'
               u'phosphorescent screen that are detected by a'
               u'photomultiplier.'),
              'detector type',
              [u'conversion dynode', u'detector type']),
    Component(u'focal plane array', u'MS:1000113',
              (u'An array of detectors for spatially disperse ion beams in'
               u'which all ions simultaneously impinge on the detector plane.'),
              'detector type',
              [u'focal plane collector', u'detector type']),
    Component(u'Acquity UPLC FLR', u'MS:1000819',
              (u'Acquity UPLC Fluorescence Detector.'),
              'detector type',
              [u'Waters instrument model', u'fluorescence detector', u'instrument model', u'detector type']),
    Component(u'Acquity UPLC PDA', u'MS:1000818',
              (u'Acquity UPLC Photodiode Array Detector.'),
              'detector type',
              [u'Waters instrument model', u'photodiode array detector', u'instrument model', u'array detector', u'detector type']),
])
# [[[end]]]


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('mass analyzer type', 'analyzer_types', term_cls_name="Component", writer=cog.out)
# ]]]
analyzer_types = TermSet([
    Component(u'fourier transform ion cyclotron resonance mass spectrometer', u'MS:1000079',
              (u'A mass spectrometer based on the principle of ion cyclotron'
               u'resonance in which an ion in a magnetic field moves in a'
               u'circular orbit at a frequency characteristic of its m/z'
               u'value. Ions are coherently excited to a larger radius orbit'
               u'using a pulse of radio frequency energy and their image'
               u'charge is detected on receiver plates as a time domain'
               u'signal. Fourier transformation of the time domain signal'
               u'results in a frequency domain signal which is converted to a'
               u'mass spectrum based in the inverse relationship between'
               u'frequency and m/z.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'magnetic sector', u'MS:1000080',
              (u'A device that produces a magnetic field perpendicular to a'
               u'charged particle beam that deflects the beam to an extent'
               u'that is proportional to the particle momentum per unit'
               u'charge. For a monoenergetic beam, the deflection is'
               u'proportional to m/z.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'quadrupole', u'MS:1000081',
              (u'A mass spectrometer that consists of four parallel rods'
               u'whose centers form the corners of a square and whose'
               u'opposing poles are connected. The voltage applied to the'
               u'rods is a superposition of a static potential and a'
               u'sinusoidal radio frequency potential. The motion of an ion'
               u'in the x and y dimensions is described by the Matthieu'
               u'equation whose solutions show that ions in a particular m/z'
               u'range can be transmitted along the z axis.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'time-of-flight', u'MS:1000084',
              (u'Instrument that separates ions by m/z in a field-free region'
               u'after acceleration to a fixed acceleration energy.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'electrostatic energy analyzer', u'MS:1000254',
              (u'A device consisting of conducting parallel plates,'
               u'concentric cylinders or concentric spheres that separates'
               u'charged particles according to their kinetic energy by means'
               u'of an electric field that is constant in time.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'ion trap', u'MS:1000264',
              (u'A device for spatially confining ions using electric and'
               u'magnetic fields alone or in combination.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'stored waveform inverse fourier transform', u'MS:1000284',
              (u'A technique to create excitation waveforms for ions in FT-'
               u'ICR mass spectrometer or Paul ion trap. An excitation'
               u'waveform in the time-domain is generated by taking the'
               u'inverse Fourier transform of an appropriate frequency-domain'
               u'programmed excitation spectrum, in which the resonance'
               u'frequencies of ions to be excited are included. This'
               u'technique may be used for selection of precursor ions in MS2'
               u'experiments.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'cyclotron', u'MS:1000288',
              (u'A device that uses an oscillating electric field and'
               u'magnetic field to accelerate charged particles.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'orbitrap', u'MS:1000484',
              (u'An ion trapping device that consists of an outer barrel-like'
               u'electrode and a coaxial inner spindle-like electrode that'
               u'form an electrostatic field with quadro-logarithmic'
               u'potential distribution. The frequency of harmonic'
               u'oscillations of the orbitally trapped ions along the axis of'
               u'the electrostatic field is independent of the ion velocity'
               u'and is inversely proportional to the square root of m/z so'
               u'that the trap can be used as a mass analyzer.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'quadrupole ion trap', u'MS:1000082',
              (u'Quadrupole Ion Trap mass analyzer captures the ions in a'
               u'three dimensional ion trap and then selectively ejects them'
               u'by varying the RF and DC potentials.'),
              'mass analyzer type',
              [u'ion trap', u'mass analyzer type']),
    Component(u'linear ion trap', u'MS:1000291',
              (u'A two dimensional Paul ion trap in which ions are confined'
               u'in the axial dimension by means of an electric field at the'
               u'ends of the trap.'),
              'mass analyzer type',
              [u'ion trap', u'mass analyzer type']),
    Component(u'axial ejection linear ion trap', u'MS:1000078',
              (u'A linear ion trap mass spectrometer where ions are ejected'
               u'along the axis of the analyzer.'),
              'mass analyzer type',
              [u'linear ion trap', u'ion trap', u'mass analyzer type']),
    Component(u'radial ejection linear ion trap', u'MS:1000083',
              (u'A linear ion trap mass spectrometer where ions are ejected'
               u'along the radius of the analyzer.'),
              'mass analyzer type',
              [u'linear ion trap', u'ion trap', u'mass analyzer type']),
])
# [[[end]]]


# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('inlet type', term_cls_name="Component", writer=cog.out)
# ]]]
inlet_types = TermSet([
    Component(u'continuous flow fast atom bombardment', u'MS:1000055',
              (u'Fast atom bombardment ionization in which the analyte in'
               u'solution is entrained in a flowing liquid matrix.'),
              'inlet type',
              [u'inlet type']),
    Component(u'direct inlet', u'MS:1000056',
              (u'The sample is directly inserted into the ion source, usually'
               u'on the end of a heatable probe.'),
              'inlet type',
              [u'inlet type']),
    Component(u'electrospray inlet', u'MS:1000057',
              (u'Inlet used for introducing the liquid sample into an'
               u'electrospray ionization source.'),
              'inlet type',
              [u'inlet type']),
    Component(u'flow injection analysis', u'MS:1000058',
              (u'Sample is directly injected or infused into the ionization'
               u'source.'),
              'inlet type',
              [u'inlet type']),
    Component(u'inductively coupled plasma', u'MS:1000059',
              (u'A gas discharge ion source in which the energy to the plasma'
               u'is supplied by electromagnetic induction.'),
              'inlet type',
              [u'inlet type']),
    Component(u'infusion', u'MS:1000060',
              (u'The continuous flow of solution of a sample into the'
               u'ionization source.'),
              'inlet type',
              [u'inlet type']),
    Component(u'jet separator', u'MS:1000061',
              (u'A device that separates carrier gas from gaseous analyte'
               u'molecules on the basis of diffusivity.'),
              'inlet type',
              [u'inlet type']),
    Component(u'membrane separator', u'MS:1000062',
              (u'A device to separate carrier molecules from analyte'
               u'molecules on the basis of ease of diffusion across a'
               u'semipermeable membrane.'),
              'inlet type',
              [u'inlet type']),
    Component(u'moving belt', u'MS:1000063',
              (u'Continuous moving surface in the form of a belt which passes'
               u'through an ion source carrying analyte molecules.'),
              'inlet type',
              [u'inlet type']),
    Component(u'moving wire', u'MS:1000064',
              (u'Continuous moving surface in the form of a wire which passes'
               u'through an ion source carrying analyte molecules.'),
              'inlet type',
              [u'inlet type']),
    Component(u'open split', u'MS:1000065',
              (u'A division of flowing stream of liquid into two streams.'),
              'inlet type',
              [u'inlet type']),
    Component(u'particle beam', u'MS:1000066',
              (u'Method for generating ions from a solution of an analyte.'),
              'inlet type',
              [u'inlet type']),
    Component(u'reservoir', u'MS:1000067',
              (u'A sample inlet method involving a reservoir.'),
              'inlet type',
              [u'inlet type']),
    Component(u'septum', u'MS:1000068',
              (u'A disc composed of a flexible material that seals the'
               u'entrance to the reservoir. Can also be entrance to the'
               u'vacuum chamber.'),
              'inlet type',
              [u'inlet type']),
    Component(u'thermospray inlet', u'MS:1000069',
              (u'A method for generating gas phase ions from a solution of an'
               u'analyte by rapid heating of the sample.'),
              'inlet type',
              [u'inlet type']),
    Component(u'direct insertion probe', u'MS:1000248',
              (u'A device for introducing a solid or liquid sample into a'
               u'mass spectrometer ion source for desorption ionization.'),
              'inlet type',
              [u'inlet type']),
    Component(u'direct liquid introduction', u'MS:1000249',
              (u'The delivery of a liquid sample into a mass spectrometer for'
               u'spray or desorption ionization.'),
              'inlet type',
              [u'inlet type']),
    Component(u'membrane inlet', u'MS:1000396',
              (u'A semi-permeable membrane separator that permits the passage'
               u'of gas sample directly to the mass spectrometer ion source.'),
              'inlet type',
              [u'inlet type']),
    Component(u'nanospray inlet', u'MS:1000485',
              (u'Nanospray Inlet.'),
              'inlet type',
              [u'electrospray inlet', u'inlet type']),
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
instrument_models = TermSet([
    InstrumentModel(u'SCIEX instrument model', u'MS:1000121',
                    (u'The brand of instruments from the joint venture between'
                     u'Applied Biosystems and MDS Analytical Technologies (formerly'
                     u'MDS SCIEX). Previously branded as \\"Applied Biosystems|MDS'
                     u'SCIEX\\".'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Bruker Daltonics instrument model', u'MS:1000122',
                    (u"Bruker Daltonics' instrument model."),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Shimadzu instrument model', u'MS:1000124',
                    (u'Shimadzu corporation instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Waters instrument model', u'MS:1000126',
                    (u'Waters Corporation instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Thermo Fisher Scientific instrument model', u'MS:1000483',
                    (u'Thermo Fisher Scientific instrument model. The company has'
                     u'gone through several names including Thermo Finnigan, Thermo'
                     u'Scientific.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Hitachi instrument model', u'MS:1000488',
                    (u'Hitachi instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Varian instrument model', u'MS:1000489',
                    (u'Varian instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Agilent instrument model', u'MS:1000490',
                    (u'Agilent instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Dionex instrument model', u'MS:1000491',
                    (u'Dionex instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'Applied Biosystems instrument model', u'MS:1000495',
                    (u'Applied Biosystems instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'LECO instrument model', u'MS:1001800',
                    (u'LECO instrument model.'),
                    'instrument model',
                    [u'instrument model']),
    InstrumentModel(u'4000 QTRAP', u'MS:1000139',
                    (u'Applied Biosystems/MDS SCIEX Q 4000 TRAP MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 150EX', u'MS:1000143',
                    (u'Applied Biosystems/MDS SCIEX API 150EX MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 150EX Prep', u'MS:1000144',
                    (u'Applied Biosystems/MDS SCIEX API 150EX Prep MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 2000', u'MS:1000145',
                    (u'Applied Biosystems/MDS SCIEX API 2000 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 3000', u'MS:1000146',
                    (u'Applied Biosystems/MDS SCIEX API 3000 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 4000', u'MS:1000147',
                    (u'Applied Biosystems/MDS SCIEX API 4000 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'proteomics solution 1', u'MS:1000186',
                    (u'Applied Biosystems/MDS SCIEX Proteomics Solution 1 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Q TRAP', u'MS:1000187',
                    (u'Applied Biosystems/MDS SCIEX Q TRAP MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QSTAR', u'MS:1000190',
                    (u'Applied Biosystems/MDS SCIEX QSTAR MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'SymBiot I', u'MS:1000194',
                    (u'Applied Biosystems/MDS SCIEX SymBiot I MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'SymBiot XVI', u'MS:1000195',
                    (u'Applied Biosystems/MDS SCIEX SymBiot XVI MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'3200 QTRAP', u'MS:1000651',
                    (u'SCIEX or Applied Biosystems|MDS SCIEX QTRAP 3200.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'4800 Plus MALDI TOF/TOF', u'MS:1000652',
                    (u'SCIEX or Applied Biosystems|MDS SCIEX 4800 Plus MALDI TOF-'
                     u'TOF Analyzer.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 3200', u'MS:1000653',
                    (u'SCIEX or Applied Biosystems|MDS SCIEX API 3200 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 5000', u'MS:1000654',
                    (u'SCIEX or Applied Biosystems|MDS SCIEX API 5000 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QSTAR Elite', u'MS:1000655',
                    (u'SCIEX or Applied Biosystems|MDS SCIEX QSTAR Elite.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QSTAR Pulsar', u'MS:1000656',
                    (u'Applied Biosystems|MDS SCIEX QSTAR Pulsar.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QSTAR XL', u'MS:1000657',
                    (u'Applied Biosystems|MDS SCIEX QSTAR XL.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QTRAP 5500', u'MS:1000931',
                    (u'Applied Biosystems|MDS SCIEX QTRAP 5500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'TripleTOF 5600', u'MS:1000932',
                    (u'SCIEX TripleTOF 5600, a quadrupole - quadrupole - time-of-'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'5800 TOF/TOF', u'MS:1001482',
                    (u'SCIEX 5800 TOF-TOF Analyzer.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'TripleTOF 6600', u'MS:1002533',
                    (u'SCIEX TripleTOF 6600, a quadrupole - quadrupole - time-of-'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'2000 QTRAP', u'MS:1002577',
                    (u'SCIEX 2000 QTRAP.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'2500 QTRAP', u'MS:1002578',
                    (u'SCIEX 2500 QTRAP.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'3500 QTRAP', u'MS:1002579',
                    (u'SCIEX 3500 QTRAP.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QTRAP 4500', u'MS:1002580',
                    (u'SCIEX QTRAP 4500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QTRAP 6500', u'MS:1002581',
                    (u'SCIEX QTRAP 6500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'QTRAP 6500+', u'MS:1002582',
                    (u'SCIEX QTRAP 6500+.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'TripleTOF 4600', u'MS:1002583',
                    (u'SCIEX TripleTOF 4600 time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'TripleTOF 5600+', u'MS:1002584',
                    (u'SCIEX TripleTOF 5600+ time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 100', u'MS:1002585',
                    (u'Applied Biosystems/MDS SCIEX API 100 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 100LC', u'MS:1002586',
                    (u'Applied Biosystems/MDS SCIEX API 100LC MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 165', u'MS:1002587',
                    (u'Applied Biosystems/MDS SCIEX API 165 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 300', u'MS:1002588',
                    (u'Applied Biosystems/MDS SCIEX API 300 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 350', u'MS:1002589',
                    (u'Applied Biosystems/MDS SCIEX API 350 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'API 365', u'MS:1002590',
                    (u'Applied Biosystems/MDS SCIEX API 365 MS.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Triple Quad 3500', u'MS:1002591',
                    (u'SCIEX Triple Quad 3500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Triple Quad 4500', u'MS:1002592',
                    (u'SCIEX Triple Quad 4500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Triple Quad 5500', u'MS:1002593',
                    (u'SCIEX Triple Quad 5500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Triple Quad 6500', u'MS:1002594',
                    (u'SCIEX Triple Quad 6500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Triple Quad 6500+', u'MS:1002595',
                    (u'SCIEX Triple Quad 6500+.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'X500R QTOF', u'MS:1002674',
                    (u'SCIEX X500R QTOF, a quadrupole - quadrupole - time-of-flight'
                     u'mass spectrometer.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Triple Quad 7500', u'MS:1003144',
                    (u'SCIEX Triple Quad 7500.'),
                    'instrument model',
                    [u'SCIEX instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics HCT Series', u'MS:1000697',
                    (u"Bruker Daltonics' HCT Series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics esquire series', u'MS:1001533',
                    (u"Bruker Daltonics' esquire series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics flex series', u'MS:1001534',
                    (u"Bruker Daltonics' flex series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics BioTOF series', u'MS:1001535',
                    (u"Bruker Daltonics' BioTOF series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics micrOTOF series', u'MS:1001536',
                    (u"Bruker Daltonics' micrOTOF series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics amaZon series', u'MS:1001545',
                    (u"Bruker Daltonics' amaZon series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics maXis series', u'MS:1001547',
                    (u"Bruker Daltonics' maXis series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics solarix series', u'MS:1001548',
                    (u"Bruker Daltonics' solarix: ESI quadrupole ion trap, APCI,"
                     u'APPI, ETD, PTR.'),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics apex series', u'MS:1001556',
                    (u"Bruker Daltonics' apex series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics SCION series', u'MS:1002293',
                    (u"Bruker Daltonics' SCION series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics EVOQ series', u'MS:1002294',
                    (u"Bruker Daltonics' EVOQ series."),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Bruker Daltonics timsTOF series', u'MS:1003123',
                    (u'Bruker Daltonics timsTOF series'),
                    'instrument model',
                    [u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'Shimadzu MALDI-TOF instrument model', u'MS:1000602',
                    (u'Shimadzu MALDI-TOF instrument model.'),
                    'instrument model',
                    [u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'Shimadzu Scientific Instruments instrument model', u'MS:1000603',
                    (u'Shimadzu Scientific Instruments instrument model.'),
                    'instrument model',
                    [u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'Auto Spec Ultima NT', u'MS:1000150',
                    (u'Waters magnetic sector based AutoSpec Ultima NT MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'GCT', u'MS:1000159',
                    (u'Waters oa-ToF based GCT.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'IsoPrime', u'MS:1000164',
                    (u'Waters IsoPrime MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'IsoProbe', u'MS:1000165',
                    (u'Waters IsoProbe MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'IsoProbe T', u'MS:1000166',
                    (u'Waters IsoProbe T MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'M@LDI L', u'MS:1000170',
                    (u'Waters oa-ToF based MALDI L.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'M@LDI LR', u'MS:1000171',
                    (u'Waters oa-ToF based MALDI LR.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'NG-5400', u'MS:1000180',
                    (u'Waters NG-5400 MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Platform ICP', u'MS:1000184',
                    (u'Waters Platform ICP MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Q-Tof micro', u'MS:1000188',
                    (u'Waters oa-ToF based Q-Tof micro.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Q-Tof Ultima', u'MS:1000189',
                    (u'Waters oa-ToF based Q-Tof Ultima.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'quattro micro', u'MS:1000191',
                    (u'Waters (triple) quadrupole based micro.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Quattro Ultima', u'MS:1000192',
                    (u'Waters (triple) quadrupole based Ultima.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Q-Tof Premier', u'MS:1000632',
                    (u'Waters oa-ToF based Q-Tof Premier.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Acquity UPLC PDA', u'MS:1000818',
                    (u'Acquity UPLC Photodiode Array Detector.'),
                    'instrument model',
                    [u'Waters instrument model', u'photodiode array detector', u'instrument model', u'array detector', u'detector type']),
    InstrumentModel(u'Acquity UPLC FLR', u'MS:1000819',
                    (u'Acquity UPLC Fluorescence Detector.'),
                    'instrument model',
                    [u'Waters instrument model', u'fluorescence detector', u'instrument model', u'detector type']),
    InstrumentModel(u'ACQUITY UPLC', u'MS:1001761',
                    (u'Waters LC-system ACQUITY UPLC.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'ACQUITY UPLC H-Class', u'MS:1001762',
                    (u'Waters LC-system ACQUITY UPLC H-Class.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'ACQUITY UPLC H-Class Bio', u'MS:1001763',
                    (u'Waters LC-system ACQUITY UPLC H-Class Bio.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'ACQUITY UPLC I-Class', u'MS:1001764',
                    (u'Waters LC-system ACQUITY UPLC I-Class.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'ACQUITY UPLC Systems with 2D Technology', u'MS:1001765',
                    (u'Waters LC-system ACQUITY UPLC Systems with 2D Technology.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'nanoACQUITY UPLC', u'MS:1001766',
                    (u'Waters LC-system nanoACQUITY UPLC.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'nanoACQUITY UPLC System with 1D Technology', u'MS:1001767',
                    (u'Waters LC-system nanoACQUITY UPLC System with 1D Technology.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'nanoACQUITY UPLC with HDX Technology', u'MS:1001768',
                    (u'Waters LC-system nanoACQUITY UPLC with HDX Technology.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'TRIZAIC UPLC nanoTile', u'MS:1001769',
                    (u'Waters LC-system TRIZAIC UPLC nanoTile.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'GCT Premier', u'MS:1001770',
                    (u'Waters oa-ToF based GCT Premier.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'MALDI Synapt G2 HDMS', u'MS:1001771',
                    (u'Waters oa-ToF based MALDI Synapt G2 HDMS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'MALDI Synapt G2 MS', u'MS:1001772',
                    (u'Waters oa-ToF based MALDI Synapt G2 MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'MALDI Synapt G2-S HDMS', u'MS:1001773',
                    (u'Waters oa-ToF based MALDI Synapt G2 MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'MALDI Synapt G2-S MS', u'MS:1001774',
                    (u'Waters oa-ToF based MALDI Synapt G2-S MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'MALDI Synapt HDMS', u'MS:1001775',
                    (u'Waters oa-ToF based MALDI Synapt HDMS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'MALDI Synapt MS', u'MS:1001776',
                    (u'Waters oa-ToF based MALDI Synapt MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Synapt G2 HDMS', u'MS:1001777',
                    (u'Waters oa-ToF based Synapt G2 HDMS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Synapt G2 MS', u'MS:1001778',
                    (u'Waters oa-ToF based Synapt G2 MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Synapt G2-S HDMS', u'MS:1001779',
                    (u'Waters oa-ToF based Synapt G2-S HDMS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Synapt G2-S MS', u'MS:1001780',
                    (u'Waters oa-ToF based Synapt G2-S MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Synapt HDMS', u'MS:1001781',
                    (u'Waters oa-ToF based Synapt HDMS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Synapt MS', u'MS:1001782',
                    (u'Waters oa-ToF based Synapt MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo G2 Q-Tof', u'MS:1001783',
                    (u'Waters oa-ToF based Xevo G2 Q-Tof.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo G2 Tof', u'MS:1001784',
                    (u'Waters oa-ToF based Xevo G2 Tof.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo Q-Tof', u'MS:1001785',
                    (u'Waters oa-ToF based Xevo Q-Tof.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'3100', u'MS:1001786',
                    (u'Waters quadrupole based 3100.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Acquity SQD', u'MS:1001787',
                    (u'Waters quadrupole based Acquity SQD.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Acquity TQD', u'MS:1001788',
                    (u'Waters quadrupole based Acquity TQD.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Quattro micro GC', u'MS:1001789',
                    (u'Waters (triple) quadrupole based Quattro micro GC.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo TQ MS', u'MS:1001790',
                    (u'Waters quadrupole based Xevo TQ MS.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo TQD', u'MS:1001791',
                    (u'Waters quadrupole based Xevo TQD.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo TQ-S', u'MS:1001792',
                    (u'Waters quadrupole based Xevo TQ-S.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'SQ Detector 2', u'MS:1002274',
                    (u'Waters quadrupole based SQ Detector 2.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo G2-S Tof', u'MS:1002275',
                    (u'Waters oa-ToF based Xevo G2-S Tof.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo G2-S QTof', u'MS:1002276',
                    (u'Waters oa-ToF based Xevo G2-S QTof.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'AutoSpec Premier', u'MS:1002277',
                    (u'Waters AutoSpec Premier magnetic sector instrument.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'SYNAPT G2-Si', u'MS:1002726',
                    (u'Waters Corporation SYNAPT G2-Si orthogonal acceleration'
                     u'time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'MALDI SYNAPT G2-Si', u'MS:1002727',
                    (u'Waters Corporation MALDI SYNAPT G2-Si orthogonal'
                     u'acceleration time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Vion IMS QTof', u'MS:1002728',
                    (u'Waters Corporation Vion IMS QTof orthogonal acceleration'
                     u'time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo G2 XS Tof', u'MS:1002729',
                    (u'Waters Corporation Xevo G2 XS Tof orthogonal acceleration'
                     u'time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo TQ-XS', u'MS:1002730',
                    (u'Waters Corporation Xevo TQ-XS triple quadrupole mass'
                     u'spectrometer.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Xevo TQ-S micro', u'MS:1002731',
                    (u'Waters Corporation Xevo TQ-S micro triple quadrupole mass'
                     u'spectrometer.'),
                    'instrument model',
                    [u'Waters instrument model', u'instrument model']),
    InstrumentModel(u'Thermo Finnigan instrument model', u'MS:1000125',
                    (u'ThermoFinnigan from Thermo Electron Corporation instrument'
                     u'model.'),
                    'instrument model',
                    [u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Thermo Electron instrument model', u'MS:1000492',
                    (u'Thermo Electron Corporation instrument model.'),
                    'instrument model',
                    [u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Finnigan MAT instrument model', u'MS:1000493',
                    (u'Finnigan MAT instrument model.'),
                    'instrument model',
                    [u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Thermo Scientific instrument model', u'MS:1000494',
                    (u'Thermo Scientific instrument model.'),
                    'instrument model',
                    [u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'IonSpec instrument model', u'MS:1000123',
                    (u'IonSpec corporation instrument model.'),
                    'instrument model',
                    [u'Varian instrument model', u'instrument model']),
    InstrumentModel(u'1200 series LC/MSD SL', u'MS:1000467',
                    (u'The 1200 Series LC/MSD SL ion trap belongs to the Agilent'
                     u'LC/MSD ion trap family. It provides fast polarity switching'
                     u'and multisignal data acquisition capabilities in a single'
                     u'run while also providing 5 stages of automated data'
                     u'dependent MS2 and 11 stages of manual MS2.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6110 Quadrupole LC/MS', u'MS:1000468',
                    (u'The 6110 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with an entry level'
                     u'single quadrupole mass spectrometer from the 6100 Series of'
                     u'Agilent quadrupole mass spectrometers. 6110 Quadrupole mass'
                     u'spectrometer has m/z range of 10-1500 and 2500 u/s scan'
                     u'speed. It proves useful for wide range of SIM quantitative'
                     u'applications.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6120A Quadrupole LC/MS', u'MS:1000469',
                    (u'The 6120A Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a single quadrupole'
                     u'mass spectrometer from the 6100 Series of Agilent mass'
                     u'spectrometers. 6120 quadrupole mass spectrometer has m/z'
                     u'range of 10-1500, 2500 u/s scan speed and utilizes multiple'
                     u'signal acquisition.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6130 Quadrupole LC/MS', u'MS:1000470',
                    (u'The 6130 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a single quadrupole'
                     u'mass spectrometer from the 6100 series of Agilent mass'
                     u'spectrometers. The 6130 quadrupole mass spectrometer has m/z'
                     u'range of 2-3000, 2500 u/s scan speed in standard mode and'
                     u'5250 u/s speed in fast-scan mode. It also uses multiple'
                     u'signal acquisition.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6140 Quadrupole LC/MS', u'MS:1000471',
                    (u'The 6140 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a single quadrupole'
                     u'mass spectrometer from the 6100 Series of Agilent quadrupole'
                     u'mass spectrometers. 6140 Quadrupole mass spectrometer has'
                     u'm/z range of 10-1350, 2500 u/s scan speed in standard mode'
                     u'and 10000 u/s speed in fast-scan mode. It also uses multiple'
                     u'signal acquisition.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6210 Time-of-Flight LC/MS', u'MS:1000472',
                    (u'The 6210 Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer. This time of flight mass'
                     u'spectrometer has a m/z range of 50-12000, mass accuracy of'
                     u'less than 2 ppm and resolution greater than 13,000 at m/z'
                     u'2722. It has multiple ion sources and can be used with'
                     u'multimode ion sources.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6310 Ion Trap LC/MS', u'MS:1000473',
                    (u'The 6310 Ion Trap LC/MS is a Agilent liquid chromatography'
                     u'instrument combined with a 6300 series Agilent ion trap. It'
                     u'has a mass range of 50-2200 between 0.6 to 0.35 resolution'
                     u'and mass range of 200-4000 with resolution of 3-4. The scan'
                     u'speed varies from 1650-27000 for the respective mass ranges.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6320 Ion Trap LC/MS', u'MS:1000474',
                    (u'The 6320 Ion Trap LC/MS is a Agilent liquid chromatography'
                     u'instrument combined with a 6300 series Agilent ion trap. It'
                     u'has a mass range of 50-2200 between 0.6 to 0.25 resolution'
                     u'and mass range of 200-4000 with resolution of less than 3.'
                     u'The scan speed varies from 1650-27000 for the respective'
                     u'mass ranges.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6330 Ion Trap LC/MS', u'MS:1000475',
                    (u'The 6330 Ion Trap LC/MS is a Agilent liquid chromatography'
                     u'instrument combined with a 6300 series Agilent ion trap. It'
                     u'has a mass range of 50-2200 between 0.6 to 0.25 resolution'
                     u'and mass range of 200-4000 with resolution of less than 3.'
                     u'The scan speed varies from 1650-27000 for the respective'
                     u'mass ranges.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6340 Ion Trap LC/MS', u'MS:1000476',
                    (u'The 6340 Ion Trap LC/MS is a Agilent liquid chromatography'
                     u'instrument combined with a 6300 series Agilent ion trap. It'
                     u'has a mass range of 50-2200 between 0.6 to 0.25 resolution'
                     u'and mass range of 200-4000 with resolution of less than 3.'
                     u'The scan speed varies from 1650-27000 for the respective'
                     u'mass ranges.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6410 Triple Quadrupole LC/MS', u'MS:1000477',
                    (u'The 6410 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer. Mass range of the mass'
                     u'spectrometer is 15-1650 m/z, resolution is at three settings'
                     u'of 0.7 u (unit), 1.2 u (wide) and 2.5 u (widest). The mass'
                     u'accuracy for 6410 mass spectrometer is 0.1 across the mass'
                     u'range. The collision cell is a hexapole with linear'
                     u'acceleration.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'1200 series LC/MSD VL', u'MS:1000478',
                    (u'The LC/MSD VL ion trap is part of the family of Agilent ion'
                     u'trap mass spectrometers. It has ESI, APCI and APPI ion'
                     u'sources and is a useful ion trap when the amount of sample'
                     u'is not the limiting factor.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6220 Time-of-Flight LC/MS', u'MS:1000675',
                    (u'The 6220 Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer. This time of flight mass'
                     u'spectrometer has a m/z range of 50-12000, mass accuracy of'
                     u'less than 2 ppm and resolution greater than 13,000 at m/z'
                     u'2722. It has multiple ion sources and can be used with'
                     u'multimode ion sources.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6510 Quadrupole Time-of-Flight LC/MS', u'MS:1000676',
                    (u'The 6510 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer. This time of flight mass'
                     u'spectrometer has a m/z range of 50-12000, mass accuracy of'
                     u'less than 2 ppm and resolution greater than 13,000 at m/z'
                     u'2722. It has multiple ion sources and can be used with'
                     u'multimode ion sources.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6520A Quadrupole Time-of-Flight LC/MS', u'MS:1000677',
                    (u'The 6520A Quadrupole Time-of-Flight LC/MS is a Agilent'
                     u'liquid chromatography instrument combined with a Agilent'
                     u'time of flight mass spectrometer. This time of flight mass'
                     u'spectrometer has a m/z range of 50-12000, mass accuracy of'
                     u'less than 2 ppm and resolution greater than 26,000 at m/z'
                     u'2722. It has multiple ion sources and can be used with'
                     u'multimode ion sources.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6420 Triple Quadrupole LC/MS', u'MS:1002444',
                    (u'The 6420 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6460 Triple Quadrupole LC/MS', u'MS:1002445',
                    (u'The 6460 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer. It is similar to the 6420 but'
                     u'adds Agilent Jet Stream (AJS) technology to increase'
                     u'sensitivity.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6490 Triple Quadrupole LC/MS', u'MS:1002446',
                    (u'The 6490 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer. It is similar to the 6420 but'
                     u'adds the Agilent iFunnel technology to increase sensitivity.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6550 iFunnel Q-TOF LC/MS', u'MS:1002783',
                    (u'The 6550 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6550A iFunnel Q-TOF LC/MS', u'MS:1002784',
                    (u'The 6550A Quadrupole Time-of-Flight LC/MS is a Agilent'
                     u'liquid chromatography instrument combined with a Agilent'
                     u'time of flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6520B Q-TOF LC/MS', u'MS:1002785',
                    (u'The 6520B Quadrupole Time-of-Flight LC/MS is a Agilent'
                     u'liquid chromatography instrument combined with a Agilent'
                     u'time of flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6530A Q-TOF LC/MS', u'MS:1002786',
                    (u'The 6530A Quadrupole Time-of-Flight LC/MS is a Agilent'
                     u'liquid chromatography instrument combined with a Agilent'
                     u'time of flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6530B Q-TOF LC/MS', u'MS:1002787',
                    (u'The 6530B Quadrupole Time-of-Flight LC/MS is a Agilent'
                     u'liquid chromatography instrument combined with a Agilent'
                     u'time of flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6538 Q-TOF LC/MS', u'MS:1002788',
                    (u'The 6538 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6540 Q-TOF LC/MS', u'MS:1002789',
                    (u'The 6540 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6542 Q-TOF LC/MS', u'MS:1002790',
                    (u'The 6542 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6545 Q-TOF LC/MS', u'MS:1002791',
                    (u'The 6545 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6560 Q-TOF LC/MS', u'MS:1002792',
                    (u'The 6560 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6570 Q-TOF LC/MS', u'MS:1002793',
                    (u'The 6570 Quadrupole Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6120B Quadrupole LC/MS', u'MS:1002794',
                    (u'The 6120B Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a single quadrupole'
                     u'mass spectrometer from the 6100 Series of Agilent mass'
                     u'spectrometers.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6150 Quadrupole LC/MS', u'MS:1002795',
                    (u'The 6150 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a single quadrupole'
                     u'mass spectrometer from the 6100 Series of Agilent mass'
                     u'spectrometers.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6224 Time-of-Flight LC/MS', u'MS:1002796',
                    (u'The 6224 Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6230A Time-of-Flight LC/MS', u'MS:1002797',
                    (u'The 6230A Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6230B Time-of-Flight LC/MS', u'MS:1002798',
                    (u'The 6230B Time-of-Flight LC/MS is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent time of'
                     u'flight mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6430 Triple Quadrupole LC/MS', u'MS:1002799',
                    (u'The 6430 Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6495A Triple Quadrupole LC/MS', u'MS:1002800',
                    (u'The 6495A Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'6495B Triple Quadrupole LC/MS', u'MS:1002801',
                    (u'The 6495B Quadrupole LC/MS system is a Agilent liquid'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'7000A Triple Quadrupole GC/MS', u'MS:1002802',
                    (u'The 7000A Quadrupole GC/MS system is a Agilent gas'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'7000B Triple Quadrupole GC/MS', u'MS:1002803',
                    (u'The 7000B Quadrupole GC/MS system is a Agilent gas'
                     u'chromatography instrument combined with a Agilent triple'
                     u'quadrupole mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'7800 Quadrupole ICP-MS', u'MS:1002804',
                    (u'The 7800 Quadrupole ICP-MS system is a Agilent inductively'
                     u'couple plasma instrument combined with a Agilent quadrupole'
                     u'mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'8800 Triple Quadrupole ICP-MS', u'MS:1002805',
                    (u'The 8800 Quadrupole ICP-MS system is a Agilent inductively'
                     u'couple plasma instrument combined with a Agilent quadrupole'
                     u'mass spectrometer.'),
                    'instrument model',
                    [u'Agilent instrument model', u'instrument model']),
    InstrumentModel(u'4700 Proteomics Analyzer', u'MS:1000140',
                    (u'Applied Biosystems/MDS SCIEX 4700 Proteomics Analyzer MS.'),
                    'instrument model',
                    [u'Applied Biosystems instrument model', u'instrument model']),
    InstrumentModel(u'Voyager-DE PRO', u'MS:1000203',
                    (u'Applied Biosystems/MDS SCIEX Voyager-DE PRO MS.'),
                    'instrument model',
                    [u'Applied Biosystems instrument model', u'instrument model']),
    InstrumentModel(u'Voyager-DE STR', u'MS:1000204',
                    (u'Applied Biosystems/MDS SCIEX Voyager-DE STR MS.'),
                    'instrument model',
                    [u'Applied Biosystems instrument model', u'instrument model']),
    InstrumentModel(u'4800 Proteomics Analyzer', u'MS:1000658',
                    (u'Applied Biosystems|MDS SCIEX 4800 Proteomics Analyzer.'),
                    'instrument model',
                    [u'Applied Biosystems instrument model', u'instrument model']),
    InstrumentModel(u'Pegasus HRT', u'MS:1001801',
                    (u'LECO high resolution time-of-flight GC mass spectrometer.'),
                    'instrument model',
                    [u'LECO instrument model', u'instrument model']),
    InstrumentModel(u'Citius HRT', u'MS:1001802',
                    (u'LECO high resolution time-of-flight LC mass spectrometer.'),
                    'instrument model',
                    [u'LECO instrument model', u'instrument model']),
    InstrumentModel(u'Pegasus', u'MS:1001803',
                    (u'LECO GC time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'LECO instrument model', u'instrument model']),
    InstrumentModel(u'TruTOF', u'MS:1001804',
                    (u'LECO bench-top GC time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'LECO instrument model', u'instrument model']),
    InstrumentModel(u'Pegasus 4D', u'MS:1001945',
                    (u'LECO nominal mass resolution time-of-flight GCxGC mass'
                     u'spectrometer.'),
                    'instrument model',
                    [u'LECO instrument model', u'instrument model']),
    InstrumentModel(u'Pegasus III', u'MS:1002278',
                    (u'LECO nominal mass resolution time-of-flight GC mass'
                     u'spectrometer.'),
                    'instrument model',
                    [u'LECO instrument model', u'instrument model']),
    InstrumentModel(u'Pegasus BT', u'MS:1002719',
                    (u'LECO bench-top GC time-of-flight mass spectrometer.'),
                    'instrument model',
                    [u'LECO instrument model', u'instrument model']),
    InstrumentModel(u'HCT', u'MS:1000160',
                    (u"Bruker Daltonics' HCT: ESI Q-TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    [u'Bruker Daltonics HCT Series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'HCTplus', u'MS:1000161',
                    (u"Bruker Daltonics' HCTplus: ESI Q-TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    [u'Bruker Daltonics HCT Series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'HCTultra', u'MS:1000698',
                    (u"Bruker Daltonics' HCTultra: ESI TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    [u'Bruker Daltonics HCT Series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'HCTultra PTM', u'MS:1000699',
                    (u"Bruker Daltonics' HCTultra PTM: ESI TOF, Nanospray, APCI,"
                     u'APPI, PTR.'),
                    'instrument model',
                    [u'Bruker Daltonics HCT Series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'HCTultra ETD II', u'MS:1000700',
                    (u"Bruker Daltonics' HCTultra ETD II: ESI Q-TOF, Nanospray,"
                     u'APCI, APPI, ETD.'),
                    'instrument model',
                    [u'Bruker Daltonics HCT Series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'esquire 4000', u'MS:1000156',
                    (u"Bruker Daltonics' esquire 4000: linear ion trap, ESI, MALDI,"
                     u'Nanospray, APCI, APPI.'),
                    'instrument model',
                    [u'Bruker Daltonics esquire series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'esquire 6000', u'MS:1000157',
                    (u"Bruker Daltonics' esquire 6000: linear ion trap, ESI, MALDI,"
                     u'Nanospray, APCI, APPI.'),
                    'instrument model',
                    [u'Bruker Daltonics esquire series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'autoflex II', u'MS:1000148',
                    (u"Bruker Daltonics' autoflex II: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'autoflex TOF/TOF', u'MS:1000149',
                    (u"Bruker Daltonics' autoflex TOF/TOF MS: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'microflex', u'MS:1000177',
                    (u"Bruker Daltonics' microflex: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'OmniFlex', u'MS:1000183',
                    (u"Bruker Daltonics' OmniFlex: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'ultraflex', u'MS:1000201',
                    (u"Bruker Daltonics' ultraflex: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'ultraflex TOF/TOF', u'MS:1000202',
                    (u"Bruker Daltonics' ultraflex TOF/TOF: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'autoflex III smartbeam', u'MS:1000696',
                    (u"Bruker Daltonics' autoflex III smartbeam: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'microflex LT', u'MS:1000701',
                    (u"Bruker Daltonics' microflex LT: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'ultraflex III TOF/TOF', u'MS:1000705',
                    (u"Bruker Daltonics' ultraflex III TOF/TOF: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'microflex LRF', u'MS:1001543',
                    (u"Bruker Daltonics' microflex LRF: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'ultrafleXtreme', u'MS:1001544',
                    (u"Bruker Daltonics' ultrafleXtreme: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'microflex II', u'MS:1001550',
                    (u"Bruker Daltonics' microflex II: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'autoflex II TOF/TOF', u'MS:1001553',
                    (u"Bruker Daltonics' autoflex II TOF/TOF: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'autoflex III TOF/TOF smartbeam', u'MS:1001554',
                    (u"Bruker Daltonics' autoflex III TOF/TOF smartbeam: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'autoflex', u'MS:1001555',
                    (u"Bruker Daltonics' autoflex: MALDI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'rapifleX', u'MS:1003122',
                    (u"Bruker Daltonics' rapiflex: MALDI TOF/TOF."),
                    'instrument model',
                    [u'Bruker Daltonics flex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'BioTOF II', u'MS:1000151',
                    (u"Bruker Daltonics' BioTOF II: ESI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics BioTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'BioTOF-Q', u'MS:1000152',
                    (u"Bruker Daltonics' BioTOF-Q: ESI Q-TOF."),
                    'instrument model',
                    [u'Bruker Daltonics BioTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'BioTOF', u'MS:1001537',
                    (u"Bruker Daltonics' BioTOF: ESI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics BioTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'BioTOF III', u'MS:1001538',
                    (u"Bruker Daltonics' BioTOF III: ESI TOF."),
                    'instrument model',
                    [u'Bruker Daltonics BioTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'UltroTOF-Q', u'MS:1001539',
                    (u"Bruker Daltonics' UltroTOF-Q: ESI Q-TOF (MALDI optional)."),
                    'instrument model',
                    [u'Bruker Daltonics BioTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'microTOF LC', u'MS:1000178',
                    (u"Bruker Daltonics' microTOF LC: ESI TOF, Nanospray, APCI,"
                     u'APPI.'),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'micrOTOF', u'MS:1000702',
                    (u"Bruker Daltonics' micrOTOF: ESI TOF, APCI, APPI."),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'micrOTOF-Q', u'MS:1000703',
                    (u"Bruker Daltonics' micrOTOF-Q: ESI Q-TOF, Nanospray, APCI,"
                     u'APPI.'),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'micrOTOF-Q II', u'MS:1000704',
                    (u"Bruker Daltonics' micrOTOF-Q II: ESI Q-TOF, Nanospray, APCI,"
                     u'APPI.'),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'micrOTOF II', u'MS:1001540',
                    (u"Bruker Daltonics' micrOTOF II: ESI TOF, Nanospray, APCI,"
                     u'APPI.'),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'impact', u'MS:1002077',
                    (u"Bruker Daltonics' impact: ESI Q-TOF, Nanospray, APCI, APPI,"
                     u'GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'compact', u'MS:1002280',
                    (u"Bruker Daltonics' compact: ESI Q-TOF, Nanospray, APCI, APPI,"
                     u'GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'micrOTOF-Q III', u'MS:1002299',
                    (u"Bruker Daltonics' micrOTOF-Q III: ESI Q-TOF, Nanospray,"
                     u'APCI, APPI, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'impact II', u'MS:1002666',
                    (u"Bruker Daltonics' impact II."),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'impact HD', u'MS:1002667',
                    (u"Bruker Daltonics' impact HD."),
                    'instrument model',
                    [u'Bruker Daltonics micrOTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'amaZon ETD', u'MS:1001542',
                    (u"Bruker Daltonics' amaZon ETD: ESI quadrupole ion trap,"
                     u'Nanospray, APCI, APPI, ETD, PTR.'),
                    'instrument model',
                    [u'Bruker Daltonics amaZon series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'amaZon X', u'MS:1001546',
                    (u"Bruker Daltonics' amaZon X: ESI quadrupole ion trap, APCI,"
                     u'APPI, ETD, PTR.'),
                    'instrument model',
                    [u'Bruker Daltonics amaZon series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'amaZon Speed ETD', u'MS:1002300',
                    (u"Bruker Daltonics' amaZon Speed ETD: ESI quadrupole ion trap,"
                     u'Nanospray, APCI, APPI, ETD, PTR, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    [u'Bruker Daltonics amaZon series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'amaZon Speed', u'MS:1002301',
                    (u"Bruker Daltonics' amaZon ETD: ESI quadrupole ion trap,"
                     u'Nanospray, APCI, APPI, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    [u'Bruker Daltonics amaZon series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'maXis', u'MS:1001541',
                    (u"Bruker Daltonics' maXis: ESI Q-TOF, Nanospray, APCI, APPI."),
                    'instrument model',
                    [u'Bruker Daltonics maXis series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'maXis 4G', u'MS:1002279',
                    (u"Bruker Daltonics' maXis 4G: ESI Q-TOF, Nanospray, APCI,"
                     u'APPI, GC-APCI, CaptiveSpray.'),
                    'instrument model',
                    [u'Bruker Daltonics maXis series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'maXis II', u'MS:1003004',
                    (u"Bruker Daltonics' maXis II."),
                    'instrument model',
                    [u'Bruker Daltonics maXis series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'solariX', u'MS:1001549',
                    (u"Bruker Daltonics' solariX: ESI, MALDI, Qh-FT_ICR."),
                    'instrument model',
                    [u'Bruker Daltonics solarix series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'apex IV', u'MS:1000141',
                    (u"Bruker Daltonics' apex IV: ESI, MALDI, Nanospray, APCI,"
                     u'APPI, Qh-FT_ICR.'),
                    'instrument model',
                    [u'Bruker Daltonics apex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'apex Q', u'MS:1000142',
                    (u"Bruker Daltonics' apex Q: ESI, MALDI, Nanospray, APCI, APPI,"
                     u'Qh-FT_ICR.'),
                    'instrument model',
                    [u'Bruker Daltonics apex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'apex ultra', u'MS:1000695',
                    (u"Bruker Daltonics' apex ultra: ESI, MALDI, Nanospray, APCI,"
                     u'APPI, Qh-FT_ICR.'),
                    'instrument model',
                    [u'Bruker Daltonics apex series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'SCION SQ', u'MS:1002295',
                    (u"Bruker Daltonics' SCION SQ: GC-single quadrupole."),
                    'instrument model',
                    [u'Bruker Daltonics SCION series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'SCION TQ', u'MS:1002296',
                    (u"Bruker Daltonics' SCION TQ: GC-triple quadrupole."),
                    'instrument model',
                    [u'Bruker Daltonics SCION series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'EVOQ Elite', u'MS:1002297',
                    (u"Bruker Daltonics' EVOQ Elite: LC-triple quadrupole."),
                    'instrument model',
                    [u'Bruker Daltonics EVOQ series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'EVOQ Qube', u'MS:1002298',
                    (u"Bruker Daltonics' EVOQ Qube: LC-triple quadrupole."),
                    'instrument model',
                    [u'Bruker Daltonics EVOQ series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'timsTOF Pro', u'MS:1003005',
                    (u"Bruker Daltonics' timsTOF Pro."),
                    'instrument model',
                    [u'Bruker Daltonics timsTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'timsTOF fleX', u'MS:1003124',
                    (u"Bruker Daltonics' timsTOF fleX"),
                    'instrument model',
                    [u'Bruker Daltonics timsTOF series', u'Bruker Daltonics instrument model', u'instrument model']),
    InstrumentModel(u'AXIMA CFR MALDI-TOF', u'MS:1000607',
                    (u'Shimadzu Biotech AXIMA CFR MALDI-TOF MS.'),
                    'instrument model',
                    [u'Shimadzu MALDI-TOF instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'AXIMA-QIT', u'MS:1000608',
                    (u'Shimadzu Biotech AXIMA-QIT MS.'),
                    'instrument model',
                    [u'Shimadzu MALDI-TOF instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'AXIMA-CFR plus', u'MS:1000609',
                    (u'Shimadzu Biotech AXIMA-CFR plus MS.'),
                    'instrument model',
                    [u'Shimadzu MALDI-TOF instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'AXIMA Performance MALDI-TOF/TOF', u'MS:1000610',
                    (u'Shimadzu Biotech AXIMA Performance MALDI-TOF/TOF MS.'),
                    'instrument model',
                    [u'Shimadzu MALDI-TOF instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'AXIMA Confidence MALDI-TOF', u'MS:1000611',
                    (u'Shimadzu Biotech AXIMA Confidence MALDI-TOF (curved field'
                     u'reflectron) MS.'),
                    'instrument model',
                    [u'Shimadzu MALDI-TOF instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'AXIMA Assurance Linear MALDI-TOF', u'MS:1000612',
                    (u'Shimadzu Biotech AXIMA Assurance Linear MALDI-TOF MS.'),
                    'instrument model',
                    [u'Shimadzu MALDI-TOF instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'Shimadzu MALDI-7090', u'MS:1002382',
                    (u'Shimadzu MALDI-7090: MALDI-TOF-TOF.'),
                    'instrument model',
                    [u'Shimadzu MALDI-TOF instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-IT-TOF', u'MS:1000604',
                    (u'Shimadzu Scientific Instruments LCMS-IT-TOF MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-2010EV', u'MS:1000605',
                    (u'Shimadzu Scientific Instruments LCMS-2010EV MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-2010A', u'MS:1000606',
                    (u'Shimadzu Scientific Instruments LCMS-2010A MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-9030', u'MS:1002998',
                    (u'Shimadzu Scientific Instruments LCMS-9030 Q-TOF MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-8060', u'MS:1002999',
                    (u'Shimadzu Scientific Instruments LCMS-8060 MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-8050', u'MS:1003000',
                    (u'Shimadzu Scientific Instruments LCMS-8050 MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-8045', u'MS:1003001',
                    (u'Shimadzu Scientific Instruments LCMS-8045 MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-8040', u'MS:1003002',
                    (u'Shimadzu Scientific Instruments LCMS-8040 MS.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'LCMS-2020', u'MS:1003003',
                    (u'Shimadzu Scientific Instruments LCMS-2020.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'GCMS-QP2010SE', u'MS:1003152',
                    (u'Shimadzu Scientific Instruments GCMS-QP2010SE.'),
                    'instrument model',
                    [u'Shimadzu Scientific Instruments instrument model', u'Shimadzu instrument model', u'instrument model']),
    InstrumentModel(u'DELTA plusAdvantage', u'MS:1000153',
                    (u'ThermoFinnigan DELTA plusAdvantage MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'DELTAplusXP', u'MS:1000154',
                    (u'ThermoFinnigan DELTAplusXP MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LCQ Advantage', u'MS:1000167',
                    (u'ThermoFinnigan LCQ Advantage MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LCQ Classic', u'MS:1000168',
                    (u'ThermoFinnigan LCQ Classic MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LCQ Deca XP Plus', u'MS:1000169',
                    (u'ThermoFinnigan LCQ Deca XP Plus MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'neptune', u'MS:1000179',
                    (u'ThermoFinnigan NEPTUNE MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'PolarisQ', u'MS:1000185',
                    (u'ThermoFinnigan PolarisQ MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Surveyor MSQ', u'MS:1000193',
                    (u'ThermoFinnigan Surveyor MSQ MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TEMPUS TOF', u'MS:1000196',
                    (u'ThermoFinnigan TEMPUS TOF MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TRACE DSQ', u'MS:1000197',
                    (u'ThermoFinnigan TRACE DSQ MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TRITON', u'MS:1000198',
                    (u'ThermoFinnigan TRITON MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Quantum', u'MS:1000199',
                    (u'ThermoFinnigan TSQ Quantum MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LCQ Deca', u'MS:1000554',
                    (u'ThermoFinnigan LCQ Deca.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'GC Quantum', u'MS:1000558',
                    (u'GC Quantum.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LCQ Fleet', u'MS:1000578',
                    (u'LCQ Fleet.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'DSQ', u'MS:1000634',
                    (u'ThermoFinnigan DSQ GC-MS.'),
                    'instrument model',
                    [u'Thermo Finnigan instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'MAT253', u'MS:1000172',
                    (u'ThermoFinnigan MAT253 MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'MAT900XP', u'MS:1000173',
                    (u'ThermoFinnigan MAT900XP MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'MAT900XP Trap', u'MS:1000174',
                    (u'ThermoFinnigan MAT900XP Trap MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'MAT95XP', u'MS:1000175',
                    (u'ThermoFinnigan MAT95XP MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'MAT95XP Trap', u'MS:1000176',
                    (u'ThermoFinnigan MAT95XP Trap MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'SSQ 7000', u'MS:1000748',
                    (u'ThermoFinnigan SSQ 7000 MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ 7000', u'MS:1000749',
                    (u'ThermoFinnigan TSQ 7000 MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ', u'MS:1000750',
                    (u'ThermoFinnigan TSQ MS.'),
                    'instrument model',
                    [u'Finnigan MAT instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ', u'MS:1000447',
                    (u'Finnigan LTQ MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ FT', u'MS:1000448',
                    (u'Finnigan LTQ FT MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap', u'MS:1000449',
                    (u'Finnigan LTQ Orbitrap MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LXQ', u'MS:1000450',
                    (u'Finnigan LXQ MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap Discovery', u'MS:1000555',
                    (u'LTQ Orbitrap Discovery.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap XL', u'MS:1000556',
                    (u'LTQ Orbitrap XL.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ FT Ultra', u'MS:1000557',
                    (u'LTQ FT Ultra.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Surveyor PDA', u'MS:1000622',
                    (u'Surveyor PDA.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Accela PDA', u'MS:1000623',
                    (u'Accela PDA.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'ITQ 700', u'MS:1000635',
                    (u'Thermo Scientific ITQ 700 GC-MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'ITQ 900', u'MS:1000636',
                    (u'Thermo Scientific ITQ 900 GC-MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'ITQ 1100', u'MS:1000637',
                    (u'Thermo Scientific ITQ 1100 GC-MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ XL ETD', u'MS:1000638',
                    (u'Thermo Scientific LTQ XL MS with ETD.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap XL ETD', u'MS:1000639',
                    (u'Thermo Scientific LTQ Orbitrap XL MS with ETD.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'DFS', u'MS:1000640',
                    (u'Thermo Scientific DFS HR GC-MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'DSQ II', u'MS:1000641',
                    (u'Thermo Scientific DSQ II GC-MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'MALDI LTQ XL', u'MS:1000642',
                    (u'Thermo Scientific MALDI LTQ XL MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'MALDI LTQ Orbitrap', u'MS:1000643',
                    (u'Thermo Scientific MALDI LTQ Orbitrap MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Quantum Access', u'MS:1000644',
                    (u'Thermo Scientific TSQ Quantum Access MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Element XR', u'MS:1000645',
                    (u'Thermo Scientific Element XR HR-ICP-MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Element 2', u'MS:1000646',
                    (u'Thermo Scientific Element 2 HR-ICP-MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Element GD', u'MS:1000647',
                    (u'Thermo Scientific Element GD Glow Discharge MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'GC IsoLink', u'MS:1000648',
                    (u'Thermo Scientific GC IsoLink Isotope Ratio MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Exactive', u'MS:1000649',
                    (u'Thermo Scientific Exactive MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Quantum Ultra AM', u'MS:1000743',
                    (u'Thermo Scientific TSQ Quantum Ultra AM.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Quantum Ultra', u'MS:1000751',
                    (u'Thermo Scientific TSQ Quantum Ultra.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ XL', u'MS:1000854',
                    (u'Thermo Scientific LTQ XL MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Velos', u'MS:1000855',
                    (u'Thermo Scientific LTQ Velos MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Velos ETD', u'MS:1000856',
                    (u'Thermo Scientific LTQ Velos MS with ETD.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Vantage', u'MS:1001510',
                    (u'TSQ Vantage.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap Velos', u'MS:1001742',
                    (u'Finnigan LTQ Orbitrap Velos MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'ISQ', u'MS:1001908',
                    (u'Thermo Scientific ISQ single quadrupole MS with the'
                     u'ExtractraBrite source.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Velos Plus', u'MS:1001909',
                    (u'Thermo Scientific second generation Velos.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap Elite', u'MS:1001910',
                    (u'Thermo Scientific LTQ Orbitrap Elite, often just referred to'
                     u'as the Orbitrap Elite.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Q Exactive', u'MS:1001911',
                    (u'Thermo Scientific Q Exactive.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap Fusion', u'MS:1002416',
                    (u'Thermo Scientific Orbitrap Fusion.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap Fusion ETD', u'MS:1002417',
                    (u'Thermo Scientific Orbitrap Fusion with ETD.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Quantiva', u'MS:1002418',
                    (u'Thermo Scientific TSQ Quantiva MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Endura', u'MS:1002419',
                    (u'Thermo Scientific TSQ Endura MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Q Exactive HF', u'MS:1002523',
                    (u'Thermo Scientific Q Exactive.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ 8000 Evo', u'MS:1002525',
                    (u'Thermo Scientific TSQ 8000 Evo MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Exactive Plus', u'MS:1002526',
                    (u'Thermo Scientific Exactive Plus MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Q Exactive Plus', u'MS:1002634',
                    (u'Thermo Scientific Q Exactive Plus.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap Fusion Lumos', u'MS:1002732',
                    (u'Thermo Scientific Orbitrap Fusion Lumos mass spectrometer'
                     u'with Tribrid architecture consisting of quadrupole mass'
                     u'filter, linear ion trap and Orbitrap mass analyzers.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap Classic', u'MS:1002835',
                    (u'Thermo Fisher Scientific LTQ Orbitrap Classic.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Altis', u'MS:1002874',
                    (u'Thermo Scientific TSQ Altis Triple Quadrupole MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ Quantis', u'MS:1002875',
                    (u'Thermo Scientific TSQ Quantis Triple Quadrupole MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'TSQ 9000', u'MS:1002876',
                    (u'Thermo Scientific TSQ 9000 Triple Quadrupole MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Q Exactive HF-X', u'MS:1002877',
                    (u'Thermo Scientific Q Exactive HF-X Hybrid Quadrupole Orbitrap'
                     u'MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap Exploris 480', u'MS:1003028',
                    (u'Thermo Scientific Orbitrap Exploris 480 Quadrupole Orbitrap'
                     u'MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap Eclipse', u'MS:1003029',
                    (u'Thermo Scientific Orbitrap Eclipse mass spectrometer with'
                     u'Tribrid architecture consisting of quadrupole mass filter,'
                     u'linear ion trap and Orbitrap mass analyzers.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap Exploris 240', u'MS:1003094',
                    (u'Thermo Scientific Orbitrap Exploris 240 Quadrupole Orbitrap'
                     u'MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap Exploris 120', u'MS:1003095',
                    (u'Thermo Scientific Orbitrap Exploris 120 Quadrupole Orbitrap'
                     u'MS.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'LTQ Orbitrap Velos Pro', u'MS:1003096',
                    (u'Thermo Scientific LTQ Orbitrap Velos Pro, often just'
                     u'referred to as the Orbitrap Velos Pro.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'Orbitrap ID-X', u'MS:1003112',
                    (u'Thermo Scientific Orbitrap ID-X mass spectrometer with'
                     u'Tribrid architecture consisting of quadrupole mass filter,'
                     u'linear ion trap and Orbitrap mass analyzers.'),
                    'instrument model',
                    [u'Thermo Scientific instrument model', u'Thermo Fisher Scientific instrument model', u'instrument model']),
    InstrumentModel(u'explorer', u'MS:1000158',
                    (u'IonSpec Explorer MS.'),
                    'instrument model',
                    [u'IonSpec instrument model', u'Varian instrument model', u'instrument model']),
    InstrumentModel(u'HiRes ESI', u'MS:1000162',
                    (u'IonSpec HiResESI MS.'),
                    'instrument model',
                    [u'IonSpec instrument model', u'Varian instrument model', u'instrument model']),
    InstrumentModel(u'HiRes MALDI', u'MS:1000163',
                    (u'IonSpec HiResMALDI MS.'),
                    'instrument model',
                    [u'IonSpec instrument model', u'Varian instrument model', u'instrument model']),
    InstrumentModel(u'OMEGA', u'MS:1000181',
                    (u'IonSpec OMEGA MS.'),
                    'instrument model',
                    [u'IonSpec instrument model', u'Varian instrument model', u'instrument model']),
    InstrumentModel(u'OMEGA-2001', u'MS:1000182',
                    (u'IonSpec OMEGA-2001 MS.'),
                    'instrument model',
                    [u'IonSpec instrument model', u'Varian instrument model', u'instrument model']),
    InstrumentModel(u'ultima', u'MS:1000200',
                    (u'IonSpec Ultima MS.'),
                    'instrument model',
                    [u'IonSpec instrument model', u'Varian instrument model', u'instrument model']),
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
