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
    Component(u'multiphoton ionization', u'MS:1000227',
              (u'Photoionization of an atom or molecule in which in two or'
               u'more photons are absorbed.'),
              'ionization type',
              [u'ionization type']),
    Component(u'fast ion bombardment', u'MS:1000446',
              (u'The ionization of any species by the interaction of a'
               u'focused beam of ions having a translational energy of'
               u'several thousand eV with a solid sample.'),
              'ionization type',
              [u'ionization type']),
    Component(u'resonance enhanced multiphoton ionization', u'MS:1000276',
              (u'Multiphoton ionization in which the ionization cross section'
               u'is significantly enhanced because the energy of the incident'
               u'photons is resonant with an intermediate excited state of'
               u'the neutral species.'),
              'ionization type',
              [u'ionization type']),
    Component(u'pyrolysis mass spectrometry', u'MS:1000274',
              (u'A mass spectrometry technique in which the sample is heated'
               u'to the point of decomposition and the gaseous decomposition'
               u'products are introduced into the ion source.'),
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
    Component(u'Negative Ion chemical ionization', u'MS:1000271',
              (u'Chemical ionization that results in the formation of'
               u'negative ions.'),
              'ionization type',
              [u'ionization type']),
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
    Component(u'flowing afterglow', u'MS:1000255',
              (u'An ion source immersed in a flow of helium or other inert'
               u'buffer gas that carries the ions through a meter-long'
               u'reactor at pressures around 100 Pa.'),
              'ionization type',
              [u'ionization type']),
    Component(u'desorption ionization', u'MS:1000247',
              (u'The formation of ions from a solid or liquid material after'
               u'the rapid vaporization of that sample.'),
              'ionization type',
              [u'ionization type']),
    Component(u'atmospheric pressure ionization', u'MS:1000240',
              (u'Any ionization process in which ions are formed in the gas'
               u'phase at atmospheric pressure.'),
              'ionization type',
              [u'ionization type']),
    Component(u'spark ionization', u'MS:1000404',
              (u'The formation of ions from a solid material by an'
               u'intermittent electrical discharge.'),
              'ionization type',
              [u'ionization type']),
    Component(u'thermal ionization', u'MS:1000407',
              (u'The ionization of a neutral species through contact with a'
               u'high temperature surface.'),
              'ionization type',
              [u'ionization type']),
    Component(u'surface ionization', u'MS:1000406',
              (u'The ionization of a neutral species when it interacts with a'
               u'solid surface with an appropriate work function and'
               u'temperature.'),
              'ionization type',
              [u'ionization type']),
    Component(u'plasma desorption ionization', u'MS:1000400',
              (u'The ionization of material in a solid sample by bombarding'
               u'it with ionic or neutral atoms formed as a result of the'
               u'fission of a suitable nuclide, typically 252Cf. Synonymous'
               u'with fission fragment ionization.'),
              'ionization type',
              [u'ionization type']),
    Component(u'soft ionization', u'MS:1000403',
              (u'The formation of gas-phase ions without extensive'
               u'fragmentation.'),
              'ionization type',
              [u'ionization type']),
    Component(u'secondary ionization', u'MS:1000402',
              (u'The process in which ions are ejected from a sample surface'
               u'as a result of bombardment by a primary beam of atoms or'
               u'ions.'),
              'ionization type',
              [u'ionization type']),
    Component(u'vertical ionization', u'MS:1000408',
              (u'A process in which an electron is removed from or added to a'
               u'molecule without a change in the positions of the atoms. The'
               u'resulting ion is typically in an excited vibrational state.'),
              'ionization type',
              [u'ionization type']),
    Component(u'autodetachment', u'MS:1000383',
              (u'The formation of a neutral when a negative ion in a discrete'
               u'state with an energy greater than the detachment threshold'
               u'loses an electron spontaneously without further interaction'
               u'with an energy source.'),
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
    Component(u'chemi-ionization', u'MS:1000386',
              (u'The reaction of a neutral molecule with an internally'
               u'excited molecule to form an ion. Note that this term is not'
               u'synonymous with chemical ionization.'),
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
    Component(u'surface-assisted laser desorption ionization', u'MS:1000405',
              (u'The formation of gas-phase ions from molecules that are'
               u'deposited on a particular surface substrate that is'
               u'irradiated with a pulsed laser. See also matrix-assisted'
               u'laser desorption ionization.'),
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
    Component(u'atmospheric pressure matrix-assisted laser desorption ionization', u'MS:1000239',
              (u'Matrix-assisted laser desorption ionization in which the'
               u'sample target is at atmospheric pressure and the ions formed'
               u'by the pulsed laser are sampled through a small aperture'
               u'into the mass spectrometer.'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'atmospheric pressure chemical ionization', u'MS:1000070',
              (u'Chemical ionization that takes place at atmospheric pressure'
               u'as opposed to the reduced pressure is normally used for'
               u'chemical ionization.'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'desorption electrospray ionization', u'MS:1002011',
              (u'Combination of electrospray and desorption ionization method'
               u'that ionizes gases, liquids and solids in open air under'
               u'atmospheric pressure." [DOI:10.1126/science.1104404'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
    Component(u'atmospheric pressure photoionization', u'MS:1000382',
              (u'Atmospheric pressure chemical ionization in which the'
               u'reactant ions are generated by photo-ionization.'),
              'ionization type',
              [u'atmospheric pressure ionization', u'ionization type']),
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
    Component(u'photomultiplier', u'MS:1000116',
              (u'A detector for conversion of the ion/electron signal into'
               u'photon(s) which are then amplified and detected.'),
              'detector type',
              [u'detector type']),
    Component(u'multi-collector', u'MS:1000115',
              (u'A detector system commonly used in inductively coupled'
               u'plasma mass spectrometers.'),
              'detector type',
              [u'detector type']),
    Component(u'faraday cup', u'MS:1000112',
              (u'A conducting cup or chamber that intercepts a charged'
               u'particle beam and is electrically connected to a current'
               u'measuring device.'),
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
    Component(u'electron multiplier', u'MS:1000253',
              (u'A device to amplify the current of a beam or packet of'
               u'charged particles or photons by incidence upon the surface'
               u'of an electrode to produce secondary electrons. The'
               u'secondary electrons are then accelerated to other electrodes'
               u'or parts of a continuous electrode to produce further'
               u'secondary electrons.'),
              'detector type',
              [u'detector type']),
    Component(u'fluorescence detector', u'MS:1002308',
              (u'A detector using a fluorescent signal after excitation with'
               u'light.'),
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
    Component(u'array detector', u'MS:1000345',
              (u'Detector comprising several ion collection elements,'
               u'arranged in a line or grid where each element is an'
               u'individual detector.'),
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
    Component(u'postacceleration detector', u'MS:1000351',
              (u'A detector in which the charged particles are accelerated to'
               u'a high velocity and impinge on a conversion dynode, emitting'
               u'secondary electrons. The electrons are accelerated onto a'
               u'phosphor screen, which emits photons that are in turn'
               u'detected using a photomultiplier or other photon detector.'),
              'detector type',
              [u'detector type']),
    Component(u'point collector', u'MS:1000350',
              (u'A detector in which the ion beam is focused onto a point and'
               u'the individual ions arrive sequentially.'),
              'detector type',
              [u'detector type']),
    Component(u'inductive detector', u'MS:1000624',
              (u'Inductive detector.'),
              'detector type',
              [u'detector type']),
    Component(u'electron multiplier tube', u'MS:1000111',
              (u'A device to amplify the current of a beam or packet of'
               u'charged particles or photons by incidence upon the surface'
               u'of an electrode to produce secondary electrons.'),
              'detector type',
              [u'electron multiplier', u'detector type']),
    Component(u'Acquity UPLC FLR', u'MS:1000819',
              (u'Acquity UPLC Fluorescence Detector.'),
              'detector type',
              [u'Waters instrument model', u'fluorescence detector', u'instrument model', u'detector type']),
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
    Component(u'focal plane array', u'MS:1000113',
              (u'An array of detectors for spatially disperse ion beams in'
               u'which all ions simultaneously impinge on the detector plane.'),
              'detector type',
              [u'focal plane collector', u'detector type']),
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
    Component(u'ion trap', u'MS:1000264',
              (u'A device for spatially confining ions using electric and'
               u'magnetic fields alone or in combination.'),
              'mass analyzer type',
              [u'mass analyzer type']),
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
    Component(u'electrostatic energy analyzer', u'MS:1000254',
              (u'A device consisting of conducting parallel plates,'
               u'concentric cylinders or concentric spheres that separates'
               u'charged particles according to their kinetic energy by means'
               u'of an electric field that is constant in time.'),
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
    Component(u'magnetic sector', u'MS:1000080',
              (u'A device that produces a magnetic field perpendicular to a'
               u'charged particle beam that deflects the beam to an extent'
               u'that is proportional to the particle momentum per unit'
               u'charge. For a monoenergetic beam, the deflection is'
               u'proportional to m/z.'),
              'mass analyzer type',
              [u'mass analyzer type']),
    Component(u'time-of-flight', u'MS:1000084',
              (u'Instrument that separates ions by m/z in a field-free region'
               u'after acceleration to a fixed acceleration energy.'),
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
    Component(u'continuous flow fast atom bombardment', u'MS:1000055',
              (u'Fast atom bombardment ionization in which the analyte in'
               u'solution is entrained in a flowing liquid matrix.'),
              'inlet type',
              [u'inlet type']),
    Component(u'reservoir', u'MS:1000067',
              (u'A sample inlet method involving a reservoir.'),
              'inlet type',
              [u'inlet type']),
    Component(u'particle beam', u'MS:1000066',
              (u'Method for generating ions from a solution of an analyte.'),
              'inlet type',
              [u'inlet type']),
    Component(u'open split', u'MS:1000065',
              (u'A division of flowing stream of liquid into two streams.'),
              'inlet type',
              [u'inlet type']),
    Component(u'moving wire', u'MS:1000064',
              (u'Continuous moving surface in the form of a wire which passes'
               u'through an ion source carrying analyte molecules.'),
              'inlet type',
              [u'inlet type']),
    Component(u'moving belt', u'MS:1000063',
              (u'Continuous moving surface in the form of a belt which passes'
               u'through an ion source carrying analyte molecules.'),
              'inlet type',
              [u'inlet type']),
    Component(u'membrane separator', u'MS:1000062',
              (u'A device to separate carrier molecules from analyte'
               u'molecules on the basis of ease of diffusion across a'
               u'semipermeable membrane.'),
              'inlet type',
              [u'inlet type']),
    Component(u'jet separator', u'MS:1000061',
              (u'A device that separates carrier gas from gaseous analyte'
               u'molecules on the basis of diffusivity.'),
              'inlet type',
              [u'inlet type']),
    Component(u'infusion', u'MS:1000060',
              (u'The continuous flow of solution of a sample into the'
               u'ionization source.'),
              'inlet type',
              [u'inlet type']),
    Component(u'thermospray inlet', u'MS:1000069',
              (u'A method for generating gas phase ions from a solution of an'
               u'analyte by rapid heating of the sample.'),
              'inlet type',
              [u'inlet type']),
    Component(u'septum', u'MS:1000068',
              (u'A disc composed of a flexible material that seals the'
               u'entrance to the reservoir. Can also be entrance to the'
               u'vacuum chamber.'),
              'inlet type',
              [u'inlet type']),
    Component(u'direct liquid introduction', u'MS:1000249',
              (u'The delivery of a liquid sample into a mass spectrometer for'
               u'spray or desorption ionization.'),
              'inlet type',
              [u'inlet type']),
    Component(u'direct insertion probe', u'MS:1000248',
              (u'A device for introducing a solid or liquid sample into a'
               u'mass spectrometer ion source for desorption ionization.'),
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

all_components_by_name = {c.name: c for c in all_components}


def component(name):
    try:
        return all_components_by_name[name]
    except KeyError:
        return Component(name, name, name, name, [name])


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

    def __init__(self, id=None, groups=None):
        self.id = id or uid()
        self.groups = sorted(groups or [], key=lambda x: x.order)
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


__all__ = [
    "InstrumentInformation", "ComponentGroup", "Component",
    "all_components_by_name", "ionization_types", "detector_types",
    "analyzer_types", "inlet_types", "component"
]
