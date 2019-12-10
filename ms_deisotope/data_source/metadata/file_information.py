"""Defines types for describing different kinds of mass spectrometry
data files and their contents, and a database of controlled vocabulary
terms for them.
"""

import os
import hashlib
import warnings

try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping

from six import string_types as basestring

from .cv import Term, TermSet

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError


class IDFormat(Term):
    """Describes a named spectrum identifier format, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`IDFormat` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


class FileFormat(Term):
    """Describes a named mass spectrometry data file format, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`FileFormat` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


class FileContent(Term):
    """Describes a named mass spectrometry data file content type,
    either using a controlled-vocabulary term or user-defined name.

    A :class:`FileContent` is equal to its name and its controlled
    vocabulary identifier.
    """
    pass


id_formats = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('native spectrum identifier format',
#             "id_formats", term_cls_name="IDFormat", writer=cog.out)
# ]]]
id_formats = TermSet([
    IDFormat(u'Agilent MassHunter nativeID format', u'MS:1001508',
             (u'Native format defined by scanId=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Shimadzu Biotech QTOF nativeID format', u'MS:1002898',
             (u'Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'spectrum from database string nativeID format', u'MS:1001532',
             (u'Native format defined by databasekey=xsd:string.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format', u'spectra data details', u'search input details']),
    IDFormat(u'spectrum from ProteinScape database nativeID format', u'MS:1001531',
             (u'Native format defined by databasekey=xsd:long.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format', u'spectra data details', u'search input details']),
    IDFormat(u'Bruker U2 nativeID format', u'MS:1000823',
             (u'Native format defined by declaration=xsd:nonNegativeInteger'
              u'collection=xsd:nonNegativeInteger'
              u'scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'no nativeID format', u'MS:1000824',
             (u'No nativeID format indicates that the file tagged with this'
              u'term does not contain spectra that can have a nativeID'
              u'format.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Mascot query number', u'MS:1001528',
             (u'Native format defined by query=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format', u'spectrum identification result details']),
    IDFormat(u'spectrum from database integer nativeID format', u'MS:1001526',
             (u'Native format defined by databasekey=xsd:long.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Shimadzu Biotech nativeID format', u'MS:1000929',
             (u'Native format defined by source=xsd:string'
              u'start=xsd:nonNegativeInteger end=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'UIMF nativeID format', u'MS:1002532',
             (u'Native format defined by frame=xsd:nonNegativeInteger'
              u'scan=xsd:nonNegativeInteger'
              u'frameType=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Bruker TDF nativeID format', u'MS:1002818',
             (u'Native format defined by frame=xsd:nonNegativeInteger'
              u'scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Bruker Container nativeID format', u'MS:1002303',
             (u'Native identifier (UUID).'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'SCIEX TOF/TOF nativeID format', u'MS:1001480',
             (u'Native format defined by jobRun=xsd:nonNegativeInteger'
              u'spotLabel=xsd:string spectrum=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Thermo nativeID format', u'MS:1000768',
             (u'Native format defined by'
              u'controllerType=xsd:nonNegativeInteger'
              u'controllerNumber=xsd:positiveInteger'
              u'scan=xsd:positiveInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Waters nativeID format', u'MS:1000769',
             (u'Native format defined by function=xsd:positiveInteger'
              u'process=xsd:nonNegativeInteger scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Scaffold nativeID format', u'MS:1001562',
             (u'Scaffold native ID format.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Bruker FID nativeID format', u'MS:1000773',
             (u'Native format defined by file=xsd:IDREF.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Bruker BAF nativeID format', u'MS:1000772',
             (u'Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'Bruker/Agilent YEP nativeID format', u'MS:1000771',
             (u'Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'WIFF nativeID format', u'MS:1000770',
             (u'Native format defined by sample=xsd:nonNegativeInteger'
              u'period=xsd:nonNegativeInteger cycle=xsd:nonNegativeInteger'
              u'experiment=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'spectrum identifier nativeID format', u'MS:1000777',
             (u'Native format defined by spectrum=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'scan number only nativeID format', u'MS:1000776',
             (u'Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'single peak list nativeID format', u'MS:1000775',
             (u'Native format defined by file=xsd:IDREF.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'multiple peak list nativeID format', u'MS:1000774',
             (u'Native format defined by index=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
    IDFormat(u'SCIEX TOF/TOF T2D nativeID format', u'MS:1001559',
             (u'Native format defined by file=xsd:IDREF.'),
             'native spectrum identifier format',
             [u'native spectrum identifier format']),
])
# [[[end]]]


file_formats = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('mass spectrometer file format',
#             "file_formats", term_cls_name="FileFormat", writer=cog.out)
# ]]]
file_formats = TermSet([
    FileFormat(u'Agilent MassHunter format', u'MS:1001509',
               (u'A data file format found in an Agilent MassHunter directory'
                u'which contains raw data acquired by an Agilent mass'
                u'spectrometer.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'msalign format', u'MS:1002899',
               (u'msalign file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'SCiLS Lab format', u'MS:1002385',
               (u'SCiLS Lab file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'chrom format', u'MS:1002966',
               (u'The Lipid Data Analyzer native chrom format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Mascot MGF format', u'MS:1001062',
               (u'Mascot MGF file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Andi-MS format', u'MS:1002441',
               (u'AIA Analytical Data Interchange file format for mass'
                u'spectrometry data.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bruker FID format', u'MS:1000825',
               (u'Bruker FID file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Proteinscape spectra', u'MS:1001527',
               (u'Spectra from Bruker/Protagen Proteinscape database.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'MS2 format', u'MS:1001466',
               (u'MS2 file format for MS2 spectral data." [PMID:15317041,'
                u'DOI:10.1002/rcm.1603'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'PerSeptive PKS format', u'MS:1001245',
               (u'PerSeptive peak list file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'SCIEX API III format', u'MS:1001246',
               (u'PE SCIEX peak list file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bruker XML format', u'MS:1001247',
               (u'Bruker data exchange XML format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Shimadzu Biotech database entity', u'MS:1000930',
               (u'Shimadzu Biotech format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Shimadzu Biotech LCD format', u'MS:1003009',
               (u'Shimadzu Biotech LCD file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'mz5 format', u'MS:1001881',
               (u'mz5 file format, modelled after mzML.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'mzMLb format', u'MS:1002838',
               (u'mzMLb file format, mzML encapsulated within HDF5." [PSI:PI'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'UIMF format', u'MS:1002531',
               (u'SQLite-based file format created at Pacific Northwest'
                u'National Lab. It stores an intermediate analysis of ion-'
                u'mobility mass spectrometry data.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bruker TDF format', u'MS:1002817',
               (u'Bruker TDF raw file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Waters raw format', u'MS:1000526',
               (u'Waters data file format found in a Waters RAW directory,'
                u'generated from an MS acquisition.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'feature format', u'MS:1002900',
               (u'TopFD feature file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bruker Container format', u'MS:1002302',
               (u'Bruker Container raw file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Phenyx XML format', u'MS:1001463',
               (u'Phenyx open XML file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'intermediate analysis format', u'file format']),
    FileFormat(u'mzML format', u'MS:1000584',
               (u'Proteomics Standards Inititative mzML file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bioworks SRF format', u'MS:1000742',
               (u'Thermo Finnigan SRF file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'intermediate analysis format', u'file format']),
    FileFormat(u'parameter file', u'MS:1000740',
               (u'Parameter file used to configure the acquisition of raw data'
                u'on the instrument.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Andromeda:apl file format', u'MS:1002996',
               (u'Peak list file format of the Andromeda search engine.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'ISB mzXML format', u'MS:1000566',
               (u'Institute of Systems Biology mzXML file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bruker/Agilent YEP format', u'MS:1000567',
               (u'Bruker/Agilent YEP file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'PSI mzData format', u'MS:1000564',
               (u'Proteomics Standards Inititative mzData file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Micromass PKL format', u'MS:1000565',
               (u'Micromass PKL file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'ABI WIFF format', u'MS:1000562',
               (u'Applied Biosystems WIFF file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Thermo RAW format', u'MS:1000563',
               (u'Thermo Scientific RAW file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'SCIEX TOF/TOF database', u'MS:1001481',
               (u'Applied Biosystems/MDS Analytical Technologies TOF/TOF'
                u'instrument database.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'DTA format', u'MS:1000613',
               (u'SEQUEST DTA file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'ProteinLynx Global Server mass spectrum XML format', u'MS:1000614',
               (u'Peak list file format used by ProteinLynx Global Server.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'text format', u'MS:1001369',
               (u'Simple text file format of \\"m/z [intensity]\\" values for a'
                u'PMF (or single MS2) search.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'SCIEX TOF/TOF T2D format', u'MS:1001560',
               (u'Applied Biosystems/MDS Analytical Technologies TOF/TOF'
                u'instrument export format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'MS1 format', u'MS:1002597',
               (u'MS1 file format for MS1 spectral data." [PMID:15317041'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bruker U2 format', u'MS:1000816',
               (u'Bruker HyStar U2 file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
    FileFormat(u'Bruker BAF format', u'MS:1000815',
               (u'Bruker BAF raw file format.'),
               'mass spectrometer file format',
               [u'mass spectrometer file format', u'file format']),
])
# [[[end]]]


content_keys = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('data file content',
#             "content_keys", term_cls_name="FileContent", writer=cog.out)
# ]]]
content_keys = TermSet([
    FileContent(u'electromagnetic radiation spectrum', u'MS:1000804',
                (u'A plot of the relative intensity of electromagnetic'
                 u'radiation as a function of the wavelength.'),
                'data file content',
                [u'data file content', u'spectrum type']),
    FileContent(u'absorption spectrum', u'MS:1000806',
                (u'A plot of the relative intensity of electromagnetic'
                 u'radiation absorbed by atoms or molecules when excited.'),
                'data file content',
                [u'data file content', u'spectrum type']),
    FileContent(u'emission spectrum', u'MS:1000805',
                (u'A plot of the relative intensity of electromagnetic'
                 u'radiation emitted by atoms or molecules when excited.'),
                'data file content',
                [u'data file content', u'spectrum type']),
    FileContent(u'mass spectrum', u'MS:1000294',
                (u'A plot of the relative abundance of a beam or other'
                 u'collection of ions as a function of the mass-to-charge ratio'
                 u'(m/z).'),
                'data file content',
                [u'data file content', u'spectrum type']),
    FileContent(u'PDA spectrum', u'MS:1000620',
                (u'OBSOLETE Spectrum generated from a photodiode array detector'
                 u'(ultraviolet/visible spectrum).'),
                'data file content',
                [u'data file content', u'spectrum type']),
    FileContent(u'mass chromatogram', u'MS:1000810',
                (u'A plot of the relative abundance of a beam or other'
                 u'collection of ions as a function of the retention time.'),
                'data file content',
                [u'data file content', u'chromatogram type']),
    FileContent(u'electromagnetic radiation chromatogram', u'MS:1000811',
                (u'The measurement of electromagnetic properties as a function'
                 u'of the retention time.'),
                'data file content',
                [u'data file content', u'chromatogram type']),
    FileContent(u'MSn spectrum', u'MS:1000580',
                (u'MSn refers to multi-stage MS2 experiments designed to record'
                 u'product ion spectra where n is the number of product ion'
                 u'stages (progeny ions). For ion traps, sequential MS/MS'
                 u'experiments can be undertaken where n > 2 whereas for a'
                 u'simple triple quadrupole system n=2. Use the term ms level'
                 u'(MS:1000511) for specifying n.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'CRM spectrum', u'MS:1000581',
                (u'Spectrum generated from MSn experiment with three or more'
                 u'stages of m/z separation and in which a particular multi-'
                 u'step reaction path is monitored.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'SIM spectrum', u'MS:1000582',
                (u'Spectrum obtained with the operation of a mass spectrometer'
                 u'in which the abundances of one ion or several ions of'
                 u'specific m/z values are recorded rather than the entire mass'
                 u'spectrum (Selected Ion Monitoring).'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'SRM spectrum', u'MS:1000583',
                (u'Spectrum obtained when data are acquired from specific'
                 u'product ions corresponding to m/z values of selected'
                 u'precursor ions a recorded via two or more stages of mass'
                 u'spectrometry. The precursor/product ion pair is called a'
                 u'transition pair. Data can be obtained for a single'
                 u'transition pair or multiple transition pairs. Multiple time'
                 u'segments of different transition pairs can exist in a single'
                 u'file. Single precursor ions can have multiple product ions'
                 u'consitituting multiple transition pairs. Selected reaction'
                 u'monitoring can be performed as tandem mass spectrometry in'
                 u'time or tandem mass spectrometry in space.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'e/2 mass spectrum', u'MS:1000328',
                (u'A mass spectrum obtained using a sector mass spectrometer in'
                 u'which the electric sector field E is set to half the value'
                 u'required to transmit the main ion-beam. This spectrum'
                 u'records the signal from doubly charged product ions of'
                 u'charge-stripping reactions.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'constant neutral gain spectrum', u'MS:1000325',
                (u'A spectrum formed of all product ions that have been'
                 u'produced by gain of a pre-selected neutral mass following'
                 u'the reaction with and addition of the gas in a collision'
                 u'cell.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'constant neutral loss spectrum', u'MS:1000326',
                (u'A spectrum formed of all product ions that have been'
                 u'produced with a selected m/z decrement from any precursor'
                 u'ions. The spectrum shown correlates to the precursor ion'
                 u'spectrum. See also neutral loss spectrum.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'charge inversion mass spectrum', u'MS:1000322',
                (u'The measurement of the relative abundance of ions that'
                 u'result from a charge inversion reaction as a function of'
                 u'm/z.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'product ion spectrum', u'MS:1000343',
                (u'OBSOLETE A mass spectrum recorded from any spectrometer in'
                 u'which the appropriate m/z separation scan function is set to'
                 u'record the product ion or ions of selected precursor ions.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'precursor ion spectrum', u'MS:1000341',
                (u'Spectrum generated by scanning precursor m/z while'
                 u'monitoring a fixed product m/z.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'MS1 spectrum', u'MS:1000579',
                (u'Mass spectrum created by a single-stage MS experiment or the'
                 u'first stage of a multi-stage experiment.'),
                'data file content',
                [u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'total ion current chromatogram', u'MS:1000235',
                (u'Chromatogram obtained by plotting the total ion current'
                 u'detected in each of a series of mass spectra recorded as a'
                 u'function of retention time.'),
                'data file content',
                [u'mass chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'selected reaction monitoring chromatogram', u'MS:1001473',
                (u'Chromatogram created by creating an array of the'
                 u'measurements of a selectively monitored reaction at each'
                 u'time point.'),
                'data file content',
                [u'mass chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'selected ion monitoring chromatogram', u'MS:1001472',
                (u'Chromatogram created by creating an array of the'
                 u'measurements of a selectively monitored ion at each time'
                 u'point.'),
                'data file content',
                [u'mass chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'consecutive reaction monitoring chromatogram', u'MS:1001474',
                (u'OBSOLETE Chromatogram created by creating an array of the'
                 u'measurements of a series of monitored reactions at each time'
                 u'point.'),
                'data file content',
                [u'mass chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'basepeak chromatogram', u'MS:1000628',
                (u'Chromatogram created by creating an array of the most'
                 u'intense peaks at each time point.'),
                'data file content',
                [u'mass chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'selected ion current chromatogram', u'MS:1000627',
                (u'Chromatogram created by creating an array of the'
                 u'measurements of a specific single ion current at each time'
                 u'point.'),
                'data file content',
                [u'mass chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'absorption chromatogram', u'MS:1000812',
                (u'The measurement of light absorbed by the sample as a'
                 u'function of the retention time.'),
                'data file content',
                [u'electromagnetic radiation chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'emission chromatogram', u'MS:1000813',
                (u'The measurement of light emitted by the sample as a function'
                 u'of the retention time.'),
                'data file content',
                [u'electromagnetic radiation chromatogram', u'data file content', u'chromatogram type']),
    FileContent(u'time-delayed fragmentation spectrum', u'MS:1000790',
                (u'MSn spectrum in which the product ions are collected after a'
                 u'time delay, which allows the observation of lower energy'
                 u'fragmentation processes after precursor ion activation.'),
                'data file content',
                [u'MSn spectrum', u'mass spectrum', u'data file content', u'spectrum type']),
    FileContent(u'enhanced multiply charged spectrum', u'MS:1000789',
                (u'MS1 spectrum that is enriched in multiply-charged ions'
                 u'compared to singly-charged ions.'),
                'data file content',
                [u'MS1 spectrum', u'mass spectrum', u'data file content', u'spectrum type']),
])
# [[[end]]]


spectrum_representation = []
# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('spectrum representation',
#             "spectrum_representation", term_cls_name="FileContent", writer=cog.out)
# ]]]
spectrum_representation = TermSet([
    FileContent(u'centroid spectrum', u'MS:1000127',
                (u'Processing of profile data to produce spectra that contains'
                 u'discrete peaks of zero width. Often used to reduce the size'
                 u'of dataset.'),
                'spectrum representation',
                [u'spectrum representation']),
    FileContent(u'profile spectrum', u'MS:1000128',
                (u'A profile mass spectrum is created when data is recorded'
                 u'with ion current (counts per second) on one axis and'
                 u'mass/charge ratio on another axis.'),
                'spectrum representation',
                [u'spectrum representation']),
])
# [[[end]]]


content_keys = content_keys + spectrum_representation

id_formats_by_name = {k.name: k for k in id_formats}
file_formats_by_name = {k.name: k for k in file_formats}
content_keys_by_name = {k.name: k for k in content_keys}


MS_MS1_Spectrum = content_keys.get('MS1 spectrum')
MS_MSn_Spectrum = content_keys.get('MSn spectrum')


def id_format(name):
    '''Translate a given name or identifier into a :class:`IDFormat`
    instance.

    If no match is found in the database of known :class:`IDFormat`
    types, a new dummy :class:`IDFormat` is returned with all fields
    set to the value of ``name``

    Returns
    -------
    IDFormat
    '''
    try:
        return id_formats[name]
    except KeyError:
        if name is None:
            return None
        return IDFormat(name, name, name, name, [name])


def file_format(name):
    '''Translate a given name or identifier into a :class:`FileFormat`
    instance.

    If no match is found in the database of known :class:`FileFormat`
    types, a new dummy :class:`FileFormat` is returned with all fields
    set to the value of ``name``

    Returns
    -------
    FileFormat
    '''
    try:
        return file_formats[name]
    except KeyError:
        if name is None:
            return None
        return FileFormat(name, name, name, name, [name])


def content_key(name):
    '''Translate a given name or identifier into a :class:`FileContent`
    instance.

    If no match is found in the database of known :class:`FileContent`
    types, a new dummy :class:`FileContent` is returned with all fields
    set to the value of ``name``

    Returns
    -------
    FileContent
    '''
    try:
        return content_keys[name]
    except KeyError:
        if name is None:
            return None
        return FileContent(name, name, name, name, [name])


class FileInformation(MutableMapping):
    """Describes the type of data found in this file and the
    source files that contributed to it.

    Implements the :class:`MutableMapping` Interface

    Attributes
    ----------
    contents : dict
        A mapping between controlled vocabullary names or user-defined
        names and an optional value. For standard controlled names see
        :data:`content_keys`
    source_files : list of :class:`.SourceFile` objects
        The set of files which either define the current file, or were used
        to create the current file if recorded.
    """

    def __init__(self, contents=None, source_files=None):
        if contents is None:
            contents = {}
        if source_files is None:
            source_files = []
        self.contents = dict(contents)
        self.source_files = list(source_files)

    def add_file(self, source, check=True):
        """Add a new file to :attr:`source_files`

        If ``source`` is a string, it will be interpreted as a path
        and an instance of :class:`SourceFile` will be created using
        :meth:`SourceFile.from_path`. Otherwise, it is assumed to be
        an instance of :class:`SourceFile`.

        Parameters
        ----------
        source : str or :class:`SourceFile`
            Either the path to a file to be added to the source file
            collection, or an instance of :class:`SourceFile`
        check : bool, optional
            Whether or not to check and validate that a path points to
            a real file

        Raises
        ------
        ValueError
            If a path fails to validate as real
        """
        if isinstance(source, basestring):
            source = os.path.realpath(source)
            if check:
                if not os.path.exists(source):
                    raise ValueError(
                        "Source File %r does not exist" % (source,))
            source = SourceFile.from_path(source)
        elif not isinstance(source, SourceFile):
            raise TypeError("Must pass an object of type %r, could not coerce %r" % (SourceFile, type(source)))
        self.source_files.append(source)

    def add_content(self, key, value=None):
        """Adds a new key-value pair to :attr:`contents` with
        an optional value

        Parameters
        ----------
        key : str or :attr:`content`
            The content name, either a CV-term or a user-defined name
        value : object, optional
            The optional value, which should be any type of object whose
            meaning makes sense given the definition of ``key``
        """
        self.contents[key] = value

    def remove_content(self, key):
        """Remove a key from :attr:`content`

        Parameters
        ----------
        key : str or :class:`FileContent`
            The content key to remove
        """
        self.contents.pop(key, None)

    def get_content(self, key):
        '''Retrieve the value of ``key`` from :attr:`contents`.

        This method is aliased to :meth:`__getitem__`

        Parameters
        ----------
        key : str or :class:`FileContent`

        Returns
        -------
        object
        '''
        return self.contents[key]

    def has_content(self, key):
        '''Check if ``key`` is found in :attr:`content`

        Parameters
        ----------
        key : str or :class:`FileContent`

        Returns
        -------
        bool
        '''
        return key in self.contents

    def __getitem__(self, key):
        return self.get_content(key)

    def __contains__(self, key):
        return self.has_content(key)

    def __setitem__(self, key, value):
        self.add_content(key, value)

    def __delitem__(self, key):
        self.remove_content(key)

    def __repr__(self):
        template = "FileInformation(%s, %s)"
        return template % (self.contents, self.source_files)

    def __len__(self):
        return len(self.contents)

    def __iter__(self):
        return iter(self.contents)

    def copy(self):
        '''Create a deep copy of this object

        Returns
        -------
        FileInformation
        '''
        return self.__class__(
            self.contents.copy(), [f.copy() for f in self.source_filesr])


format_parameter_map = {
    "thermo raw": (id_formats_by_name.get("Thermo nativeID format"),
                   file_formats_by_name.get("Thermo RAW format")),
    "agilent d": (id_formats_by_name.get("Agilent MassHunter nativeID format"),
                  file_formats_by_name.get("Agilent MassHunter format")),
    'mgf': (id_formats_by_name.get("no nativeID format"),
            # id_formats_by_name.get("multiple peak list nativeID format"),
            file_formats_by_name.get('Mascot MGF format')),
}


class SourceFile(object):
    """Represents a single raw data file which either defines or contributed data to another
    data file, the "reference file"

    Attributes
    ----------
    file_format : :class:`~.FileFormat`
        The name of a data file format. See :data:`file_formats`
    id : str
        The unique identifier for this file, among files which contributed to
        the reference file
    id_format : :class:`~.IDFormat`
        The name of a formal identifier schema. See :data:`~.id_formats`
    location : str
        The directory path to this file on the machine it was last read on to contribute to
        or define the reference file
    name : str
        The base name of this file
    parameters : dict
        A set of key-value pairs associated with this file, either encoding extra metadata
        annotations, or precomputed hash checksums
    """

    _checksum_translation_map = {
        'sha1': 'sha1',
        'SHA-1': 'sha1',
        'md5': 'md5',
        'MD5': 'md5'
    }

    @classmethod
    def _resolve_checksum_hash_type(cls, hash_type):
        try:
            return cls._checksum_translation_map[hash_type]
        except KeyError:
            try:
                return cls._checksum_translation_map[hash_type.lower().replace("-", '')]
            except KeyError:
                raise KeyError(hash_type)

    @classmethod
    def from_path(cls, path):
        '''Construct a new :class:`SourceFile` from a path to a real file on the local
        file system.

        Parameters
        ----------
        path: str
            The path to the file to describe

        Returns
        -------
        SourceFile
        '''
        path = os.path.realpath(path)
        name = os.path.basename(path)
        location = os.path.dirname(path)
        idfmt, file_fmt = cls.guess_format(path)
        source = cls(name, location, name, idfmt, file_fmt)
        return source

    def __init__(self, name, location, id=None, id_format=None, file_format=None, parameters=None):
        if id is None:
            id = name
        if parameters is None:
            parameters = {}
        self.name = name
        self.location = location.replace("file:///", '')
        self.id = id
        self.parameters = dict(parameters)
        self._clean_parameters()
        self.id_format = id_format or self.infer_id_format_from_paramters()
        self.file_format = file_format or self.infer_file_format_from_parameters()

    @property
    def id_format(self):
        return self._id_format

    @id_format.setter
    def id_format(self, value):
        if value is None:
            self._id_format = None
        else:
            self._id_format = id_format(str(value))

    @property
    def file_format(self):
        return self._file_format

    @file_format.setter
    def file_format(self, value):
        if value is None:
            self._file_format = None
        else:
            self._file_format = file_format(str(value))

    def _clean_parameters(self):
        self.parameters.pop("location", None)
        self.parameters.pop("id", None)
        self.parameters.pop("name", None)

    def infer_id_format_from_paramters(self):
        try:
            fmt = list(set(self.parameters) & set(id_formats))[0]
            self.parameters.pop(fmt)
            return fmt
        except IndexError:
            return None

    def infer_file_format_from_parameters(self):
        try:
            fmt = list(set(self.parameters) & set(file_formats))[0]
            self.parameters.pop(fmt)
            return fmt
        except IndexError:
            return None

    @property
    def path(self):
        return os.path.join(self.location, self.name)

    def is_resolvable(self):
        return os.path.exists(self.path)

    @staticmethod
    def guess_format(path):
        if not os.path.exists(path):
            return None, None
        if os.path.isdir(path):
            if os.path.exists(os.path.join(path, 'AcqData')):
                return format_parameter_map['agilent d']

        parts = os.path.splitext(path)
        if len(parts) > 1:
            is_compressed = False
            ext = parts[1]
            if ext.lower() == '.gz':
                is_compressed = True
                parts = os.path.splitext(parts[0])
                ext = parts[1]
            if ext.lower() == '.mzml':
                fmt = "mzML format"
                id_fmt = "no nativeID format"
                hit = False
                if is_compressed:
                    from .._compression import get_opener
                    fh = get_opener(path)
                else:
                    fh = open(path, 'rb')
                with fh:
                    from ..xml_reader import iterparse_until
                    for sf_tag in iterparse_until(fh, 'sourceFile', 'run'):
                        for param in sf_tag.getchildren():
                            if "nativeID" in param.attrib['name']:
                                id_fmt = param.attrib['name']
                                hit = True
                                break
                        if hit:
                            break
                return id_fmt, fmt
            elif ext.lower() == '.mzxml':
                fmt = "ISB mzXML format"
                id_fmt = "no nativeID format"
                return fmt, id_fmt
            elif ext.lower() == '.mgf':
                fmt = "Mascot MGF format"
                id_fmt = "no nativeID format"
                return fmt, id_fmt
        with open(path, 'rb') as fh:
            lead_bytes = fh.read(30)
            # looking for pattern matching b'\x01\xa1F\x00i\x00n\x00n\x00i\x00g\x00a\x00n\x00'
            decoded = lead_bytes.decode("utf-16")[1:9]
            if decoded == "Finnigan":
                return format_parameter_map['thermo raw']
        return None, None

    def __repr__(self):
        template = "SourceFile(%r, %r, %r, %s, %s%s)"
        if self.parameters:
            tail = ", %r" % self.parameters
        else:
            tail = ''
        return template % (
            self.name, self.location,
            self.id, self.id_format, self.file_format,
            tail
        )

    def _compute_checksum(self, hash_type='sha1', buffer_size=2**16):
        from .._compression import get_opener
        hasher = hashlib.new(hash_type)
        buffer_size = int(buffer_size)
        with get_opener(self.path, 'rb') as fh:
            content_buffer = fh.read(buffer_size)
            while content_buffer:
                hasher.update(content_buffer)
                content_buffer = fh.read(buffer_size)
        return hasher.hexdigest()

    def checksum(self, hash_type='sha1'):
        hash_type = self._resolve_checksum_hash_type(hash_type)
        return self._compute_checksum(hash_type)

    def add_checksum(self, hash_type='sha1'):
        hash_type = self._resolve_checksum_hash_type(hash_type)
        checksum = self.checksum(hash_type)
        if hash_type == 'sha1':
            self.parameters['SHA-1'] = checksum
        elif hash_type == "md5":
            self.parameters['MD5'] = checksum

    def validate_checksum(self):
        if not os.path.exists(self.path):
            FileNotFoundError("%s not found" % (self.path,))
        if 'SHA-1' in self.parameters:
            checksum = self.checksum('sha1')
            return self.parameters['SHA-1'] == checksum
        elif 'MD5' in self.parameters:
            checksum = self.checksum("md5")
            return self.parameters['MD5'] == checksum
        else:
            warnings.warn(
                "%r did not have a reference checksum. Could not validate" % (self,))
            return True

    def has_checksum(self, hash_type=None):
        if hash_type is None:
            return ("SHA-1" in self.parameters) or ("MD5" in self.parameters)
        elif self._resolve_checksum_hash_type(hash_type) == 'sha1':
            return ("SHA-1" in self.parameters)
        elif self._resolve_checksum_hash_type(hash_type) == 'md5':
            return "MD5" in self.parameters

    def __eq__(self, other):
        if other is None:
            return False
        if self.path == other.path:
            if self.is_resolvable() and other.is_resolvable():
                return self.checksum() == other.checksum()
            else:
                if self.is_resolvable():
                    for hash_type in ['SHA-1', 'MD5']:
                        if other.has_checksum(hash_type):
                            return self.checksum(hash_type) == other.parameters[hash_type]
                elif other.is_resolvable():
                    for hash_type in ['SHA-1', 'MD5']:
                        if self.has_checksum(hash_type):
                            return other.checksum(hash_type) == self.parameters[hash_type]
                else:
                    for hash_type in ['SHA-1', 'MD5']:
                        if self.has_checksum(hash_type) and other.has_checksum(hash_type):
                            return other.parameters[hash_type] == self.parameters[hash_type]
        return False

    def __ne__(self, other):
        return not self == other

    def copy(self):
        return self.__class__(self.name, self.location, self.id,
                              self.id_format, self.file_format,
                              parameters=self.parameters.copy())


__all__ = [
    "IDFormat", "FileFormat", "FileContent",
    "id_formats", "file_formats", "content_keys",
    "id_format", "file_format", "content_key",
    "FileInformation", "SourceFile"
]
