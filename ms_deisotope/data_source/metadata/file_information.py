"""Defines types for describing different kinds of mass spectrometry
data files and their contents, and a database of controlled vocabulary
terms for them.
"""

import os
import re
import hashlib
import warnings

from collections import OrderedDict

try:
    from collections.abc import MutableMapping, Mapping
except ImportError:
    from collections import MutableMapping, Mapping

from six import string_types as basestring

from .cv import Term, TermSet

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError


class IDParserBase(object):
    """A base class for ID parsing and formatting.
    """
    def parse(self, text):
        """Parse a string looking for fields defined by this ID format.

        Parameters
        ----------
        text : str
            The string to parse

        Returns
        -------
        dict:
            The parsed fields of the ID string.
        """
        raise NotImplementedError()

    def format(self, fields):
        """Format a set of fields as a nativeID string.

        Parameters
        ----------
        fields : dict
            The fields to populate the nativeID from.

        Returns
        -------
        str
        """
        raise NotImplementedError()

    def __call__(self, text):
        """Parse a string looking for fields defined by this ID format.

        Parameters
        ----------
        text : str
            The string to parse

        Returns
        -------
        dict:
            The parsed fields of the ID string.

        See Also
        --------
        :meth:`parse`
        """
        return self.parse(text)


class IDFormat(Term, IDParserBase):
    """Describes a named spectrum identifier format, either
    using a controlled-vocabulary term or user-defined name.

    A :class:`IDFormat` is equal to its name and its controlled
    vocabulary identifier.

    An instance may also be used to parse a string in its format
    using its :meth:`parse` method, creating a :class:`dict` of
    its fields.

    Attributes
    ----------
    parser: NativeIDParser
        A parser for the specified nativeID format type.
    """

    def __init__(self, name, id, description, category, specialization):
        super(IDFormat, self).__init__(
            name, id, description, category, specialization)
        self.parser = NativeIDParser.from_term(self)

    def parse(self, text):
        """Parse a string looking for fields defined by this term's nativeID
        format

        Parameters
        ----------
        text : str
            The string to parse

        Returns
        -------
        dict:
            The parsed fields of the ID string.
        """
        if self.parser is None:
            raise ValueError(
                "This IDFormat doesn't have a pattern!")
        return self.parser.parse(text)

    def format(self, fields):
        """Format a set of fields as a nativeID string.

        Parameters
        ----------
        fields : dict
            The fields to populate the nativeID from.

        Returns
        -------
        str
        """
        return self.parser.format(fields)

    def __call__(self, text):
        return self.parse(text)

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


type_pat = re.compile("([A-Za-z]+)=xsd:(%s+)" % '|'.join(
    {'IDREF', "long", 'nonNegativeInteger', 'positiveInteger', 'string'}))

xsd_to_regex = {
    "IDREF": r"(\S+)",
    "long": r"(-?\d+)",
    "nonNegativeInteger": r"(\d+)",
    "positiveInteger": r"(\d+)",
    "string": r"(\S+)",
}

xsd_to_type = {
    "IDREF": str,
    "long": int,
    "nonNegativeInteger": int,
    "positiveInteger": int,
    "string": str,
}


class NativeIDParser(IDParserBase):
    """A parser for a single nativeID format.

    These may be automatically derived from the CV-param defining them by parsing the
    XSD string included, but no guarantee is available.
    """

    def __init__(self, parser, tokens, name):
        self.parser = parser
        self.tokens = OrderedDict(tokens)
        self.name = name

    @classmethod
    def from_term(cls, term):
        """Construct a :class:`NativeIDParser` from :class:`IDFormat` term.

        Parameters
        ----------
        term : IDFormat
            The nativeID format specification to build a parser for

        Returns
        -------
        :class:`NativeIDParser`:
            The constructed parser, or :const:`None` if no regular expression could be
            constructed.
        """
        if "Native format defined by" in term.description:
            tokens = []
            desc = term.description.split(
                "Native format defined by", 1)[1].rstrip()
            for mat in type_pat.finditer(desc):
                tokens.append(mat.groups())
            parser = re.compile(
                ''.join([r"(%s)=%s\s?" % (k, xsd_to_regex[v]) for k, v in tokens]))
            return cls(parser, tokens, term.name)
        return None

    def parse(self, string):
        """Parse a string according to this parser's pattern,
        returning the type-cast fields as a :class:`dict`.

        Parameters
        ----------
        string : str
            The string to parse

        Returns
        -------
        dict
            The fields of the scan ID

        Raises
        ------
        ValueError:
            If the string does not conform to the expected pattern
        """
        match = self.parser.match(string)
        if match is None:
            return OrderedDict()
        groups = match.groups()
        n = len(groups)
        i = 0
        fields = OrderedDict()
        while i < n:
            k = groups[i]
            v = groups[i + 1]
            i += 2
            try:
                v = int(v)
            except ValueError:
                pass
            fields[k] = v
        return fields

    def format(self, fields):
        """Format a set of fields as a nativeID string.

        Parameters
        ----------
        fields : dict
            The fields to populate the nativeID from.

        Returns
        -------
        str
        """
        parts = []
        for key in self.tokens:
            parts.append("%s=%s" % (key, fields[key]))
        return ' '.join(parts)


class MultipleIDFormats(Mapping, IDParserBase):
    '''Represent an ambiguous group of multiple nativeID formats.

    Implements the :class:`~collections.abc.Mapping` interface.

    Attributes
    ----------
    id_formats : OrderedDict
        A mapping of format name to :class:`IDFormat` instances
    '''
    def __init__(self, id_formats):
        self.id_formats = id_formats

    def parse(self, text):
        fields = OrderedDict()
        for name, parser in self.id_formats.items():
            fields = parser.parse(text)
            if not fields:
                continue
            else:
                fields['id_format'] = name
                break
        return fields

    def format(self, fields):
        format_name = fields.get('id_format')
        id_format = self.id_formats[format_name]
        return id_format.format(fields)

    def __iter__(self):
        return iter(self.id_formats)

    def __getitem__(self, key):
        return self.id_formats[key]

    def __len__(self):
        return len(self.id_formats)

    def __repr__(self):
        template = "{self.__class__.__name__}({self.id_formats})"
        return template.format(self=self)


id_formats = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('native spectrum identifier format',
#             "id_formats", term_cls_name="IDFormat", writer=cog.out)
# ]]]
# CV Version: 4.1.55
id_formats = TermSet([
    IDFormat('Thermo nativeID format', 'MS:1000768',
             ('Native format defined by '
              'controllerType=xsd:nonNegativeInteger '
              'controllerNumber=xsd:positiveInteger '
              'scan=xsd:positiveInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Waters nativeID format', 'MS:1000769',
             ('Native format defined by function=xsd:positiveInteger '
              'process=xsd:nonNegativeInteger scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('WIFF nativeID format', 'MS:1000770',
             ('Native format defined by sample=xsd:nonNegativeInteger '
              'period=xsd:nonNegativeInteger cycle=xsd:nonNegativeInteger '
              'experiment=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Bruker/Agilent YEP nativeID format', 'MS:1000771',
             ('Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Bruker BAF nativeID format', 'MS:1000772',
             ('Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Bruker FID nativeID format', 'MS:1000773',
             ('Native format defined by file=xsd:IDREF.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('multiple peak list nativeID format', 'MS:1000774',
             ('Native format defined by index=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('single peak list nativeID format', 'MS:1000775',
             ('Native format defined by file=xsd:IDREF.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('scan number only nativeID format', 'MS:1000776',
             ('Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('spectrum identifier nativeID format', 'MS:1000777',
             ('Native format defined by spectrum=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Bruker U2 nativeID format', 'MS:1000823',
             ('Native format defined by declaration=xsd:nonNegativeInteger '
              'collection=xsd:nonNegativeInteger '
              'scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('no nativeID format', 'MS:1000824',
             ('No nativeID format indicates that the file tagged with this '
              'term does not contain spectra that can have a nativeID '
              'format.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Shimadzu Biotech nativeID format', 'MS:1000929',
             ('Native format defined by source=xsd:string '
              'start=xsd:nonNegativeInteger end=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('SCIEX TOF/TOF nativeID format', 'MS:1001480',
             ('Native format defined by jobRun=xsd:nonNegativeInteger '
              'spotLabel=xsd:string spectrum=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Agilent MassHunter nativeID format', 'MS:1001508',
             ('Native format defined by scanId=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('spectrum from database integer nativeID format', 'MS:1001526',
             ('Native format defined by databasekey=xsd:long.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Mascot query number', 'MS:1001528',
             ('Native format defined by query=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format', 'spectrum identification result details']),
    IDFormat('spectrum from ProteinScape database nativeID format', 'MS:1001531',
             ('Native format defined by databasekey=xsd:long.'),
             'native spectrum identifier format',
             ['native spectrum identifier format', 'spectra data details', 'search input details']),
    IDFormat('spectrum from database string nativeID format', 'MS:1001532',
             ('Native format defined by databasekey=xsd:string.'),
             'native spectrum identifier format',
             ['native spectrum identifier format', 'spectra data details', 'search input details']),
    IDFormat('SCIEX TOF/TOF T2D nativeID format', 'MS:1001559',
             ('Native format defined by file=xsd:IDREF.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Scaffold nativeID format', 'MS:1001562',
             ('Scaffold native ID format.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Bruker Container nativeID format', 'MS:1002303',
             ('Native identifier (UUID).'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('UIMF nativeID format', 'MS:1002532',
             ('Native format defined by frame=xsd:nonNegativeInteger '
              'scan=xsd:nonNegativeInteger '
              'frameType=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Bruker TDF nativeID format', 'MS:1002818',
             ('Native format defined by frame=xsd:nonNegativeInteger '
              'scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
    IDFormat('Shimadzu Biotech QTOF nativeID format', 'MS:1002898',
             ('Native format defined by scan=xsd:nonNegativeInteger.'),
             'native spectrum identifier format',
             ['native spectrum identifier format']),
])
# [[[end]]]


file_formats = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('mass spectrometer file format',
#             "file_formats", term_cls_name="FileFormat", writer=cog.out)
# ]]]
# CV Version: 4.1.55
file_formats = TermSet([
    FileFormat('Waters raw format', 'MS:1000526',
               ('Waters data file format found in a Waters RAW directory, '
                'generated from an MS acquisition.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('ABI WIFF format', 'MS:1000562',
               ('Applied Biosystems WIFF file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Thermo RAW format', 'MS:1000563',
               ('Thermo Scientific RAW file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('PSI mzData format', 'MS:1000564',
               ('Proteomics Standards Inititative mzData file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Micromass PKL format', 'MS:1000565',
               ('Micromass PKL file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('ISB mzXML format', 'MS:1000566',
               ('Institute of Systems Biology mzXML file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Bruker/Agilent YEP format', 'MS:1000567',
               ('Bruker/Agilent YEP file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('mzML format', 'MS:1000584',
               ('Proteomics Standards Inititative mzML file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('DTA format', 'MS:1000613',
               ('SEQUEST DTA file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('ProteinLynx Global Server mass spectrum XML format', 'MS:1000614',
               ('Peak list file format used by ProteinLynx Global Server.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('parameter file', 'MS:1000740',
               ('Parameter file used to configure the acquisition of raw data '
                'on the instrument.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Bioworks SRF format', 'MS:1000742',
               ('Thermo Finnigan SRF file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'intermediate analysis format', 'file format']),
    FileFormat('Bruker BAF format', 'MS:1000815',
               ('Bruker BAF raw file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Bruker U2 format', 'MS:1000816',
               ('Bruker HyStar U2 file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Bruker FID format', 'MS:1000825',
               ('Bruker FID file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Shimadzu Biotech database entity', 'MS:1000930',
               ('Shimadzu Biotech format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Mascot MGF format', 'MS:1001062',
               ('Mascot MGF file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('PerSeptive PKS format', 'MS:1001245',
               ('PerSeptive peak list file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('SCIEX API III format', 'MS:1001246',
               ('PE SCIEX peak list file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Bruker XML format', 'MS:1001247',
               ('Bruker data exchange XML format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('text format', 'MS:1001369',
               ('Simple text file format of \\"m/z [intensity]\\" values for a '
                'PMF (or single MS2) search.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Phenyx XML format', 'MS:1001463',
               ('Phenyx open XML file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'intermediate analysis format', 'file format']),
    FileFormat('MS2 format', 'MS:1001466',
               ('MS2 file format for MS2 spectral data." [PMID:15317041, '
                'DOI:10.1002/rcm.1603'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('SCIEX TOF/TOF database', 'MS:1001481',
               ('Applied Biosystems/MDS Analytical Technologies TOF/TOF '
                'instrument database.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Agilent MassHunter format', 'MS:1001509',
               ('A data file format found in an Agilent MassHunter directory '
                'which contains raw data acquired by an Agilent mass '
                'spectrometer.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Proteinscape spectra', 'MS:1001527',
               ('Spectra from Bruker/Protagen Proteinscape database.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('SCIEX TOF/TOF T2D format', 'MS:1001560',
               ('Applied Biosystems/MDS Analytical Technologies TOF/TOF '
                'instrument export format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('mz5 format', 'MS:1001881',
               ('mz5 file format, modelled after mzML.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Bruker Container format', 'MS:1002302',
               ('Bruker Container raw file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('SCiLS Lab format', 'MS:1002385',
               ('SCiLS Lab file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Andi-MS format', 'MS:1002441',
               ('AIA Analytical Data Interchange file format for mass '
                'spectrometry data.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('UIMF format', 'MS:1002531',
               ('SQLite-based file format created at Pacific Northwest '
                'National Lab. It stores an intermediate analysis of ion- '
                'mobility mass spectrometry data.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('MS1 format', 'MS:1002597',
               ('MS1 file format for MS1 spectral data." [PMID:15317041'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Bruker TDF format', 'MS:1002817',
               ('Bruker TDF raw file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('mzMLb format', 'MS:1002838',
               ('mzMLb file format, mzML encapsulated within HDF5." [PSI:PI'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('msalign format', 'MS:1002899',
               ('msalign file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('feature format', 'MS:1002900',
               ('TopFD feature file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('chrom format', 'MS:1002966',
               ('The Lipid Data Analyzer native chrom format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Andromeda:apl file format', 'MS:1002996',
               ('Peak list file format of the Andromeda search engine.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
    FileFormat('Shimadzu Biotech LCD format', 'MS:1003009',
               ('Shimadzu Biotech LCD file format.'),
               'mass spectrometer file format',
               ['mass spectrometer file format', 'file format']),
])
# [[[end]]]


content_keys = []

# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('data file content',
#             "content_keys", term_cls_name="FileContent", writer=cog.out)
# ]]]
# CV Version: 4.1.55
content_keys = TermSet([
    FileContent('mass spectrum', 'MS:1000294',
                ('A plot of the relative abundance of a beam or other '
                 'collection of ions as a function of the mass-to-charge ratio '
                 '(m/z).'),
                'data file content',
                ['data file content', 'spectrum type']),
    FileContent('PDA spectrum', 'MS:1000620',
                ('OBSOLETE Spectrum generated from a photodiode array detector '
                 '(ultraviolet/visible spectrum).'),
                'data file content',
                ['data file content', 'spectrum type']),
    FileContent('electromagnetic radiation spectrum', 'MS:1000804',
                ('A plot of the relative intensity of electromagnetic '
                 'radiation as a function of the wavelength.'),
                'data file content',
                ['data file content', 'spectrum type']),
    FileContent('emission spectrum', 'MS:1000805',
                ('A plot of the relative intensity of electromagnetic '
                 'radiation emitted by atoms or molecules when excited.'),
                'data file content',
                ['data file content', 'spectrum type']),
    FileContent('absorption spectrum', 'MS:1000806',
                ('A plot of the relative intensity of electromagnetic '
                 'radiation absorbed by atoms or molecules when excited.'),
                'data file content',
                ['data file content', 'spectrum type']),
    FileContent('ion current chromatogram', 'MS:1000810',
                ('Representation of the current of ions versus time.'),
                'data file content',
                ['data file content', 'chromatogram type']),
    FileContent('electromagnetic radiation chromatogram', 'MS:1000811',
                ('Representation of electromagnetic properties versus time.'),
                'data file content',
                ['data file content', 'chromatogram type']),
    FileContent('charge inversion mass spectrum', 'MS:1000322',
                ('The measurement of the relative abundance of ions that '
                 'result from a charge inversion reaction as a function of '
                 'm/z.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('constant neutral gain spectrum', 'MS:1000325',
                ('A spectrum formed of all product ions that have been '
                 'produced by gain of a pre-selected neutral mass following '
                 'the reaction with and addition of the gas in a collision '
                 'cell.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('constant neutral loss spectrum', 'MS:1000326',
                ('A spectrum formed of all product ions that have been '
                 'produced with a selected m/z decrement from any precursor '
                 'ions. The spectrum shown correlates to the precursor ion '
                 'spectrum. See also neutral loss spectrum.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('e/2 mass spectrum', 'MS:1000328',
                ('A mass spectrum obtained using a sector mass spectrometer in '
                 'which the electric sector field E is set to half the value '
                 'required to transmit the main ion-beam. This spectrum '
                 'records the signal from doubly charged product ions of '
                 'charge-stripping reactions.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('precursor ion spectrum', 'MS:1000341',
                ('Spectrum generated by scanning precursor m/z while '
                 'monitoring a fixed product m/z.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('product ion spectrum', 'MS:1000343',
                ('OBSOLETE A mass spectrum recorded from any spectrometer in '
                 'which the appropriate m/z separation scan function is set to '
                 'record the product ion or ions of selected precursor ions.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('MS1 spectrum', 'MS:1000579',
                ('Mass spectrum created by a single-stage MS experiment or the '
                 'first stage of a multi-stage experiment.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('MSn spectrum', 'MS:1000580',
                ('MSn refers to multi-stage MS2 experiments designed to record '
                 'product ion spectra where n is the number of product ion '
                 'stages (progeny ions). For ion traps, sequential MS/MS '
                 'experiments can be undertaken where n > 2 whereas for a '
                 'simple triple quadrupole system n=2. Use the term ms level '
                 '(MS:1000511) for specifying n.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('CRM spectrum', 'MS:1000581',
                ('Spectrum generated from MSn experiment with three or more '
                 'stages of m/z separation and in which a particular multi- '
                 'step reaction path is monitored.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('SIM spectrum', 'MS:1000582',
                ('Spectrum obtained with the operation of a mass spectrometer '
                 'in which the abundances of one ion or several ions of '
                 'specific m/z values are recorded rather than the entire mass '
                 'spectrum (Selected Ion Monitoring).'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('SRM spectrum', 'MS:1000583',
                ('Spectrum obtained when data are acquired from specific '
                 'product ions corresponding to m/z values of selected '
                 'precursor ions a recorded via two or more stages of mass '
                 'spectrometry. The precursor/product ion pair is called a '
                 'transition pair. Data can be obtained for a single '
                 'transition pair or multiple transition pairs. Multiple time '
                 'segments of different transition pairs can exist in a single '
                 'file. Single precursor ions can have multiple product ions '
                 'consitituting multiple transition pairs. Selected reaction '
                 'monitoring can be performed as tandem mass spectrometry in '
                 'time or tandem mass spectrometry in space.'),
                'data file content',
                ['mass spectrum', 'data file content', 'spectrum type']),
    FileContent('total ion current chromatogram', 'MS:1000235',
                ('Representation of the total ion current detected in each of '
                 'a series of mass spectra versus time.'),
                'data file content',
                ['ion current chromatogram', 'data file content', 'chromatogram type']),
    FileContent('selected ion current chromatogram', 'MS:1000627',
                ('Representation of an array of the measurements of a specific '
                 'single ion current versus time.'),
                'data file content',
                ['ion current chromatogram', 'data file content', 'chromatogram type']),
    FileContent('basepeak chromatogram', 'MS:1000628',
                ('Representation of an array of the most intense peaks versus '
                 'time.'),
                'data file content',
                ['ion current chromatogram', 'data file content', 'chromatogram type']),
    FileContent('selected ion monitoring chromatogram', 'MS:1001472',
                ('Representation of an array of the measurements of a '
                 'selectively monitored ion versus time.'),
                'data file content',
                ['ion current chromatogram', 'data file content', 'chromatogram type']),
    FileContent('selected reaction monitoring chromatogram', 'MS:1001473',
                ('Representation of an array of the measurements of a '
                 'selectively monitored reaction versus time.'),
                'data file content',
                ['ion current chromatogram', 'data file content', 'chromatogram type']),
    FileContent('consecutive reaction monitoring chromatogram', 'MS:1001474',
                ('OBSOLETE Representation of an array of the measurements of a '
                 'series of monitored reactions versus time.'),
                'data file content',
                ['ion current chromatogram', 'data file content', 'chromatogram type']),
    FileContent('absorption chromatogram', 'MS:1000812',
                ('Representation of light absorbed by the sample versus time.'),
                'data file content',
                ['electromagnetic radiation chromatogram', 'data file content', 'chromatogram type']),
    FileContent('emission chromatogram', 'MS:1000813',
                ('Representation of light emitted by the sample versus time.'),
                'data file content',
                ['electromagnetic radiation chromatogram', 'data file content', 'chromatogram type']),
    FileContent('enhanced multiply charged spectrum', 'MS:1000789',
                ('MS1 spectrum that is enriched in multiply-charged ions '
                 'compared to singly-charged ions.'),
                'data file content',
                ['MS1 spectrum', 'mass spectrum', 'data file content', 'spectrum type']),
    FileContent('time-delayed fragmentation spectrum', 'MS:1000790',
                ('MSn spectrum in which the product ions are collected after a '
                 'time delay, which allows the observation of lower energy '
                 'fragmentation processes after precursor ion activation.'),
                'data file content',
                ['MSn spectrum', 'mass spectrum', 'data file content', 'spectrum type']),
])
# [[[end]]]


spectrum_representation = []
# [[[cog
# import cog
# from ms_deisotope.data_source.metadata.cv import render_list
# render_list('spectrum representation',
#             "spectrum_representation", term_cls_name="FileContent", writer=cog.out)
# ]]]
# CV Version: 4.1.55
spectrum_representation = TermSet([
    FileContent('centroid spectrum', 'MS:1000127',
                ('Processing of profile data to produce spectra that contains '
                 'discrete peaks of zero width. Often used to reduce the size '
                 'of dataset.'),
                'spectrum representation',
                ['spectrum representation']),
    FileContent('profile spectrum', 'MS:1000128',
                ('A profile mass spectrum is created when data is recorded '
                 'with ion current (counts per second) on one axis and '
                 'mass/charge ratio on another axis.'),
                'spectrum representation',
                ['spectrum representation']),
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
        self._id_format = None

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
            raise TypeError("Must pass an object of type %r, could not coerce %r" % (
                SourceFile, type(source)))
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
            self.contents.copy(), [f.copy() for f in self.source_files])

    def _find_native_id_format(self):
        options = list()
        unique = set()
        sf = None  # type: SourceFile
        for sf in self.source_files:
            i = sf.id_format
            if i is None:
                continue
            elif i not in unique:
                unique.add(i)
                options.append(i)
        n = len(options)
        if n == 0:
            return id_format("MS:1000824")
        elif n == 1:
            return options[0]
        else:
            no_format = id_format("MS:1000824")
            if no_format in options:
                options.remove(no_format)
                n = len(options)
            if n == 1:
                return options[0]
            else:
                return MultipleIDFormats(OrderedDict([(fmt.name, fmt) for fmt in options]))

    @property
    def id_format(self):
        if self._id_format is None:
            self._id_format = self._find_native_id_format()
        return self._id_format


format_parameter_map = {
    "thermo raw": (id_formats_by_name.get("Thermo nativeID format"),
                   file_formats_by_name.get("Thermo RAW format")),
    "agilent d": (id_formats_by_name.get("Agilent MassHunter nativeID format"),
                  file_formats_by_name.get("Agilent MassHunter format")),
    'mgf': (id_formats_by_name.get("no nativeID format"),
            file_formats_by_name.get('Mascot MGF format')),
    'waters': (id_formats_by_name.get('Waters nativeID format'),
               file_formats_by_name.get("Waters raw format")),
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

        id_fmt = "no nativeID format"
        if os.path.isdir(path):
            if os.path.exists(os.path.join(path, 'AcqData')):
                return format_parameter_map['agilent d']
            elif os.path.basename(path).replace(os.sep, '').endswith("raw"):
                return format_parameter_map['waters']
            else:
                warnings.warn("Could not determine source file type from directory path %r" % (path, ))
                return id_fmt, None

        parts = os.path.splitext(path)
        if len(parts) > 1:
            is_compressed = False
            ext = parts[1]
            if ext.lower() == '.gz':
                is_compressed = True
                parts = os.path.splitext(parts[0])
                ext = parts[1]
            if ext.lower() == '.mzml':
                fmt = file_formats['MS:1000584']
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
                id_fmt = "scan number only nativeID format"
                return id_fmt, fmt
            elif ext.lower() == '.mgf':
                fmt = file_formats['MS:1001062']
                id_fmt = "no nativeID format"
                return id_fmt, fmt
            elif ext.lower() == '.mzmlb':
                fmt = file_formats['MS:1002838']
                id_fmt = "no nativeID format"
                from ms_deisotope.data_source.mzmlb import determine_if_available, MzMLbLoader
                # TODO: Try to open the file and get the nativeID format information from
                # the XML buffer, either just opening the file fully or by openinhg it at
                # the HDF5 level and grab the first 4k characters from the XML buffer and
                #  doing the same as the mzML case.
                if determine_if_available():
                    handle = MzMLbLoader(path, use_index=False)
                    fid = handle.file_description()
                    id_fmt = fid.id_format
                return id_fmt, fmt
        with open(path, 'rb') as fh:
            lead_bytes = fh.read(30)
            # looking for pattern matching b'\x01\xa1F\x00i\x00n\x00n\x00i\x00g\x00a\x00n\x00'
            decoded = lead_bytes.decode("utf-16")[1:9]
            if decoded == "Finnigan":
                return format_parameter_map['thermo raw']
        return id_format, None

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
        with get_opener(self.path) as fh:
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
