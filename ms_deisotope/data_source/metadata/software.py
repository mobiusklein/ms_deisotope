from .cv import Term, render_list


def __generate_list_code():
    '''Prints the code to generate these static lists
    '''
    render_list('software', list_name='software_names', term_cls_name="SoftwareName")


class SoftwareName(Term):
    pass


class Software(object):

    @classmethod
    def is_name(cls, name):
        try:
            software_names_by_name[name]
            return True
        except KeyError:
            return False

    def __init__(self, name=None, id=None, version=None, **options):
        if name is None:
            name, options = self._resolve_name_from_kwargs(options)
        if name is None and id is not None:
            name = id
        self.name = name
        self.id = id or name
        self.version = version
        self.options = options

    def _resolve_name_from_kwargs(self, options):
        names = dict()
        not_names = dict()
        for key, value in options.items():
            if self.is_name(key):
                names[key] = value
            else:
                not_names[key] = value
        options = not_names
        if len(names) == 1:
            name = list(names.keys())[0]
        elif len(names) == 0:
            name = None
        else:
            raise ValueError("Multiple possible names found")
        return name, options

    def __str__(self):
        return self.name

    def __repr__(self):
        template = '{self.__class__.__name__}({self.name}, {self.id}, {self.version})'
        return template.format(self=self)

    def __eq__(self, other):
        try:
            return self.name == other.name
        except AttributeError:
            return str(self) == str(other)


software_names = [
    SoftwareName(u'SCiLS software', u'MS:1002383', 'software', [u'software']),
    SoftwareName(u'acquisition software', u'MS:1001455',
                 'software', [u'software']),
    SoftwareName(u'data processing software',
                 u'MS:1001457', 'software', [u'software']),
    SoftwareName(u'analysis software', u'MS:1001456',
                 'software', [u'software']),
    SoftwareName(u'SCIEX software', u'MS:1000690', 'software', [u'software']),
    SoftwareName(u'Applied Biosystems software',
                 u'MS:1000691', 'software', [u'software']),
    SoftwareName(u'Bruker software', u'MS:1000692', 'software', [u'software']),
    SoftwareName(u'Thermo Finnigan software',
                 u'MS:1000693', 'software', [u'software']),
    SoftwareName(u'Waters software', u'MS:1000694', 'software', [u'software']),
    SoftwareName(u'SRM software', u'MS:1000871', 'software', [u'software']),
    SoftwareName(u'peptide attribute calculation software',
                 u'MS:1000873', 'software', [u'software']),
    SoftwareName(u'Agilent software', u'MS:1000689',
                 'software', [u'software']),
    SoftwareName(u'quantitation software name', u'MS:1001139',
                 'software', [u'software', u'quantification information']),
    SoftwareName(u'LECO software', u'MS:1001798', 'software', [u'software']),
    SoftwareName(u'custom unreleased software tool',
                 u'MS:1000799', 'software', [u'software']),
    SoftwareName(u'BSI software', u'MS:1001949', 'software', [u'software']),
    SoftwareName(u'Shimadzu Corporation software',
                 u'MS:1001557', 'software', [u'software']),
    SoftwareName(u'SCiLS Lab', u'MS:1002384', 'software', [
                 u'SCiLS software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Analyst', u'MS:1000551', 'software', [
                 u'SCIEX software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'apexControl', u'MS:1000706', 'software', [
                 u'Bruker software', u'acquisition software', u'software']),
    SoftwareName(u'MALDI Solutions LC-MALDI', u'MS:1002381', 'software',
                 [u'acquisition software', u'analysis software', u'data processing software', u'Shimadzu Corporation software', u'software']),
    SoftwareName(u'dpControl', u'MS:1000720', 'software', [
                 u'Bruker software', u'acquisition software', u'software']),
    SoftwareName(u'esquireControl', u'MS:1000721', 'software', [
                 u'Bruker software', u'acquisition software', u'software']),
    SoftwareName(u'HCTcontrol', u'MS:1000725', 'software', [
                 u'Bruker software', u'acquisition software', u'software']),
    SoftwareName(u'micrOTOFcontrol', u'MS:1000726', 'software', [
                 u'Bruker software', u'acquisition software', u'software']),
    SoftwareName(u'spControl', u'MS:1000737', 'software', [
                 u'Bruker software', u'acquisition software', u'software']),
    SoftwareName(u'ChromaTOF HRT software', u'MS:1001877', 'software', [
                 u'acquisition software', u'analysis software', u'data processing software', u'LECO software', u'software']),
    SoftwareName(u'MALDI Solutions Microbial Identification', u'MS:1001878', 'software', [
                 u'acquisition software', u'analysis software', u'data processing software', u'MALDI Solutions', u'software', u'Shimadzu Corporation software']),
    SoftwareName(u'6300 Series Ion Trap Data Analysis Software', u'MS:1000688', 'software', [
                 u'Agilent software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ChromaTOF software', u'MS:1001799', 'software', [
                 u'acquisition software', u'analysis software', u'data processing software', u'LECO software', u'software']),
    SoftwareName(u'MassHunter Data Acquisition', u'MS:1000678', 'software', [
                 u'Agilent software', u'acquisition software', u'software']),
    SoftwareName(u'MassHunter Easy Access', u'MS:1000679', 'software', [
                 u'Agilent software', u'acquisition software', u'software']),
    SoftwareName(u'GPS Explorer', u'MS:1000661', 'software', [
                 u'SCIEX software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Voyager Biospectrometry Workstation System', u'MS:1000539', 'software', [
                 u'Applied Biosystems software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Xcalibur', u'MS:1000532', 'software', [
                 u'Thermo Finnigan software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MassLynx', u'MS:1000534', 'software', [
                 u'Waters software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'4700 Explorer', u'MS:1000537', 'software', [
                 u'Applied Biosystems software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Data Explorer', u'MS:1000536', 'software', [
                 u'Applied Biosystems software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'4000 Series Explorer Software', u'MS:1000659', 'software', [
                 u'SCIEX software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'FlexControl', u'MS:1000540', 'software', [
                 u'Bruker software', u'acquisition software', u'software']),
    SoftwareName(u'SCIEX TOF/TOF Series Explorer Software', u'MS:1001483', 'software',
                 [u'SCIEX software', u'acquisition software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MALDI Solutions', u'MS:1001558', 'software', [
                 u'acquisition software', u'analysis software', u'data processing software', u'Shimadzu Corporation software', u'software']),
    SoftwareName(u'Trapper', u'MS:1000553', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'CLINPROT', u'MS:1000708', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'CLINPROT micro', u'MS:1000709', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'BioTools', u'MS:1000707', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'preprocessing software', u'MS:1002386', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'PIA', u'MS:1002387', 'software', [
                 u'postprocessing software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'PEAKS Online', u'MS:1001947', 'software', [
                 u'quantitation software name', u'analysis software', u'data processing software', u'software', u'quantification information']),
    SoftwareName(u'PEAKS Studio', u'MS:1001946', 'software', [
                 u'quantitation software name', u'analysis software', u'data processing software', u'software', u'quantification information']),
    SoftwareName(u'DataAnalysis', u'MS:1000719', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'CompassXtract', u'MS:1000718', 'software', [
                 u'Bruker software', u'data processing software', u'software']),
    SoftwareName(u'Compass OpenAccess', u'MS:1000715', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Compass for micrOTOF', u'MS:1000714', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'CompassXport', u'MS:1000717', 'software', [
                 u'Bruker software', u'data processing software', u'software']),
    SoftwareName(u'Compass for HCT/esquire', u'MS:1000713', 'software',
                 [u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Compass', u'MS:1000712', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'flexImaging', u'MS:1000722', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProfileAnalysis', u'MS:1000728', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MSDK', u'MS:1002645', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'QuantAnalysis', u'MS:1000736', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MzWiff', u'MS:1000591', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'SQID', u'MS:1001886', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Maltcms', u'MS:1002344', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MZmine', u'MS:1002342', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'PepFinder', u'MS:1002524', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'Spectrum Mill for MassHunter Workstation', u'MS:1000687', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'METLIN', u'MS:1000686', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MassHunter Mass Profiler', u'MS:1000685', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Genespring MS', u'MS:1000684', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MassHunter BioConfirm', u'MS:1000683', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MassHunter Metabolite ID', u'MS:1000682', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MassHunter Quantitative Analysis', u'MS:1000681', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MassHunter Qualitative Analysis', u'MS:1000680', 'software', [
                 u'Agilent software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'postprocessing software', u'MS:1002414', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'XCMS', u'MS:1001582', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'UNIFY', u'MS:1001796', 'software', [
                 u'Waters software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Empower', u'MS:1001795', 'software', [
                 u'Waters software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Pro Quant', u'MS:1000670', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Pro BLAST', u'MS:1000671', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MultiQuant', u'MS:1000674', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteinPilot Software', u'MS:1000663', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'LightSight Software', u'MS:1000662', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MarkerView Software', u'MS:1000665', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TissueView Software', u'MS:1000664', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'BioAnalyst', u'MS:1000667', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MRMPilot Software', u'MS:1000666', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Pro ICAT', u'MS:1000669', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Pro ID', u'MS:1000668', 'software', [
                 u'SCIEX software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'PRIDE Converter2', u'MS:1002335', 'software', [
                 u'conversion software', u'data processing software', u'software']),
    SoftwareName(u'ProCon', u'MS:1002334', 'software', [
                 u'conversion software', u'data processing software', u'software']),
    SoftwareName(u'conversion software', u'MS:1002333', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'PEAKS Node', u'MS:1001948', 'software', [
                 u'quantitation software name', u'analysis software', u'data processing software', u'software', u'quantification information']),
    SoftwareName(u'massWolf', u'MS:1000538', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'Bioworks', u'MS:1000533', 'software', [
                 u'Thermo Finnigan software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'FlexAnalysis', u'MS:1000535', 'software', [
                 u'Bruker software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Proteome Discoverer', u'MS:1000650', 'software', [
                 u'Thermo Finnigan software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MSnbase', u'MS:1002870', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'CAMERA', u'MS:1002871', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'mzR', u'MS:1002869', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'pymzML', u'MS:1001914', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'ReAdW', u'MS:1000541', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'MzStar', u'MS:1000542', 'software', [
                 u'data processing software', u'software']),
    SoftwareName(u'Maui', u'MS:1002452', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP software', u'MS:1000752', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteoWizard software', u'MS:1000615', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'PinPoint', u'MS:1001912', 'software', [
                 u'Thermo Finnigan software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteinLynx Global Server', u'MS:1000601', 'software', [
                 u'Waters software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Proteios', u'MS:1000600', 'software', [
                 u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Tide', u'MS:1002575', 'software', [
                 u'analysis software', u'software']),
    SoftwareName(u'Morpheus', u'MS:1002661', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Mascot Distiller', u'MS:1001488', 'software', [
                 u'quantitation software name', u'analysis software', u'software', u'quantification information']),
    SoftwareName(u'IsobariQ', u'MS:1002210', 'software', [
                 u'quantitation software name', u'analysis software', u'software', u'quantification information']),
    SoftwareName(u'Ascore software', u'MS:1001984', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'ProteinScape', u'MS:1000734', 'software', [
                 u'Bruker software', u'analysis software', u'software']),
    SoftwareName(u'greylag', u'MS:1001461', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Mascot Parser', u'MS:1001478', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'SpectraST', u'MS:1001477', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'X\\!Tandem', u'MS:1001476', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'OMSSA', u'MS:1001475', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'mzidLib', u'MS:1002237', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Pepitome', u'MS:1001588', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'MaxQuant', u'MS:1001583', 'software', [
                 u'quantitation software name', u'analysis software', u'software', u'quantification information']),
    SoftwareName(u'Comet', u'MS:1002251', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Andromeda', u'MS:1002337', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Amanda', u'MS:1002336', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'SEQUEST', u'MS:1001208', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Phenyx', u'MS:1001209', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Mascot', u'MS:1001207', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Percolator', u'MS:1001490', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'ProteinExtractor', u'MS:1001487', 'software', [
                 u'Bruker software', u'analysis software', u'software']),
    SoftwareName(u'Mascot Integra', u'MS:1001489', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Byonic', u'MS:1002261', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'PeptideShaker', u'MS:1002458', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'TagRecon', u'MS:1001587', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'DirecTag', u'MS:1001586', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'MyriMatch', u'MS:1001585', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'PAnalyzer', u'MS:1002076', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'NIST MSPepSearch', u'MS:1002750', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Trans-Proteomic Pipeline', u'MS:1002285',
                 'software', [u'analysis software', u'software']),
    SoftwareName(u'Trans-Proteomic Pipeline software', u'MS:1002286',
                 'software', [u'analysis software', u'software']),
    SoftwareName(u'DTASelect', u'MS:1002598', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'MSQuant', u'MS:1001977', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'DeBunker', u'MS:1001973', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Scaffold', u'MS:1001561', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'ProLuCID', u'MS:1002596', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'xiFDR', u'MS:1002543', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Skyline mzQuantML converter', u'MS:1002546', 'software', [
                 u'quantitation software name', u'analysis software', u'software', u'quantification information']),
    SoftwareName(u'xi', u'MS:1002544', 'software', [
                 u'analysis software', u'software']),
    SoftwareName(u'MSPathFinder', u'MS:1002720', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'MetaMorpheus', u'MS:1002826', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'MS-GF', u'MS:1002047', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'ProteinProspector', u'MS:1002043', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'MS-GF+', u'MS:1002048', 'software',
                 [u'analysis software', u'software']),
    SoftwareName(u'Cliquid', u'MS:1000672', 'software',
                 [u'SCIEX software', u'software']),
    SoftwareName(u'MIDAS Workflow Designer', u'MS:1000673',
                 'software', [u'SCIEX software', u'software']),
    SoftwareName(u'Compass Security Pack', u'MS:1000716',
                 'software', [u'Bruker software', u'software']),
    SoftwareName(u'ClinProTools', u'MS:1000711', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'CLINPROT robot', u'MS:1000710', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'GENOLINK', u'MS:1000723', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'GenoTools', u'MS:1000724', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'PolyTools', u'MS:1000727', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'PROTEINEER', u'MS:1000729', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'PureDisk', u'MS:1000735', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'PROTEINEER-LC', u'MS:1000733', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'PROTEINEER spII', u'MS:1000732', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'PROTEINEER fc', u'MS:1000731', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'PROTEINEER dp', u'MS:1000730', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'WARP-LC', u'MS:1000739', 'software',
                 [u'Bruker software', u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'TargetAnalysis', u'MS:1000738', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'HyStar', u'MS:1000817', 'software',
                 [u'Bruker software', u'software']),
    SoftwareName(u'Anubis', u'MS:1002410', 'software', [
                 u'SRM software', u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'ATAQS', u'MS:1000925', 'software',
                 [u'SRM software', u'software']),
    SoftwareName(u'Skyline', u'MS:1000922', 'software', [
                 u'SRM software', u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'TIQAM', u'MS:1000923', 'software',
                 [u'SRM software', u'software']),
    SoftwareName(u'MaRiMba', u'MS:1000872', 'software',
                 [u'SRM software', u'software']),
    SoftwareName(u'MRMaid', u'MS:1002220', 'software',
                 [u'SRM software', u'software']),
    SoftwareName(u'SSRCalc', u'MS:1000874', 'software', [
                 u'peptide attribute calculation software', u'software']),
    SoftwareName(u'SILACAnalyzer', u'MS:1001831', 'software', [
                 u'quantitation software name', u'TOPP software', u'software', u'quantification information', u'analysis software', u'data processing software']),
    SoftwareName(u'Progenesis LC-MS', u'MS:1001830', 'software',
                 [u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'FindPairs', u'MS:1002063', 'software', [
                 u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'Microsoft Excel', u'MS:1002059', 'software', [
                 u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'ProteoSuite', u'MS:1002124', 'software', [
                 u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'x-Tracker', u'MS:1002123', 'software',
                 [u'quantitation software name', u'software', u'quantification information']),
    SoftwareName(u'ITRAQAnalyzer', u'MS:1002129', 'software', [
                 u'quantitation software name', u'TOPP software', u'software', u'quantification information', u'analysis software', u'data processing software']),
    SoftwareName(u'TOPP noise filter', u'MS:1002131', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP spectra filter', u'MS:1002137', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP peak picker', u'MS:1002134', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP map aligner', u'MS:1002147', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MassTraceExtractor', u'MS:1002159', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MzTabExporter', u'MS:1002158', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP IDMerger', u'MS:1002155', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP DTAExtractor', u'MS:1002154', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraMerger', u'MS:1002157', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP IDFileConverter', u'MS:1002156', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP PrecursorMassCorrector', u'MS:1002160', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP HighResPrecursorMassCorrector', u'MS:1002161', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP AdditiveSeries', u'MS:1002162', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP Decharger', u'MS:1002163', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP EICExtractor', u'MS:1002164', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP feature finder', u'MS:1002165', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP ConsensusID', u'MS:1002188', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP IDConflictResolver', u'MS:1002189', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP software adaptor', u'MS:1002180', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpecLibSearcher', u'MS:1002187', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP feature linker', u'MS:1002174', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MapRTTransformer', u'MS:1002173', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP ConsensusMapNormalizer', u'MS:1002172', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP ProteinQuantifier', u'MS:1002171', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP CompNovoCID', u'MS:1002179', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP CompNovo', u'MS:1002178', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP ProteinInference', u'MS:1002203', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FalseDiscoveryRate', u'MS:1002204', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP IDMapper', u'MS:1002191', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP IDFilter', u'MS:1002190', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP IDRTCalibration', u'MS:1002193', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP IDPosteriorErrorProbability', u'MS:1002192', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP PrecursorIonSelector', u'MS:1002195', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP PeptideIndexer', u'MS:1002194', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP OpenSwath component', u'MS:1002197', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MRMMapper', u'MS:1002196', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'OpenXQuest', u'MS:1002673', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'BaselineFilter', u'MS:1000753', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'DBImporter', u'MS:1000755', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'DBExporter', u'MS:1000754', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'FileFilter', u'MS:1000757', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'FileConverter', u'MS:1000756', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'InternalCalibration', u'MS:1000759', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'FileMerger', u'MS:1000758', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'Resampler', u'MS:1000764', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'SpectraFilter', u'MS:1000765', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOFCalibration', u'MS:1000766', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MapAligner', u'MS:1000760', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'MapNormalizer', u'MS:1000761', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'NoiseFilter', u'MS:1000762', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'PeakPicker', u'MS:1000763', 'software', [
                 u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteoWizard SeeMS', u'MS:1002209', 'software', [
                 u'ProteoWizard software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteoWizard msaccess', u'MS:1002208', 'software', [
                 u'ProteoWizard software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteoWizard msconvert', u'MS:1002205', 'software', [
                 u'ProteoWizard software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteoWizard chainsaw', u'MS:1002207', 'software', [
                 u'ProteoWizard software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'ProteoWizard idconvert', u'MS:1002206', 'software', [
                 u'ProteoWizard software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'mzidLib:Omssa2Mzid', u'MS:1002238', 'software',
                 [u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:Tandem2Mzid', u'MS:1002239', 'software', [
                 u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:Csv2Mzid', u'MS:1002240', 'software', [
                 u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:Mzidentml2Csv', u'MS:1002245', 'software', [
                 u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:FalseDiscoveryRate', u'MS:1002244', 'software', [
                 u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:InsertMetaDataFromFasta', u'MS:1002247',
                 'software', [u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:CombineSearchEngines', u'MS:1002246', 'software', [
                 u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:ProteoGrouper', u'MS:1002241', 'software', [
                 u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:Perform emPAI on mzid', u'MS:1002243',
                 'software', [u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'mzidLib:Thresholder', u'MS:1002242', 'software', [
                 u'mzidLib', u'analysis software', u'software']),
    SoftwareName(u'ASAPRatio', u'MS:1002574', 'software', [
                 u'Trans-Proteomic Pipeline software', u'analysis software', u'software']),
    SoftwareName(u'PTMProphet', u'MS:1002292', 'software', [
                 u'Trans-Proteomic Pipeline software', u'analysis software', u'software']),
    SoftwareName(u'XPRESS', u'MS:1002290', 'software', [
                 u'Trans-Proteomic Pipeline software', u'analysis software', u'software']),
    SoftwareName(u'Libra', u'MS:1002291', 'software', [
                 u'Trans-Proteomic Pipeline software', u'analysis software', u'software']),
    SoftwareName(u'PeptideProphet', u'MS:1002287', 'software', [
                 u'Trans-Proteomic Pipeline software', u'analysis software', u'software']),
    SoftwareName(u'ProteinProphet', u'MS:1002289', 'software', [
                 u'Trans-Proteomic Pipeline software', u'analysis software', u'software']),
    SoftwareName(u'iProphet', u'MS:1002288', 'software', [
                 u'Trans-Proteomic Pipeline software', u'analysis software', u'software']),
    SoftwareName(u'TOPP NoiseFilterSGolay', u'MS:1002133', 'software', [
                 u'TOPP noise filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP NoiseFilterGaussian', u'MS:1002132', 'software', [
                 u'TOPP noise filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterMarkerMower', u'MS:1002139', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterBernNorm', u'MS:1002138', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterWindowMower', u'MS:1002146', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterSqrtMower', u'MS:1002144', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterThresholdMower', u'MS:1002145', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterParentPeakMower', u'MS:1002142', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterScaler', u'MS:1002143', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterNLargest', u'MS:1002140', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP SpectraFilterNormalizer', u'MS:1002141', 'software', [
                 u'TOPP spectra filter', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP PeakPickerWavelet', u'MS:1002136', 'software', [
                 u'TOPP peak picker', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP PeakPickerHiRes', u'MS:1002135', 'software', [
                 u'TOPP peak picker', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MapAlignerIdentification', u'MS:1002148', 'software', [
                 u'TOPP map aligner', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MapAlignerPoseClustering', u'MS:1002149', 'software', [
                 u'TOPP map aligner', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MapAlignerSpectrum', u'MS:1002150', 'software', [
                 u'TOPP map aligner', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureFinderCentroided', u'MS:1002166', 'software', [
                 u'TOPP feature finder', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureFinderRaw', u'MS:1002167', 'software', [
                 u'TOPP feature finder', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureFinderIsotopeWavelet', u'MS:1002168', 'software', [
                 u'TOPP feature finder', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureFinderMetabo', u'MS:1002169', 'software', [
                 u'TOPP feature finder', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureFinderMRM', u'MS:1002170', 'software', [
                 u'TOPP feature finder', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MascotAdapter', u'MS:1002182', 'software', [
                 u'TOPP software adaptor', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP MascotAdapterOnline', u'MS:1002183', 'software', [
                 u'TOPP software adaptor', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP InspectAdapter', u'MS:1002181', 'software', [
                 u'TOPP software adaptor', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP XTandemAdapter', u'MS:1002186', 'software', [
                 u'TOPP software adaptor', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP OMSSAAdapter', u'MS:1002184', 'software', [
                 u'TOPP software adaptor', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP PepNovoAdapter', u'MS:1002185', 'software', [
                 u'TOPP software adaptor', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureLinkerUnlabeledQT', u'MS:1002177', 'software', [
                 u'TOPP feature linker', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureLinkerUnlabeled', u'MS:1002176', 'software', [
                 u'TOPP feature linker', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP FeatureLinkerLabeled', u'MS:1002175', 'software', [
                 u'TOPP feature linker', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP OpenSwathFeatureXMLToTSV', u'MS:1002201', 'software', [
                 u'TOPP OpenSwath component', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP OpenSwathDecoyGenerator', u'MS:1002200', 'software', [
                 u'TOPP OpenSwath component', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP OpenSwathRTNormalizer', u'MS:1002202', 'software', [
                 u'TOPP OpenSwath component', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP OpenSwathChromatogramExtractor', u'MS:1002199', 'software', [
                 u'TOPP OpenSwath component', u'TOPP software', u'analysis software', u'data processing software', u'software']),
    SoftwareName(u'TOPP OpenSwathAnalyzer', u'MS:1002198', 'software', [
                 u'TOPP OpenSwath component', u'TOPP software', u'analysis software', u'data processing software', u'software']),
]


software_names_by_name = {c.name: c for c in software_names}


def software_name(name):
    try:
        return software_names_by_name[name]
    except KeyError:
        return SoftwareName(name, name, name, [name])


if __name__ == '__main__':
    __generate_list_code()
