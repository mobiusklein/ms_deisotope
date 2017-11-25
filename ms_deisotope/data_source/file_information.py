import os
import hashlib
import warnings

from six import string_types as basestring

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = OSError


id_formats = {
    u'native spectrum identifier format': "Describes how the native spectrum identifiers are formated.",
    u'Agilent MassHunter nativeID format': "Native format defined by scanId=xsd:nonNegativeInteger.",
    u'spectrum from database string nativeID format': "Native format defined by databasekey=xsd:string.",
    u'spectrum from ProteinScape database nativeID format': "Native format defined by databasekey=xsd:long.",
    u'Bruker U2 nativeID format': (
        "Native format defined by declaration=xsd:nonNegativeInteger"
        " collection=xsd:nonNegativeInteger scan=xsd:nonNegativeInteger."),
    u'no nativeID format': (
        "No nativeID format indicates that the file tagged with this term "
        "does not contain spectra that can have a nativeID format."),
    u'Mascot query number': "Native format defined by query=xsd:nonNegativeInteger.",
    u'spectrum from database integer nativeID format': "Native format defined by databasekey=xsd:long.",
    u'Shimadzu Biotech nativeID format': (
        "Native format defined by source=xsd:string "
        "start=xsd:nonNegativeInteger end=xsd:nonNegativeInteger."),
    u'UIMF nativeID format': (
        "Native format defined by frame=xsd:nonNegativeInteger "
        "scan=xsd:nonNegativeInteger frameType=xsd:nonNegativeInteger."),
    u'Bruker Container nativeID format': "Native identifier (UUID).",
    u'SCIEX TOF/TOF nativeID format': (
        "Native format defined by jobRun=xsd:nonNegativeInteger "
        "spotLabel=xsd:string spectrum=xsd:nonNegativeInteger."),
    u'Thermo nativeID format': (
        "Native format defined by controllerType=xsd:nonNegativeInteger "
        "controllerNumber=xsd:positiveInteger scan=xsd:positiveInteger."),
    u'Waters nativeID format': (
        "Native format defined by function=xsd:positiveInteger "
        "process=xsd:nonNegativeInteger scan=xsd:nonNegativeInteger."),
    u'Scaffold nativeID format': "Scaffold native ID format.",
    u'Bruker FID nativeID format': "Native format defined by file=xsd:IDREF.",
    u'Bruker BAF nativeID format': "Native format defined by scan=xsd:nonNegativeInteger.",
    u'Bruker/Agilent YEP nativeID format': "Native format defined by scan=xsd:nonNegativeInteger.",
    u'WIFF nativeID format': (
        "Native format defined by sample=xsd:nonNegativeInteger period=xsd:nonNegativeInteger"
        " cycle=xsd:nonNegativeInteger experiment=xsd:nonNegativeInteger."),
    u'spectrum identifier nativeID format': "Native format defined by spectrum=xsd:nonNegativeInteger.",
    u'scan number only nativeID format': "Native format defined by scan=xsd:nonNegativeInteger.",
    u'single peak list nativeID format': "Native format defined by file=xsd:IDREF.",
    u'multiple peak list nativeID format': "Native format defined by index=xsd:nonNegativeInteger.",
    u'SCIEX TOF/TOF T2D nativeID format': "Native format defined by file=xsd:IDREF.",
}


file_formats = {
    u'Mascot MGF format',
    u'Waters raw format',
    u'DTA format',
    u'SCIEX TOF/TOF T2D format',
    u'Bioworks SRF format',
    u'Andi-MS format',
    u'SCIEX TOF/TOF database',
    u'parameter file',
    u'Bruker/Agilent YEP format',
    u'mz5 format',
    u'ISB mzXML format',
    u'ABI WIFF format',
    u'Bruker BAF format',
    u'UIMF format',
    u'mzML format',
    u'Micromass PKL format',
    u'mass spectrometer file format',
    u'Bruker XML format',
    u'Bruker U2 format',
    u'Bruker FID format',
    u'text format',
    u'PSI mzData format',
    u'Shimadzu Biotech database entity',
    u'Thermo RAW format',
    u'Agilent MassHunter format',
    u'ProteinLynx Global Server mass spectrum XML format',
    u'Bruker Container format',
    u'PerSeptive PKS format',
    u'Phenyx XML format',
    u'Proteinscape spectra',
    u'SCIEX API III format',
    u'SCiLS Lab format',
}


content_keys = [
    "data file content",
    "electromagnetic radiation spectrum",
    "total ion current chromatogram",
    "absorption spectrum",
    "selected reaction monitoring chromatogram",
    "selected ion monitoring chromatogram",
    "consecutive reaction monitoring chromatogram",
    "MSn spectrum",
    "CRM spectrum",
    "SIM spectrum",
    "SRM spectrum",
    "enhanced multiply charged spectrum",
    "constant neutral gain spectrum",
    "constant neutral loss spectrum",
    "charge inversion mass spectrum",
    "time-delayed fragmentation spectrum",
    "emission spectrum",
    "product ion spectrum",
    "precursor ion spectrum",
    "mass spectrum",
    "basepeak chromatogram",
    "selected ion current chromatogram",
    "PDA spectrum",
    "MS1 spectrum",
    "absorption chromatogram",
    "emission chromatogram",
    "mass chromatogram",
    "electromagnetic radiation chromatogram",
    "time-delayed fragmentation spectrum",
    "MSn spectrum",
    "CRM spectrum",
    "SIM spectrum",
    "SRM spectrum",
    "e/2 mass spectrum",
    "constant neutral gain spectrum",
    "constant neutral loss spectrum",
    "charge inversion mass spectrum",
    "product ion spectrum",
    "precursor ion spectrum",
    "MS1 spectrum",
    "enhanced multiply charged spectrum",
    "total ion current chromatogram",
    "selected reaction monitoring chromatogram",
    "selected ion monitoring chromatogram",
    "consecutive reaction monitoring chromatogram",
    "basepeak chromatogram",
    "selected ion current chromatogram",
    "absorption chromatogram",
    "emission chromatogram",
    "time-delayed fragmentation spectrum",
    "enhanced multiply charged spectrum"
]


class FileInformation(object):

    def __init__(self, contents=None, source_files=None):
        if contents is None:
            contents = {}
        if source_files is None:
            source_files = []
        self.contents = dict(contents)
        self.source_files = list(source_files)

    def add_file(self, source, check=True):
        if isinstance(source, basestring):
            source = os.path.realpath(source)
            if check:
                if not os.path.exists(source):
                    raise ValueError(
                        "Source File %r does not exist" % (source,))
            name = os.path.basename(source)
            location = os.path.dirname(source)
            idfmt, file_fmt = SourceFile.guess_format(source)
            source = SourceFile(name, location, name, idfmt, file_fmt)
        self.source_files.append(source)

    def add_content(self, key, value=None):
        self.contents[key] = value

    def remove_content(self, key):
        self.contents.pop(key, None)

    def get_content(self, key):
        return self.contents[key]

    def has_content(self, key):
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


format_parameter_map = {
    "thermo raw": ("Thermo nativeID format", "Thermo RAW format"),
    "agilent d": ("Agilent MassHunter nativeID format", "Agilent MassHunter format"),
}


class SourceFile(object):

    def __init__(self, name, location, id=None, id_format=None, file_format=None, parameters=None):
        if id is None:
            id = name
        if parameters is None:
            parameters = {}
        self.name = name
        self.location = location.replace("file:///", '')
        self.id = id
        self.parameters = dict(parameters)
        self.id_format = id_format or self.infer_id_format_from_paramters()
        self.file_format = file_format or self.infer_file_format_from_parameters()

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

    @staticmethod
    def guess_format(path):
        if not os.path.exists(path):
            return None, None
        if os.path.isdir(path):
            if os.path.exists(os.path.join(path, 'AcqData')):
                return format_parameter_map['agilent d']
        with open(path, 'rb') as fh:
            lead_bytes = fh.read(100)
            decoded = lead_bytes.decode("utf-16")[1:9]
            if decoded == "Finnigan":
                return format_parameter_map['thermo raw']
        return None, None

    def __repr__(self):
        template = "SourceFile(%r, %r, %r, %r, %r%s)"
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
        hasher = hashlib.new(hash_type)
        buffer_size = int(buffer_size)
        with open(self.path, 'rb') as fh:
            content_buffer = fh.read(buffer_size)
            while content_buffer:
                hasher.update(content_buffer)
                content_buffer = fh.read(buffer_size)
        return hasher.hexdigest()

    def add_checksum(self, hash_type='sha1'):
        checksum = self._compute_checksum(hash_type)
        if hash_type == 'sha1':
            self.parameters['SHA-1'] = checksum
        elif hash_type == "md5":
            self.parameters['MD5'] = checksum

    def validate_checksum(self):
        if not os.path.exists(self.path):
            FileNotFoundError("%s not found" % (self.path,))
        if 'SHA-1' in self.parameters:
            checksum = self._compute_checksum('sha1')
            return self.parameters['SHA-1'] == checksum
        elif 'MD5' in self.parameters:
            checksum = self._compute_checksum("md5")
            return self.parameters['MD5'] == checksum
        else:
            warnings.warn("%r did not have a reference checksum. Could not validate" % (self,))
            return True
