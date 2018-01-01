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
    'CRM spectrum',
    'emission spectrum',
    'emission chromatogram',
    'data file content',
    'constant neutral loss spectrum',
    'enhanced multiply charged spectrum',
    'e/2 mass spectrum',
    'product ion spectrum',
    'time-delayed fragmentation spectrum',
    'PDA spectrum',
    'consecutive reaction monitoring chromatogram',
    'MS1 spectrum',
    'SRM spectrum',
    'selected ion current chromatogram',
    'absorption chromatogram',
    'electromagnetic radiation chromatogram',
    'charge inversion mass spectrum',
    'selected reaction monitoring chromatogram',
    'basepeak chromatogram',
    'MSn spectrum',
    'mass spectrum',
    'constant neutral gain spectrum',
    'absorption spectrum',
    'SIM spectrum',
    'selected ion monitoring chromatogram',
    'total ion current chromatogram',
    'mass chromatogram',
    'electromagnetic radiation spectrum',
    'precursor ion spectrum'
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
            source = SourceFile.from_path(source)
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

    def copy(self):
        return self.__class__(
            self.contents.copy(), [f.copy() for f in self.source_filesr])


format_parameter_map = {
    "thermo raw": ("Thermo nativeID format", "Thermo RAW format"),
    "agilent d": ("Agilent MassHunter nativeID format", "Agilent MassHunter format"),
}


class SourceFile(object):

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
            ext = parts[1]
            if ext.lower() == '.mzml':
                fmt = "mzML format"
                id_fmt = "no nativeID format"
                hit = False
                with open(path, 'rb') as fh:
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
            warnings.warn("%r did not have a reference checksum. Could not validate" % (self,))
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
                              self.parameters.copy(), self.id_format,
                              self.file_format)
