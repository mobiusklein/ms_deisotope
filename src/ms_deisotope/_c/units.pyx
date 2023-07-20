# cython: language_level=3
cimport cython

from cpython.object cimport PyObject
from cpython.dict cimport PyDict_GetItem, PyDict_SetItem
from cpython.list cimport PyList_Append


from ms_deisotope._c.compat cimport PyStr_AsUTF8AndSize, PyStr_FromString


cpdef reconstruct_unit_struct(cls, value, unit_info):
    inst = cls(value)
    inst.unit_info = unit_info
    return inst


cpdef clear_unit_cv_table():
    _UNIT_CV_INTERN_TABLE.clear()


cdef str _intern_unit_or_cv(str unit_or_cv):
    """Intern `unit_or_cv` in :const:`~._UNIT_CV_INTERN_TABLE`, potentially
    keeping a reference to the object stored for the duration of the program.
    Parameters
    ----------
    unit_or_cv : object
        The value to intern
    Returns
    -------
    object:
        The object which `unit_or_cv` hash-equals in :const:`~._UNIT_CV_INTERN_TABLE`.
    """
    cdef:
        PyObject* ptmp
    if unit_or_cv is None:
        return None
    ptmp = PyDict_GetItem(_UNIT_CV_INTERN_TABLE, unit_or_cv)
    if ptmp == NULL:
        PyDict_SetItem(_UNIT_CV_INTERN_TABLE, unit_or_cv, unit_or_cv)
        return <str>PyDict_GetItem(_UNIT_CV_INTERN_TABLE, unit_or_cv)
    else:
        return <str>ptmp


class UnitInt(int):
    """Represents an integer value with a unit name.
    Behaves identically to a built-in :class:`int` type.
    Attributes
    ----------
    unit_info : :class:`str`
        The name of the unit this value posseses.
    """

    @classmethod
    def create(cls, value, unit_info):
        cdef object self = UnitInt(value)
        self.unit_info = _intern_unit_or_cv(unit_info)
        return self

    def __reduce__(self):
        return reconstruct_unit_struct, (self.__class__, int(self), self.unit_info)

    def _repr_pretty_(self, p, cycle):
        base = super(UnitInt, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


@cython.final
@cython.no_gc
cdef class UnitFloat(float):
    """Represents an float value with a unit name.
    Behaves identically to a built-in :class:`float` type.
    Attributes
    ----------
    unit_info : :class:`str`
        The name of the unit this value posseses.
    """
    cdef:
        public str unit_info

    @classmethod
    def create(cls, value, unit_info):
        cdef UnitFloat self = UnitFloat(value)
        self.unit_info = _intern_unit_or_cv(unit_info)
        return self

    def __reduce__(self):
        return reconstruct_unit_struct, (self.__class__, float(self), self.unit_info)

    def _repr_pretty_(self, p, cycle):
        base = super(UnitFloat, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


@cython.final
@cython.no_gc
cdef class UnitStr(str):
    """Represents an string value with a unit name.
    Behaves identically to a built-in :class:`str` type.
    Attributes
    ----------
    unit_info : :class:`str`
        The name of the unit this value posseses.
    """
    cdef:
        public str unit_info

    @classmethod
    def create(cls, value, unit_info):
        cdef UnitStr self = UnitStr(value)
        self.unit_info = _intern_unit_or_cv(unit_info)
        return self

    def __reduce__(self):
        return reconstruct_unit_struct, (self.__class__, str(self), self.unit_info)

    def _repr_pretty_(self, p, cycle):
        base = super(UnitStr, self).__repr__()
        if self.unit_info:
            string = "%s %s" % (base, self.unit_info)
        else:
            string = base
        p.text(string)


cpdef reconstruct_cvstr(value, accession, unit_accession):
    return CVStr.create(value, accession, unit_accession)


cdef dict _UNIT_CV_INTERN_TABLE = dict()


cdef dict CVStr_cache = {}


@cython.final
@cython.no_gc
cdef class CVStr(str):
    """A helper class to associate a controlled vocabullary accession
    number with an otherwise plain :class:`str` object
    Attributes
    ----------
    accession : str
        The accession number for this parameter, e.g. MS:1000040
    unit_accession : str
        The accession number for the unit of the value, if any
    """
    cdef:
        public str accession
        public str unit_accession

    @classmethod
    def create(cls, value, accession=None, unit_accession=None):
        cdef:
            PyObject* ptmp
            CVStr self
        accession = _intern_unit_or_cv(accession)
        unit_accession = _intern_unit_or_cv(unit_accession)
        ptmp = PyDict_GetItem(CVStr_cache, value)
        if ptmp != NULL:
            self = <CVStr>ptmp
            if self.accession is accession and self.unit_accession is unit_accession:
                return  self

        self = CVStr(value)
        self.accession = accession
        self.unit_accession = unit_accession

        PyDict_SetItem(CVStr_cache, self, self)
        return self

    def __reduce__(self):
        return reconstruct_cvstr, (str(self), self.accession, self.unit_accession)

    def _repr_pretty_(self, p, cycle):
        base = super(CVStr, self).__repr__()
        p.text(base)


cdef class UnitPrimitiveWrapperMeta(type):
    cdef:
        public type wrapped

    def __init__(self, name, bases, state):
        super().__init__(name, bases, state)
        self.wrapped = state['wrapped']

    def __instancecheck__(cls, inst):
        return isinstance(inst, cls.wrapped)

    def __call__(cls, value, unit_info=None):
        return cls.wrapped.create(value, unit_info)

    def __eq__(self, other):
        if isinstance(other, UnitPrimitiveWrapperMeta):
            return self.wrapped == other.wrapped
        else:
            return self.wrapped == other


class unitint(metaclass=UnitPrimitiveWrapperMeta):
    wrapped = UnitInt


class unitfloat(metaclass=UnitPrimitiveWrapperMeta):
    wrapped = UnitFloat


class unitstr(metaclass=UnitPrimitiveWrapperMeta):
    wrapped = UnitStr


cdef class CVPrimitiveWrapperMeta(type):
    cdef:
        public type wrapped

    def __init__(self, name, bases, state):
        super().__init__(name, bases, state)
        self.wrapped = state['wrapped']

    def __instancecheck__(cls, inst):
        return isinstance(inst, cls.wrapped)

    def __call__(cls, value, accession=None, unit_accession=None):
        return cls.wrapped.create(value, accession, unit_accession)

    def __eq__(self, other):
        if isinstance(other, UnitPrimitiveWrapperMeta):
            return self.wrapped == other.wrapped
        else:
            return self.wrapped == other

    def __hash__(self):
        return hash(self.__qualname__)

class cvstr(metaclass=CVPrimitiveWrapperMeta):
    wrapped = CVStr


@cython.freelist(10)
cdef class _XMLParam:
    cdef:
        public CVStr name
        public object value
        public str type

    cpdef bint is_empty(self):
        value = self.value
        return value == "" or value is None

    def __iter__(self):
        yield self.name
        yield self.value
        yield self.type

    def __init__(self, name, value, type):
        self.name = name
        self.value = value
        self.type = type

    @staticmethod
    cdef _XMLParam _create(CVStr name, object value, str type):
        cdef _XMLParam self = _XMLParam.__new__(_XMLParam)
        self.name = name
        self.value = value
        self.type = type
        return self

    def __repr__(self):
        return f"{self.__class__.__name__}({self.name}, {self.value}, {self.type})"


cpdef str _local_name(element):
    cdef:
        Py_ssize_t i, n
        char* cname
        str tag, result

    tag = element.tag
    if tag and tag[0] == '{':
        cname = PyStr_AsUTF8AndSize(tag, &n)
        for i in range(n - 1, 0, -1):
            if cname[i] == 125: # '}'
                result = PyStr_FromString(&cname[i + 1])
                return result
        return tag
    else:
        return tag


@cython.binding(True)
def _handle_param(self, element, **kwargs):
    """Unpacks cvParam and userParam tags into key-value pairs"""
    cdef:
        object value, attribs
        str unit_info, unit_accession, accession, param_type

    attribs = element.attrib
    unit_info = None
    unit_accesssion = None
    if 'unitCvRef' in attribs or 'unitName' in attribs:
        unit_accesssion = attribs.get('unitAccession')
        unit_name = attribs.get('unitName', unit_accesssion)
        unit_info = unit_name

    accession = attribs.get('accession')
    value = attribs.get('value', '')
    param_type = attribs.get('type')

    if param_type is not None:
        if param_type == 'xsd:int':
            value = UnitInt.create(value, unit_info)
        elif param_type == 'xsd:float':
            value = UnitFloat.create(value, unit_info)
        elif param_type == 'xsd:string':
            value = UnitStr.create(value, unit_info)
        else:
            try:
                if value and value[0].isdigit():
                    value = UnitFloat.create(value, unit_info)
                else:
                    value = UnitStr.create(value, unit_info)
            except ValueError:
                value = UnitStr.create(value, unit_info)
    else:
        try:
            if value and value[0].isdigit():
                value = UnitFloat.create(value, unit_info)
            else:
                value = UnitStr.create(value, unit_info)
        except ValueError:
            value = UnitStr.create(value, unit_info)

    return _XMLParam._create(CVStr.create(attribs['name'], accession, unit_accesssion), value, _local_name(element))


@cython.binding(True)
def _handle_param_fill_missing_value(self, element, **kwargs):
    attribs = element.attrib
    if "value" not in attribs:
        attribs['value'] = ''
    return _handle_param(self, element)


def patch_pyteomics():
    # Cannot patch here or else we can introduce issues
    # with pickling.
    # from pyteomics.auxiliary import structures
    # structures.unitint = unitint
    # structures.unitfloat = unitfloat
    # structures.unitstr = unitstr
    # structures.cvstr = cvstr
    from pyteomics import xml
    xml.unitint = unitint
    xml.unitfloat = unitfloat
    xml.unitstr = unitstr
    xml.cvstr = cvstr
    xml.XML._handle_param = _handle_param
    xml._local_name = _local_name