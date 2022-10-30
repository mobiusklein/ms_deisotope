cdef extern from "compat.h":
    char* PyStr_AsString(str string)
    char* PyStr_AsUTF8AndSize(str string, Py_ssize_t*)
    str PyStr_FromString(char* string)
    long PyInt_AsLong(object i)
    object PyInt_FromLong(long i)
    str PyStr_Format(object format, object args)
    str PyStr_FromStringAndSize(char* string, Py_ssize_t size)
    cdef int IS_PY3