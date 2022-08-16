#import copy
import ctypes
from ctypes import*
#from array import *
from .MassLynxRawReader import *

class MassLynxParameters(object):
    """Create parameter"""

    def __init__(self):
        self._mlParameters = c_void_p()              # instance variable
        self._codeHandler = MassLynxCodeHandler()
        self._stringHandler = MassLynxStringHandler()
        createParameters = MassLynxRawReader.massLynxDll.createParameters
        createParameters.argtypes = [POINTER(c_void_p)]
        self.CheckReturnCode(createParameters(self._mlParameters))

    # destroy the samplelist
    def __del__(self):
        destroyParameters = MassLynxRawReader.massLynxDll.destroyParameters
        destroyParameters.argtypes = [c_void_p]
        destroyParameters( self._mlParameters )

    # get the parameter object
    def GetParameters(self):
        return self._mlParameters

    # handle the retun code - no exceptions
    def CheckReturnCode( self, code ):
        return self._codeHandler.CheckReturnCode(code, False)

    # Set the key value pair
    def Set( self, key, value ):
        # is it a bool
        if isinstance(value, (bool)):
            if ( value ):
                value = 1
            else:
                value = 0

        # convert to string and encode
        bytes = str.encode(str(value))

        setParameterValue = MassLynxRawReader.massLynxDll.setParameterValue
        setParameterValue.argtypes = [c_void_p, c_int, c_char_p]
        setParameterValue( self._mlParameters, key, bytes )

        return self

    def Get( self, key ):
        # call the dll
        value = c_char_p()
        getParameterValue = MassLynxRawReader.massLynxDll.getParameterValue
        getParameterValue.argtypes = [c_void_p, c_int, POINTER(c_char_p)]
        getParameterValue(self._mlParameters, key, value )

        return self._stringHandler.ToString(value, False)