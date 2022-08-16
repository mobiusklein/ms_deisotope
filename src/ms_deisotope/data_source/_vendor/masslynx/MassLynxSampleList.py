import copy
import ctypes
from ctypes import*
from array import *
from .MassLynxRawReader import *

class MassLynxSampleList(object):
    """Class to create MassLynx Sample lists"""

    def __init__(self):
        self._mlSampleList = c_void_p()              # instance variable
        self._codeHandler = MassLynxCodeHandler()
        self._stringHandler = MassLynxStringHandler()
        createSampleList = MassLynxRawReader.massLynxDll.createSampleList
        createSampleList.argtypes = [POINTER(c_void_p)]
        self.CheckReturnCode(createSampleList(self._mlSampleList))

    # destroy the samplelist
    def __del__(self):
        destroySampleList = MassLynxRawReader.massLynxDll.destroySampleList
        destroySampleList.argtypes = [c_void_p]
        destroySampleList( self._mlSampleList )

    # handle the retun code - no exceptions
    def CheckReturnCode( self, code ):
        return self._codeHandler.CheckReturnCode(code, False)

    # get the error message
    def GetLastMessage(self):
        return self._codeHandler.GetLastMessage()

    # get the error code
    def GetLastCode(self):
        return self._codeHandler.GetLastCode()

    # get the sample list in csv
    def ToString( self ):
        temp = c_char_p()
        sampleListToString = MassLynxRawReader.massLynxDll.sampleListToString
        sampleListToString.argtypes = [c_void_p, POINTER(c_char_p)]
        self.CheckReturnCode( sampleListToString(self._mlSampleList, temp) )
        return self._stringHandler.ToString( temp, False )

    def AddRow( self, row ):
        addSampleListRow = MassLynxRawReader.massLynxDll.addSampleListRow
        addSampleListRow.argtypes = [c_void_p, c_void_p]
        return self.CheckReturnCode( addSampleListRow(self._mlSampleList, row.GetParameters()) )
