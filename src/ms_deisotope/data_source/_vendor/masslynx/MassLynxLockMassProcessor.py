'''
    Waters
    MassLynx Python SDK
'''

import copy
import ctypes
from ctypes import*
from array import *
from .MassLynxRawReader import *
from .MassLynxRawDefs import MassLynxBaseType

class MassLynxLockMassProcessor(object):
    """Class to enable post acquisition lock mass correction"""

    def __init__(self):
        self._codeHandler = MassLynxCodeHandler()
        self.mlLockMassProcessor = c_void_p()
        createRawProcessor = MassLynxRawReader.massLynxDll.createRawProcessor
        createRawProcessor.argtypes = [POINTER(c_void_p), c_int, c_void_p, POINTER(c_void_p)]
        self._codeHandler.CheckReturnCode(createRawProcessor(self.mlLockMassProcessor, MassLynxBaseType.LOCKMASS, c_void_p(0), c_void_p(0) ))

    # destroy the processor
    def __del__(self):
        destroyLockMassProcessor = MassLynxRawReader.massLynxDll.destroyRawProcessor
        destroyLockMassProcessor.argtypes = [c_void_p]
        destroyLockMassProcessor( self.mlLockMassProcessor )

    def SetRawData( self, source ):
        # set data from a path from a path
        if (isinstance(source, str) ):
            self._setRawPath(source)

        # set data from a reader
        elif (isinstance(source, MassLynxRawReader)):
            self._setRawReader(source)

        # did we fall through
        else:
            self._codeHandler.CheckReturnCode( 1 )

        return

    def _setRawPath( self, path ):
        bytes = str.encode(path)
        setRawPath = MassLynxRawReader.massLynxDll.setRawPath
        setRawPath.argtypes = [c_void_p, c_char_p]
        self._codeHandler.CheckReturnCode(setRawPath(self.mlLockMassProcessor, bytes) )


    def _setRawReader( self, mlReader ):
        setRawReader = MassLynxRawReader.massLynxDll.setRawReader
        setRawReader.argtypes = [c_void_p, c_void_p]
        self._codeHandler.CheckReturnCode(setRawReader(self.mlLockMassProcessor, mlReader._getReader() ))

    def SetParameters(self, parameters):
        p_params = parameters.GetParameters()
        setLockMassParameters = MassLynxRawReader.massLynxDll.setLockMassParameters
        setLockMassParameters.argtypes = [c_void_p, c_void_p]
        self._codeHandler.CheckReturnCode(
            setLockMassParameters(self.mlLockMassProcessor, p_params))

    def LockMassCorrect( self ):
         # get applied lock mass gain
        success = c_bool(0)
        lockMassCorrect =  MassLynxRawReader.massLynxDll.lockMassCorrect
        lockMassCorrect.argtypes = [c_void_p, POINTER(c_bool)]
        self._codeHandler.CheckReturnCode(lockMassCorrect(self.mlLockMassProcessor, success ))

        return success.value

    def RemoveLockMassCorrection( self ):
        # get applied lock mass gain
        removeLockMassCorrection =  MassLynxRawReader.massLynxDll.removeLockMassCorrection
        removeLockMassCorrection.argtypes = [c_void_p]
        self._codeHandler.CheckReturnCode( removeLockMassCorrection(self.mlLockMassProcessor) )

    def IsLockMassCorrected( self ):
        # get applied lock mass gain
        applied = c_int(0)
        isLockMassCorrected =  MassLynxRawReader.massLynxDll.LMP_isLockMassCorrected
        isLockMassCorrected.argtypes = [c_void_p, POINTER(c_int)]
        self._codeHandler.CheckReturnCode( isLockMassCorrected(self.mlLockMassProcessor, applied) )

        return applied.value == 1

    def CanLockMassCorrect( self ):
        # get applied lock mass gain
        canApply = c_int(0)
        canLockMassCorrect =  MassLynxRawReader.massLynxDll.LMP_canLockMassCorrect
        canLockMassCorrect.argtypes = [c_void_p, POINTER(c_int)]
        self._codeHandler.CheckReturnCode( canLockMassCorrect(self.mlLockMassProcessor, canApply) )

        return canApply.value == 1

    def GetLockMassValues( self ):
        # get applied lock mass gain
        mass = c_float(0)
        tolerance = c_float(0)
        getLockMassValues =  MassLynxRawReader.massLynxDll.getLockMassValues
        getLockMassValues.argtypes = [c_void_p, POINTER(c_float), POINTER(c_float)]
        self._codeHandler.CheckReturnCode( getLockMassValues(self.mlLockMassProcessor, mass, tolerance) )

        return mass.value, tolerance.value

    def GetLockMassCorrection( self, retentionTime ):
        # get applied lock mass gain
        gain = c_float(0)
        getLockMassCorrection =  MassLynxRawReader.massLynxDll.getLockMassCorrection
        getLockMassCorrection.argtypes = [c_void_p, c_float, POINTER(c_float)]
        self._codeHandler.CheckReturnCode( getLockMassCorrection(self.mlLockMassProcessor, retentionTime, gain) )

        return gain.value
