''' Waters
    MassLynx Python SDK
'''
import sys
import copy
import ctypes
from ctypes import*
from array import *
from .MassLynxRawReader import *
from .MassLynxRawDefs import MassLynxBaseType

import numpy as np

if sys.version_info.major >= 3:
    buf_from_mem = ctypes.pythonapi.PyMemoryView_FromMemory
    buf_from_mem.restype = ctypes.py_object
    buf_from_mem.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
else:
    buf_from_mem = ctypes.pythonapi.PyBuffer_FromMemory
    buf_from_mem.restype = ctypes.py_object


def make_nd_array(c_pointer, shape, dtype=np.float64, order='C', own_data=True):
    arr_size = np.prod(shape[:]) * np.dtype(dtype).itemsize
    if sys.version_info.major >= 3:
        buffer = buf_from_mem(c_pointer, arr_size, 0x100)
    else:
        buffer = buf_from_mem(c_pointer, arr_size)
    arr_size = np.prod(shape[:]) * np.dtype(dtype).itemsize
    arr = np.ndarray(tuple(shape[:]), dtype, buffer, order=order)
    if own_data and not arr.flags.owndata:
        return arr.copy()
    else:
        return arr


class MassLynxRawScanReader(MassLynxRawReader):
    """Read masslynx scan data"""

    def __init__(self, source ):
        super(MassLynxRawScanReader, self).__init__(source, MassLynxBaseType.SCAN)

    #@classmethod
    #def CreateFromPath( cls, path ):                      # alternative constructor - pass class to constructor
    #    return cls(MassLynxRawReader.fromPath( path, 1 ))                     # initalise with reader

    #@classmethod
    #def CreateFromReader( cls, sourceReader ):                  # alternative constructor - pass class to constructor
    #     return cls(MassLynxRawReader.fromReader( sourceReader, 1 ))                     # initalise with reader

    def ReadScan( self, whichFunction, whichScan ):
        masses = []
        intensities = []

        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()

        # read scan
        readScan = MassLynxRawReader.massLynxDll.readScan
        readScan.argtypes = [c_void_p, c_int, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        super(MassLynxRawScanReader, self).CheckReturnCode( readScan(self._getReader(),whichFunction,whichScan,pMasses,pIntensities,size) )

        # fill the array
        pM = cast(pMasses,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        masses = pM[0:size.value]
        intensities = pI[0:size.value]

        # dealocate memory
        #MassLynxRawReader.ReleaseMemory( pMasses)
        #MassLynxRawReader.ReleaseMemory( pIntensities)

        # masses = make_nd_array(pMasses, (size, ), np.float32)
        # intensities = make_nd_array(pIntensities, (size, ), np.float32)
        return masses, intensities

    def ReadScanFlags( self, whichFunction, whichScan ):
        masses = []
        intensities = []
        flags = []

        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()
        pFlags = c_void_p()

        # read scan
        readScanFlags = MassLynxRawReader.massLynxDll.readScanFlags
        readScanFlags.argtypes = [c_void_p, c_int, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        super(MassLynxRawScanReader, self).CheckReturnCode( readScanFlags(self._getReader(),whichFunction,whichScan,pMasses,pIntensities,pFlags,size) )

        # fill the array
        pM = cast(pMasses,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        masses = pM[0:size.value]
        intensities = pI[0:size.value]

        # check for flags
        if None != pFlags.value:
            pF = cast(pFlags,POINTER(c_byte))
            flags = pF[0:size.value]


        return masses, intensities, flags

    def ReadDriftScan( self, whichFunction, whichScan, whichDrift ):
        masses = []
        intensities = []

        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()
        # read scan
        readDriftScan = MassLynxRawReader.massLynxDll.readDriftScan
        readDriftScan.argtypes = [c_void_p, c_int, c_int, c_int, POINTER(c_void_p), POINTER(c_void_p), POINTER(c_int)]
        try:
            super(MassLynxRawScanReader, self).CheckReturnCode( readDriftScan(self._getReader(),whichFunction, whichScan, whichDrift, pMasses,pIntensities,size) )# fill the array
            # pM = cast(pMasses,POINTER(c_float))
            # pI = cast(pIntensities,POINTER(c_float))
            # masses = pM[0:size.value]
            # intensities = pI[0:size.value]

            masses = make_nd_array(pMasses, (size.value, ), np.float32)
            intensities = make_nd_array(pIntensities, (size.value, ), np.float32)

            # dealocate memory - causes segfaults?
            # MassLynxRawReader.ReleaseMemory( pMasses)
            # MassLynxRawReader.ReleaseMemory( pIntensities)

            return masses, intensities

        except Exception as err:
            return np.array([], dtype=np.float32), np.array([], dtype=np.float32)

    def ReadDriftScanFlagsIndex(self, whichFunction, whichScan, whichDrift):
        masses = []
        intensities = []
        flags = []

        # create the retrun values
        size = c_int(0)
        pMasses = c_void_p()
        pIntensities = c_void_p()
        pFlags = c_void_p()

        readDriftScanFlagsIndex = MassLynxRawReader.massLynxDll.readDriftScanFlagsIndex
        readDriftScanFlagsIndex.argtypes = [c_void_p, c_int, c_int, c_int, POINTER(c_void_p), POINTER(c_void_p),
                                            POINTER(c_void_p), POINTER(c_int)]

        super(MassLynxRawScanReader, self).CheckReturnCode(
            readDriftScanFlagsIndex(self._getReader(),whichFunction,whichScan,whichDrift,pMasses,pIntensities,pFlags,size))        # fill the array

        pM = cast(pMasses,POINTER(c_float))
        pI = cast(pIntensities,POINTER(c_float))

        masses = pM[0:size.value]
        intensities = pI[0:size.value]

        # check for flags
        if pFlags.value is not None:
            pF = cast(pFlags,POINTER(c_byte))
            flags = pF[0:size.value]
        # Need to manage memory here?
        return masses, intensities, flags


    #def readDaughterScan( self, whichFunction, whichScan ):
    #    try:
    #        size = self.getScanSize( whichFunction,whichScan )

    #       # get the daughter scan size
    #        daughtersize = c_int(0)

    #        # get daughter size
    #        getDaughterScanSize =  RawReader.massLynxDll.getDaughterScanSize
    #        getDaughterScanSize.argtypes = [c_void_p, c_int, c_int, POINTER(c_int)]
    #        RawReader.CheckReturnCode( getDaughterScanSize(RawReader.getReader(self),whichFunction,whichScan, daughtersize) )

       #      # create the float arrays
    #        masses = (c_float*size)()
    #        intensities = (c_float*size)()
    #        daughters = (c_float*daughtersize.value)()

    #        # read daughter size
    #        readSpectrumDaughters = RawReader.massLynxDll.readSpectrumDaughters
    #        readSpectrumDaughters.argtypes = [c_void_p, c_int, c_int,  POINTER(c_float), POINTER(c_float), POINTER(c_float)]
    #        RawReader.CheckReturnCode( readSpectrumDaughters(RawReader.getReader(self), whichFunction, whichScan, masses, intensities, daughters) )

    #    except RawReaderException as e:
    #        e.Handler()
    #        return [], [], []

    #    return list(masses), list(intensities), list(daughters)
