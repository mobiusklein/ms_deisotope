''' Waters
    MassLynx Python SDK
'''
import copy
import ctypes
from ctypes import*
from array import *
from .MassLynxRawReader import *
from .MassLynxRawDefs import MassLynxBaseType


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
            pM = cast(pMasses,POINTER(c_float))
            pI = cast(pIntensities,POINTER(c_float))

            masses = pM[0:size.value]
            intensities = pI[0:size.value]

            # dealocate memory - causes segfaults?
            # MassLynxRawReader.ReleaseMemory( pMasses)
            # MassLynxRawReader.ReleaseMemory( pIntensities)

            return masses, intensities

        except Exception:
            return [], []

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
