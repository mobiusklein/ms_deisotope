''' Waters
    MassLynx Python SDK
'''

import string
import ctypes

from ctypes import *
from .libload import proxy, register_dll

#from enum import IntEnum

#class MassLynxScanItem(IntEnum):
#    LINEAR_DETECTOR_VOLTAGE = 1
#    LINEAR_SENSITIVITY = 2
#    REFLECTRON_LENS_VOLTAGE = 3
#    REFLECTRON_DETECTOR_VOLTAGE = 4
#    REFLECTRON_SENSITIVITY = 5
#    LASER_REPETITION_RATE = 6
#    COURSE_LASER_CONTROL = 7
#    FINE_LASER_CONTROL = 8
#    LASERAIM_XPOS = 9
#    LASERAIM_YPOS = 10
#    NUM_SHOTS_SUMMED = 11
#    NUM_SHOTS_PERFORMED = 12
#    SEGMENT_NUMBER = 13
#    LCMP_TFM_WELL = 14
#    SEGMENT_TYPE = 15
#    SOURCE_REGION1 = 16
#    SOURCE_REGION2 = 17
#    REFLECTRON_FIELD_LENGTH = 18
#    REFLECTRON_LENGTH = 19
#    REFLECTRON_VOLT = 20
#    SAMPLE_PLATE_VOLT = 21
#    REFLECTRON_FIELD_LENGTH_ALT = 22
#    REFLECTRON_LENGTH_ALT = 23
#    PSD_STEP_MAJOR = 24
#    PSD_STEP_MINOR = 25
#    PSD_FACTOR_1 = 26
#    NEEDLE = 50
#    COUNTER_ELECTRODE_VOLTAGE = 51
#    SAMPLING_CONE_VOLTAGE = 52
#    SKIMMER_LENS = 53
#    SKIMMER = 54
#    PROBE_TEMPERATURE = 55
#    SOURCE_TEMPERATURE = 56
#    RF_VOLTAGE = 57
#    SOURCE_APERTURE = 58
#    SOURCE_CODE = 59
#    LM_RESOLUTION = 60
#    HM_RESOLUTION = 61
#    COLLISION_ENERGY = 62
#    ION_ENERGY = 63
#    MULTIPLIER1 = 64
#    MULTIPLIER2 = 65
#    TRANSPORTDC = 66
#    TOF_APERTURE = 67
#    ACC_VOLTAGE = 68
#    STEERING = 69
#    FOCUS = 70
#    ENTRANCE = 71
#    GUARD = 72
#    TOF = 73
#    REFLECTRON = 74
#    COLLISION_RF = 75
#    TRANSPORT_RF = 76
#    SET_MASS = 77
#    COLLISION_ENERGY2 = 78
#    SONAR_ENABLED = 81
#    QUAD_START_MASS = 82
#    QUAD_STOP_MASS = 83
#    QUAD_PEAK_WIDTH = 84
#    REFERENCE_SCAN = 100
#    USE_LOCKMASS_CORRECTION = 101
#    LOCKMASS_CORRECTION = 102
#    USETEMP_CORRECTION = 103
#    TEMP_CORRECTION = 104
#    TEMP_COEFFICIENT = 105
#    FAIMS_COMPENSATION_VOLTAGE = 106
#    TIC_TRACE_A = 107
#    TIC_TRACE_B = 108
#    RAW_EE_CV = 110
#    RAW_EE_CE = 111
#    ACCURATE_MASS = 112
#    ACCURATE_MASS_FLAGS = 113
#    SCAN_ERROR_FLAG = 114
#    DRE_TRANSMISSION = 115
#    SCAN_PUSH_COUNT = 116
#    RAW_STAT_SWAVE_NORMALISATION_FACTOR = 117
#    MIN_DRIFT_TIME_CHANNEL = 122
#    MAX_DRIFT_TIME_CHANNEL = 123
#    TOTAL_ION_CURRENT = 252
#    BASE_PEAK_MASS = 253
#    BASE_PEAK_INTENSITY = 254
#    PEAKS_IN_SCAN = 255
#    UNINITIALISED = 1000


class MassLynxException(Exception):
    code = 0
    def __init__(self, code, message):
        self.message = message
        self.code = code
    def __str__(self):
        return repr(  self.message )

    def Handler(self):
        #print( self.message )
        # can rethrow if needed
        return

# string handler
class MassLynxStringHandler(object):

    def __init__(self):
        return

    def ToString(self, chString, release):
        if (None == chString):
            return ""

        strValue  = chString.value.decode()
        if (release ):
            MassLynxRawReader.ReleaseMemory(chString)

        chString = None

        return strValue


class MassLynxCodeHandler(object):
    def __init__(self):
        self._code = 0
        self._stringHandler = MassLynxStringHandler()
        return

    # three option true, false, throw exception
    def CheckReturnCode( self, code, throw = True ):
        self._code = code;
        if (0 == code):
            return True

        if (throw):
            raise MassLynxException( self.GetLastCode(), self.GetLastMessage() )

        # get last error
        return False

    def GetLastCode(self):
        return self._code

    def GetLastMessage(self):
        # load the dll
        getErrorMessage = MassLynxRawReader.massLynxDll.getErrorMessage
        getErrorMessage.argtypes = [ c_int, POINTER(c_char_p)]

        message = (c_char_p)()
        getErrorMessage( self.GetLastCode(), message )
  #      exceptionMessage = "MassLynx Exception {} : {}".format( returnCode, message.value.decode())

        # release the memory
        return self._stringHandler.ToString(message, True)



class MassLynxRawReader(object):
    """basic functionality to read raw files"""

    # load the dll
    massLynxDll = proxy
    version = '1.0'                                 # class variable

    def __init__(self, source, mlType):

        self.mlRawReader = c_void_p()              # instance variable
        self._codeHandler = MassLynxCodeHandler()
        self._stringHandler = MassLynxStringHandler()

        # create scan reader from a path
        if (isinstance(source, str) ):
            bytes = str.encode(source)
            createRawReaderFromPath = MassLynxRawReader.massLynxDll.createRawReaderFromPath
            createRawReaderFromPath.argtypes = [ c_char_p, POINTER(c_void_p), c_int]
            self._codeHandler.CheckReturnCode(createRawReaderFromPath(bytes, self._getReader(), mlType))

        # create scan reader from a reader
        elif (isinstance(source, MassLynxRawReader)):
            createRawReaderFromReader = MassLynxRawReader.massLynxDll.createRawReaderFromReader
            createRawReaderFromReader.argtypes = [ c_void_p, POINTER(c_void_p), c_int]
            self._codeHandler.CheckReturnCode(createRawReaderFromReader(source._getReader(),self._getReader(),mlType))

        # did we fall through
        else:
            self._codeHandler.CheckReturnCode( 1 )

        return

    # destroy the reader
    def __del__(self):
        # destroy reader
        destroyRawReader = MassLynxRawReader.massLynxDll.destroyRawReader
        destroyRawReader.argtypes = [c_void_p]
        destroyRawReader( self._getReader() )

    def _getReader(self):
        return self.mlRawReader

    def ToString( self, pString, release = False ):
        return self._stringHandler.ToString( pString, release )

    def CheckReturnCode(self, code, throw = True):
        return self._codeHandler.CheckReturnCode( code , throw )

    # common util to free memory
    @staticmethod
    def ReleaseMemory( address):
        releaseMemory = MassLynxRawReader.massLynxDll.releaseMemory
        releaseMemory.argtypes = [ c_void_p]
        releaseMemory( address )