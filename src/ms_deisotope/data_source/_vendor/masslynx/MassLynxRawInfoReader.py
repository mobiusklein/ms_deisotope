'''
     Waters
    MassLynx Python SDK
'''

import ctypes
from ctypes import*

from .MassLynxRawReader import MassLynxRawReader
from .MassLynxRawDefs import MassLynxBaseType
from .MassLynxParameters import MassLynxParameters

from enum import IntEnum


class MassLynxRawInfoReader(MassLynxRawReader):

    def __init__(self, source ):
        super(MassLynxRawInfoReader, self).__init__(source, MassLynxBaseType.INFO)


    def GetNumberofFunctions( self ):
        size = c_int(0)
        getFunctionCount =  MassLynxRawReader.massLynxDll.getFunctionCount
        getFunctionCount.argtypes = [c_void_p, POINTER(c_int)]
        self.CheckReturnCode( getFunctionCount(self._getReader(),size) )

        return size.value

    def GetScansInFunction( self, whichFunction ):
        size = c_int(0)
        getScanCount = MassLynxRawReader.massLynxDll.getScanCount
        getScanCount.argtypes = [c_void_p, c_int, POINTER(c_int)]
        self.CheckReturnCode( getScanCount(self._getReader(),whichFunction,size) )

        return size.value

    def GetAcquisitionMassRange( self, whichFunction ):
        lowMass = c_float(0)
        highMass = c_float(0)
        getAcquisitionMassRange = MassLynxRawReader.massLynxDll.getAcquisitionMassRange
        getAcquisitionMassRange.argtypes = [c_void_p, c_int, c_int, POINTER(c_float), POINTER(c_float)]
        self.CheckReturnCode( getAcquisitionMassRange(self._getReader(),whichFunction, 0,lowMass,highMass) )

        return lowMass.value, highMass.value

    def GetAcquisitionTimeRange( self, whichFunction ):
        startTime = c_float(0)
        endTime = c_float(0)
        getAcquisitionTimeRange = MassLynxRawReader.massLynxDll.getAcquisitionTimeRange
        getAcquisitionTimeRange.argtypes = [c_void_p, c_int, POINTER(c_float), POINTER(c_float)]
        self.CheckReturnCode( getAcquisitionTimeRange(self._getReader(),whichFunction,startTime,endTime) )

        return startTime.value, endTime.value

    def GetFunctionType( self, whichFunction ):
        functionType = c_int(0)
        getFunctionType = MassLynxRawReader.massLynxDll.getFunctionType
        getFunctionType.argtypes = [c_void_p, c_int, POINTER(c_int)]
        self.CheckReturnCode( getFunctionType(self._getReader(),whichFunction, functionType) )

        return functionType.value

    def GetFunctionTypeString( self, functionType ):
        functionTypeString = c_char_p()
        temp = c_char_p()
        getFunctionTypeString = MassLynxRawReader.massLynxDll.getFunctionTypeString
        getFunctionTypeString.argtypes = [c_void_p, c_int, POINTER(c_char_p)]
        self.CheckReturnCode( getFunctionTypeString(self._getReader(),functionType, temp) )

        return self.ToString( temp )

    def IsContinuum( self, whichFunction ):
        continuum = c_bool(0)
        isContinuum = MassLynxRawReader.massLynxDll.isContinuum
        isContinuum.argtypes = [c_void_p, c_int, POINTER(c_bool)]
        self.CheckReturnCode( isContinuum(self._getReader(),whichFunction, continuum) )

        return continuum.value

    def GetIonMode( self, whichFunction ):
        ionMode = c_int()
        getIonMode = MassLynxRawReader.massLynxDll.getIonMode
        getIonMode.argtypes = [c_void_p, c_int, POINTER(c_int)]
        self.CheckReturnCode(getIonMode(self._getReader(),whichFunction, ionMode ))

        return ionMode.value

    def GetIonModeString( self, ionMode ):
        ionModeString = c_char_p()
        temp = c_char_p()
        getIonModeString = MassLynxRawReader.massLynxDll.getIonModeString
        getIonModeString.argtypes = [c_void_p, c_int, POINTER(c_char_p)]
        self.CheckReturnCode(getIonModeString(self._getReader(),ionMode,  temp ))
        return self.ToString(temp)

    def GetHeaderItems( self, whichItems ):
        itemString = c_char_p()
        nItems = len(whichItems )

        temp = c_char_p()
        items = (c_int * nItems)(*whichItems)
        delimiter = ctypes.create_string_buffer(1)

        getHeaderItem = MassLynxRawReader.massLynxDll.getHeaderItems
        getHeaderItem.argtypes = [c_void_p, POINTER(c_int), POINTER(c_char_p), c_int, POINTER(c_char)]
        self.CheckReturnCode( getHeaderItem( self._getReader(), items, temp, nItems, delimiter))

        itemString = self.ToString(temp)

        delim = delimiter.value.decode()
        return itemString.split(delim);

    def GetHeaderItemValue( self, whichItems ):
        nItems = len(whichItems )
        items = (c_int * nItems)(*whichItems)
        params = MassLynxParameters()
        getHeaderItemValue = MassLynxRawReader.massLynxDll.getHeaderItemValue
        getHeaderItemValue.argtypes = [c_void_p, POINTER(c_int), c_int, c_void_p]
        self.CheckReturnCode(getHeaderItemValue(self._getReader(), items, nItems, params.GetParameters()))

        return params

    def GetHeaderItem( self, whichItem ):
        whichItems =  list()
        whichItems.append( whichItem )
        values =  self.GetHeaderItems( whichItems)

        return values[0]

    # scan stats
    def GetScanItem( self, whichFunction, whichScan, whichItem ):
        whichItems = list()
        whichItems.append( whichItem )
        values =  self.GetScanItems(whichFunction, whichScan, whichItems)

        return values[0]

    def GetScanItems(self, whichFunction, whichScan, whichItems):
        nItems = len(whichItems)
        temp = c_char_p()
        items = (c_int * nItems)(*whichItems)
        delimiter = ctypes.create_string_buffer(1)

        getScanItem = MassLynxRawReader.massLynxDll.getScanItems
        getScanItem.argtypes = [c_void_p, c_int, c_int, POINTER(c_int), POINTER(c_char_p), c_int, POINTER(c_char)]
        self.CheckReturnCode( getScanItem( self._getReader(), whichFunction, whichScan, items, temp, nItems, delimiter))

        itemString = self.ToString(temp)

        delim = delimiter.value.decode()
        return itemString.split(delim);

    def GetScanItemString(self, whichItems):
         # get the array of items
        nItems = len(whichItems)
        temp = c_char_p()
        items = (c_int * nItems)(*whichItems)
        delimiter = ctypes.create_string_buffer(1)

        getScanItemNames = MassLynxRawReader.massLynxDll.getScanItemNames
        getScanItemNames.argtypes = [c_void_p, POINTER(c_int), POINTER(c_char_p), c_int, POINTER(c_char)]
        self.CheckReturnCode( getScanItemNames( self._getReader(), items, temp, nItems, delimiter))

        itemString = self.ToString(temp)

        delim = delimiter.value.decode()
        return itemString.split(delim);


    def GetItemsInFunction( self, whichFunction, nWhichScan ):
        size = c_int(0)
        pItems = c_void_p()
        getItemsInFunction = MassLynxRawReader.massLynxDll.getItemsInFunction
        getItemsInFunction.argtypes = [c_void_p, c_int, POINTER(c_void_p), POINTER(c_int)]
        self.CheckReturnCode( getItemsInFunction(self._getReader(),whichFunction,pItems,size) )

        # fill the array
        pI = cast(pItems,POINTER(c_int))
        items = pI[0:size.value]

        # dealocate memory
        MassLynxRawReader.ReleaseMemory( pItems)

        return items

    def GetRetentionTime( self, whichFunction, nWhichScan ):
        retentionTime = c_float(0)
        getRetentionTime = MassLynxRawReader.massLynxDll.getRetentionTime
        getRetentionTime.argtypes = [c_void_p, c_int, c_int, POINTER(c_float)]
        self.CheckReturnCode( getRetentionTime(self._getReader(),whichFunction,nWhichScan,retentionTime) )

        return retentionTime.value

    def GetDriftTime( self, whichFunction, nWhichDrift ):
        driftTime = c_float(0)
        getDriftTime = MassLynxRawReader.massLynxDll.getDriftTime
        getDriftTime.argtypes = [c_void_p, c_int, c_int, POINTER(c_float)]
        self.CheckReturnCode( getDriftTime(self._getReader(),whichFunction,nWhichDrift,driftTime) )

        return driftTime.value

    def GetDriftScanCount(self, whichFunction):
        scan_count = c_int(0)
        getDriftScanCount = MassLynxRawReader.massLynxDll.getDriftScanCount
        getDriftScanCount.argtypes = [c_void_p, c_int, POINTER(c_int)]
        self.CheckReturnCode(getDriftScanCount(self._getReader(), whichFunction, scan_count))
        return scan_count.value

    def GetCollisionalCrossSection(self, drift_time, mass, charge):
        ccs = c_float(0)
        getCollisionalCrossSection = MassLynxRawReader.massLynxDll.getCollisionalCrossSection
        getCollisionalCrossSection.argtypes = [
            c_void_p, c_float, c_float, c_int, POINTER(c_float)]
        self.CheckReturnCode(getCollisionalCrossSection(self._getReader(), drift_time, mass, charge, ccs))
        return ccs.value

    def CanLockMassCorrect( self ):
        canApply = c_bool(0)
        canLockMassCorrect =  MassLynxRawReader.massLynxDll.canLockMassCorrect
        canLockMassCorrect.argtypes = [c_void_p, POINTER(c_bool)]
        self.CheckReturnCode( canLockMassCorrect(self._getReader(), canApply) )

        return canApply.value


    def IsLockMassCorrected( self ):
        corrected = c_bool(0)
        mlMethod =  MassLynxRawReader.massLynxDll.isLockMassCorrected
        mlMethod.argtypes = [c_void_p, POINTER(c_bool)]
        self.CheckReturnCode( mlMethod(self._getReader(), corrected) )

        return corrected.value


