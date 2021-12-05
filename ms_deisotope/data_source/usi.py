from pyteomics.usi import USI, proxi, _PROXIBackend

from .scan import ScanDataSource, Scan, RawDataArrays, PrecursorInformation, ChargeNotProvided
from .metadata.scan_traits import IsolationWindow


USI_KEY = "_$USI"

class _PROXIScanSource(ScanDataSource):

    def _scan_arrays(self, scan):
        return RawDataArrays(scan.get("m/z array"), scan.get('intensity array'))

    def _scan_id(self, scan):
        if USI_KEY in scan:
            return scan[USI_KEY]
        if 'usi' in scan:
            return scan['usi']
        attributes = scan['attributes']
        for attrib in attributes:
            if attrib['name'] == 'scan number':
                return attrib['value']

    def _scan_index(self, scan):
        attributes = scan['attributes']
        for attrib in attributes:
            if attrib['name'] == 'scan number':
                return int(attrib['value']) - 1
        return -1

    def _scan_title(self, scan):
        if USI_KEY in scan:
            return scan[USI_KEY]
        if 'usi' in scan:
            return scan['usi']
        attributes = scan['attributes']
        for attrib in attributes:
            if attrib['name'] == 'scan number':
                return attrib['value']

    def _scan_time(self, scan):
        return None

    def _ms_level(self, scan):
        return 2

    def _precursor_information(self, scan):
        mz = 0.0
        charge = ChargeNotProvided
        attributes = scan['attributes']
        for attrib in attributes:
            if attrib['name'] == 'isolation window target m/z':
                mz = float(attrib['value'])
        return PrecursorInformation(mz, 0, charge, source=self)

    def _is_profile(self, scan):
        # TODO: Not guaranteed but no counter-examples
        return False

    def _polarity(self, scan):
        # TODO: Not guaranteed but no counter-examples
        return 1

    def _activation(self, scan):
        return None


class PROXIServiceClient(_PROXIScanSource):
    '''A :class:`~.ScanDataSource` client that provides random access to individual spectra
    by USI.

    Attributes
    ----------
    backend : str or :class:`pyteomics.usi._PROXIBackend`
        The server to fetch spectra from. If a name is given, the server will
        be inferred, else if a URL template is given, a new :class:`pyteomics.usi._PROXIBackend`
        is created from it.
    '''
    def __init__(self, backend='peptide_atlas'):
        self.backend = self._make_backend(backend)

    def _make_backend(self, backend):
        if isinstance(backend, _PROXIBackend):
            return backend
        if isinstance(backend, str):
            # We're given a URL
            if backend.startswith("http"):
                return _PROXIBackend(backend, backend)
        return backend

    def get_scan_by_id(self, scan_id):
        '''Get a spectrum by USI from the server.

        Parameters
        ----------
        scan_id : str
            The USI to fetch.


        Returns
        -------
        Scan
        '''
        payload = proxi(scan_id, self.backend)
        payload[USI_KEY] = scan_id
        return self._make_scan(payload)

    def get(self, usi):
        '''Get a spectrum by USI from the server.

        Parameters
        ----------
        scan_id : str
            The USI to fetch.


        Returns
        -------
        Scan
        '''
        return self.get_scan_by_id(usi)

    def __call__(self, usi):
        '''Get a spectrum by USI from the server.

        Parameters
        ----------
        scan_id : str
            The USI to fetch.


        Returns
        -------
        Scan
        '''
        return self.get_scan_by_id(usi)

    def __repr__(self):
        return "{self.__class__.__name__}({self.backend})".format(self=self)


peptide_atlas = PROXIServiceClient("peptide_atlas")
massive = PROXIServiceClient("massive")
pride = PROXIServiceClient("pride")

parse = USI.parse
