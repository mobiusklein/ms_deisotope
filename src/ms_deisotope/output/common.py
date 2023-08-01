# pragma: no cover
import os
import io

from typing import Any, Iterator, List, Optional, Tuple, Union, Hashable, Deque
import warnings

import numpy as np

from six import string_types as basestring
from ms_deisotope.data_source import ScanBase, Scan, ScanBunch, ScanIterator, ChargeNotProvided, PrecursorInformation, ProcessedRandomAccessScanSource
from ms_deisotope.data_source.metadata.sample import Sample
from ms_deisotope.peak_set import DeconvolutedPeak

from ms_deisotope.utils import Base, _typecast
from ms_deisotope.data_source.common import (
    _SingleScanIteratorImpl, _InterleavedGroupedScanIteratorImpl)
from ms_deisotope.data_source._compression import get_opener
from ms_deisotope.feature_map.scan_index import ExtendedScanIndex


class _ScanPackerOverridingMeta(type):

    def __new__(mcs, name, bases, namespace):
        if "get_scan_packer_type" not in namespace:
            @classmethod
            def get_scan_packer_type(cls):
                """
                Get a helper class that will pre-serialize parts of a spectrum so it can be
                transmitted between processes more easily at the cost of abstracting away the
                object model.
                """
                return None
            namespace['get_scan_packer_type'] = get_scan_packer_type
        return super().__new__(mcs, name, bases, namespace)


class ScanSerializerBase(object, metaclass=_ScanPackerOverridingMeta):
    """
    A common base class for all types which serialize :class:`~.Scan`-like objects to
    disk.

    These types support the context manager protocol, controlling their closing behavior.
    """

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def save_scan_bunch(self, bunch: ScanBunch, **kwargs):
        """
        Save a single :class:`~.ScanBunch` to the
        writing stream, recording any information relating
        the precursors and products to each other.

        May call :meth:`save_scan` multiple times.

        Parameters
        ----------
        bunch : :class:`~.ScanBunch`
            The scan bunch to save
        **kwargs
            Passed to :meth:`save_scan`
        """
        self.save_scan(bunch.precursor, **kwargs)
        for prod in bunch.products:
            self.save_scan(prod, **kwargs)

    def save_scan(self, scan: ScanBase, **kwargs):
        """
        Save a single scan to the writing stream.

        Parameters
        ----------
        scan : :class:`ScanBase`
            The scan to save
        **kwargs
            Extra information on how to save the scan
        """
        raise NotImplementedError()

    def save(self, bunch: Union[ScanBase, ScanBunch], **kwargs):
        """
        Save any scan information in `bunch`.

        This method can handle :class:`~.ScanBunch` or :class:`ScanBase`
        instances, dispatching to the appropriate logic.

        Parameters
        ----------
        bunch : :class:`~.ScanBunch` or :class:`ScanBase`
            The scan data to save. May be a collection of related scans or a single scan.
        **kwargs
            Extra information on how to save the scan

        See Also
        --------
        :meth:`save_scan`
        :meth:`save_scan_bunch`
        """
        if isinstance(bunch, ScanBunch):
            self.save_scan_bunch(bunch, **kwargs)
        else:
            self.save_scan(bunch, **kwargs)

    def close(self):
        """
        Finish writing scan data, write any pending metadata
        and close the file stream.

        May call :meth:`complete`.
        """
        self.complete()

    def complete(self):
        """
        Perform any final operations after all spectra have been
        written, adding any trailing content necessary for the format.

        May be called from :meth:`close`, but before the writing stream
        has actually been closed.
        """

    @classmethod
    def get_scan_packer_type(cls):
        """
        Get a helper class that will pre-serialize parts of a spectrum so it can be
        transmitted between processes more easily at the cost of abstracting away the
        object model.
        """
        return None


class ScanDeserializerBase(object):
    def __init__(self, *args, **kwargs):
        super(ScanDeserializerBase, self).__init__(*args, **kwargs)

    def __next__(self):
        return self.next() # pylint: disable=not-callable

    def next(self):
        raise NotImplementedError()

    def get_scan_by_id(self, scan_id):
        raise NotImplementedError()

    def get_scan_by_time(self, rt, require_ms1=False):
        raise NotImplementedError()

    def get_scan_by_index(self, index):
        raise NotImplementedError()


ScanIterator.register(ScanDeserializerBase)


class _ParamProp:
    __slots__ = ('key', )

    key: str

    def __init__(self, key: str):
        self.key = key

    def __get__(self, inst: Sample, cls):
        if inst is None:
            return self
        if self.key in inst.parameters:
            return inst.parameters[self.key]
        else:
            return None

    def __set__(self, inst: Sample, value: Any):
        inst.parameters[self.key] = _typecast(value)

    def __delete__(self, inst: Sample):
        del inst.parameters[self.key]


class SampleRun(Sample):

    uuid = _ParamProp('SampleRun-UUID')
    sample_type = _ParamProp('sample_type')
    completed = _ParamProp('completed')

    def __init__(self, name, uuid, completed=True, sample_type=None, **kwargs):
        super().__init__(id=name, name=name, **kwargs)
        self.uuid = uuid
        self.sample_type = sample_type
        self.completed = completed


class LCMSMSQueryInterfaceMixin(object):
    """A mixin class for querying an extended data index on a processed data file."""

    def require_extended_index(self) -> Optional[ExtendedScanIndex]:
        """
        Require the extended index is initialized.

        Either load the index if it exists, or build a new one.

        The new index is built within a the calling process, and may
        be considerably slower than the multiprocess index building
        tools.


        Returns
        -------
        ExtendedIndex
        """
        if not self.has_extended_index():
            try:
                if self.has_index_file():
                    self.read_index_file()
                else:
                    self.build_extended_index()
            except IOError:
                pass
            except ValueError:
                pass
            self._build_scan_id_to_rt_cache()
        return self.extended_index

    def has_extended_index(self) -> bool:
        """Check if the extended index is initialized."""
        return self.extended_index is not None

    def read_index_file(self, index_path=None):
        """
        Attempt to read the extended index from the file system.

        This updates :attr:`extended_index`.
        """
        if index_path is None:
            index_path = self._index_file_name
        with get_opener(index_path) as handle:
            if 'b' in handle.mode:
                handle = io.TextIOWrapper(handle, 'utf8')
            self.extended_index = ExtendedScanIndex.deserialize(handle)

    def has_index_file(self) -> bool:
        """
        Checks if an extended index file exists for this reader.

        Returns
        -------
        bool
        """
        try:
            return os.path.exists(self._index_file_name)
        except (TypeError, AttributeError):
            return False

    @property
    def _index_file_name(self) -> Optional[ExtendedScanIndex]:
        if isinstance(self.source_file, basestring):
            return ExtendedScanIndex.index_file_name(self.source_file)
        else:
            try:
                return ExtendedScanIndex.index_file_name(self.source_file.name)
            except AttributeError:
                return None

    def build_extended_index(self, header_only=True):
        """
        Build the extended index for this reader.

        This method builds the index in-process and may be slower
        than the multiprocessing builders.
        """
        self.reset()
        indexer = ExtendedScanIndex()
        iterator = self
        if header_only:
            iterator = self.iter_scan_headers()
        if self._has_ms1_scans():
            for bunch in iterator:
                indexer.add_scan_bunch(bunch)
        else:
            for scan in iterator:
                indexer.add_scan(scan)

        self.reset()
        self.extended_index = indexer
        try:
            with open(self._index_file_name, 'w') as handle:
                indexer.serialize(handle)
        except (IOError, OSError, AttributeError, TypeError) as err:
            warnings.warn(str(err))

    def get_index_information_by_scan_id(self, scan_id: str) -> dict:
        """Get the scan description from the extended index."""
        try:
            try:
                return self.extended_index.msn_ids[scan_id]
            except KeyError:
                return self.extended_index.ms1_ids[scan_id]
        except Exception:
            return {}

    def convert_scan_id_to_retention_time(self, scan_id: str) -> float:
        try:
            time = self._scan_id_to_rt[scan_id]
            return time
        except KeyError:
            header = self.get_scan_header_by_id(scan_id)
            return header.scan_time

    def _build_scan_id_to_rt_cache(self):
        if self.has_extended_index():
            for key in self.extended_index.ms1_ids:
                self._scan_id_to_rt[key] = self.extended_index.ms1_ids[
                    key]['scan_time']
            for key in self.extended_index.msn_ids:
                self._scan_id_to_rt[key] = self.extended_index.msn_ids[
                    key]['scan_time']

    def precursor_information(self) -> List[PrecursorInformation]:
        out = []
        for _, info in self.extended_index.msn_ids.items():
            mz = info['mz']
            prec_neutral_mass = info['neutral_mass']
            charge = info['charge']
            if charge == 'ChargeNotProvided':
                charge = ChargeNotProvided
            intensity = info['intensity']
            precursor_scan_id = info['precursor_scan_id']
            product_scan_id = info['product_scan_id']
            orphan = info.get('orphan', False)
            coisolation = info.get('coisolation', [])[:]
            defaulted = info.get('defaulted', False)
            pinfo = PrecursorInformation(
                mz, intensity, charge, precursor_scan_id,
                self, prec_neutral_mass, charge, intensity,
                product_scan_id=product_scan_id, orphan=orphan,
                coisolation=coisolation,
                defaulted=defaulted)
            if prec_neutral_mass is None or charge is None:
                continue
            out.append(pinfo)
        return out

    def ms1_peaks_above(self, mass_threshold: float = 500, intensity_threshold: float = 1000.) -> Deque[Tuple[str, DeconvolutedPeak, Hashable]]:
        """
        Read out all MS1 peaks from this spectrum source, if they satisfy the requested minimum values.

        Loop over the MS1 index, reading out peaks from the processed peak set whose
        mass exceeds ``mass_threshold`` and whose intensity exceeds ``intensity_threshold``.

        .. note::
            The final element of returned tuples is not guaranteed to be an integer, but only
            some hashable type.

        Returns
        -------
        peak_records : Deque[Tuple[:class:`str`, :class:`~.DeconvolutedPeak`, :class:`Hashable`]]
            A list of peaks associated with the MS1 scan ID they came from and an
            arbitrary hashable value to identify the pair.
        """
        accumulate = Deque()
        for ms1_id in self.extended_index.ms1_ids:
            scan = self.get_scan_by_id(ms1_id)
            for peak in scan.deconvoluted_peak_set:
                if peak.intensity < intensity_threshold or peak.neutral_mass < mass_threshold:
                    continue
                accumulate.append((ms1_id, peak, id(peak)))
        return accumulate

    def ms1_scan_times(self) -> np.ndarray:
        times = sorted(
            [bundle['scan_time'] for bundle in
             self.extended_index.ms1_ids.values()])
        return np.array(times)

    def extract_total_ion_current_chromatogram(self) -> np.ndarray:
        current = []
        for scan_id in self.extended_index.ms1_ids:
            header = self.get_scan_header_by_id(scan_id)
            current.append(header.arrays[1].sum())
        return np.array(current)

    def msms_for(self, query_mass, mass_error_tolerance=1e-5, start_time=None, end_time=None) -> List[PrecursorInformation]:
        out = []
        pinfos = self.extended_index.find_msms_by_precursor_mass(
            query_mass, mass_error_tolerance, bind=self)
        for pinfo in pinfos:
            valid = True
            product = None
            if start_time is not None or end_time is not None:
                product = self.get_scan_header_by_id(pinfo.product_scan_id)
            if start_time is not None and product.scan_time < start_time:
                valid = False
            elif end_time is not None and product.scan_time > end_time:
                valid = False
            if valid:
                out.append(pinfo)
        return out

    def _dispose(self):
        self._scan_id_to_rt.clear()
        self.extended_index.clear()
        return super(LCMSMSQueryInterfaceMixin, self)._dispose()

    def iter_scan_headers(self, iterator=None, grouped=True) -> Iterator[Scan]:
        try:
            if not self._has_ms1_scans():
                grouped = False
        except Exception:
            pass
        self.reset()
        if iterator is None:
            iterator = iter(self._source)

        _make_scan = super(LCMSMSQueryInterfaceMixin, self)._make_scan
        _validate = super(LCMSMSQueryInterfaceMixin, self)._validate

        if grouped:
            impl = _InterleavedGroupedScanIteratorImpl(
                iterator, _make_scan, _validate)
        else:
            impl = _SingleScanIteratorImpl(iterator, _make_scan, _validate)

        for x in impl:
            yield x

        self.reset()

    def get_scan_header_by_id(self, scan_id: str) -> Scan:
        """
        Retrieve the scan object for the specified scan id. If the
        scan object is still bound and in memory somewhere, a reference
        to that same object will be returned. Otherwise, a new object will
        be created.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        Scan
        """
        try:
            packed = super(LCMSMSQueryInterfaceMixin, self)._make_scan(
                self._source.get_by_id(scan_id))
            return packed
        except AttributeError as ae:
            raise AttributeError("Could not read attribute (%s) while looking up scan %s" % (
                ae, scan_id))

__all__ = [
    "ScanSerializerBase", "ScanDeserializerBase", "SampleRun", "LCMSMSQueryInterfaceMixin", "ProcessedRandomAccessScanSource"
]
