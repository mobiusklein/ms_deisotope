'''A common set of methods that are shared by
all :mod:`pyteomics`-based XML file readers.
'''

import warnings

from lxml import etree
from lxml.etree import XMLSyntaxError

from pyteomics import xml
from pyteomics.xml import unitfloat

from .common import (
    RandomAccessScanSource)
from ._compression import get_opener, test_if_file_has_fast_random_access


def in_minutes(x):
    '''Convert a time quantity to minutes

    Parameters
    ----------
    x: unitfloat
        A float representing a quantity of time annotated with a time unit

    Returns
    -------
    unitfloat:
        The time after conversion to minutes
    '''
    try:
        unit = x.unit_info
    except AttributeError:
        return x
    if unit == 'minute':
        return x
    elif unit == 'second':
        y = unitfloat(x / 60., 'minute')
        return y
    elif unit == 'hour':
        y = unitfloat(x * 60, 'minute')
        return y
    else:
        warnings.warn("Time unit %r not recognized" % unit)
    return x


class XMLReaderBase(RandomAccessScanSource):
    '''A common implementation of :mod:`pyteomics`-based XML file formats.

    Attributes
    ----------
    index: :class:`pyteomics.xml.ByteEncodingOrderedDict`
        The byte offset index used to achieve fast random access
    '''

    _parser_cls = None

    @property
    def has_fast_random_access(self):
        return test_if_file_has_fast_random_access(self.source.file)

    @classmethod
    def prebuild_byte_offset_file(cls, path):
        """Parse the file given by `path`, generating a byte offset index in
        JSON format and save it to disk for future use.

        This method is intended to provide a way to save time during repeated
        instantiation of this type over the same file by removing the need to
        do a full scan of the file to rebuild of the offset index each time.

        .. note::

            This assumes that `path` is either a path to a file in a directory
            which the invoking user has read and write access to, or that it is
            a file-like object whose `name` attribute gives a path that satisfies
            the same requirements.

        Parameters
        ----------
        path : :class:`str` or file-like
            The path to the file to index, or a file-like object with a name attribute.
        """
        return cls._parser_cls.prebuild_byte_offset_file(get_opener(path))

    @property
    def index(self):
        '''The byte offset index used to achieve fast random access.

        Maps :class:`~.ScanBase` IDs to the byte offsets, implying
        the order the scans reside in the file.

        Returns
        -------
        :class:`pyteomics.xml.ByteEncodingOrderedDict`
        '''
        return self._source._offset_index

    @property
    def source(self):
        '''The file parser that this reader consumes.
        '''
        return self._source

    @source.setter
    def source(self, value):
        self._source = value

    def close(self):
        '''Close the underlying reader.
        '''
        if self.source is not None:
            self.source.close()
            self._dispose()
            self.source = None
        super(XMLReaderBase, self).close()

    def __del__(self):
        try:
            self.close()
        except AttributeError:
            pass

    def __len__(self):
        return len(self.index)

    def reset(self):
        """Reset the object, clearing out any existing
        state.

        This resets the underlying file iterator, then
        calls :meth:`make_iterator`, and clears the scan
        cache.
        """
        self._source.reset()
        try:
            self.source.seek(0)
        except (IOError, AttributeError):
            pass
        self.make_iterator(None)
        self.initialize_scan_cache()

    def _validate(self, scan):
        raise NotImplementedError()

    def _make_default_iterator(self):
        return iter(self._source)

    def _dispose(self):
        super(XMLReaderBase, self)._dispose()
        try:
            self.index.clear()
        except Exception:
            pass

    def next(self):
        try:
            return next(self._producer)
        except XMLSyntaxError:
            raise StopIteration(
                "This iterator may need to be reset by calling `reset` to continue using it after"
                " using a random-access function like `get_by_id`")

    def __next__(self):
        return self.next()

    def get_scan_by_id(self, scan_id):
        """Retrieve the scan object for the specified scan id.

        If the scan object is still bound and in memory somewhere,
        a reference to that same object will be returned. Otherwise,
        a new object will be created.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        Scan
        """
        if isinstance(scan_id, bytes):
            scan_id = scan_id.decode('utf8')
        try:
            return self._scan_cache[scan_id]
        except KeyError:
            packed = self._make_scan(self._get_scan_by_id_raw(scan_id))
            self._scan_cache[packed.id] = packed
            return packed

    def _get_scan_by_id_raw(self, scan_id):
        try:
            return self._source.get_by_id(scan_id)
        except AttributeError as exc:
            # If the source is somehow None, then we should not try to convert
            # the error into a KeyError
            if self._source is None:
                raise exc
            # When index-driven lookup fails, pyteomics will parse the XML
            # file from the beginning until it runs out of file
            err = KeyError(scan_id)
            raise err

    def get_scan_by_time(self, time):
        """Retrieve the scan object for the specified scan time.

        This internally calls :meth:`get_scan_by_id` which will
        use its cache.

        Parameters
        ----------
        time : float
            The time to get the nearest scan from

        Returns
        -------
        Scan
        """
        scan_ids = tuple(self.index)
        lo = 0
        hi = len(scan_ids)

        best_match = None
        best_error = float('inf')

        if time == float('inf'):
            return self.get_scan_by_id(scan_ids[-1])

        while hi != lo:
            mid = (hi + lo) // 2
            sid = scan_ids[mid]
            scan = self.get_scan_by_id(sid)
            if not self._validate(scan):
                sid = scan_ids[mid - 1]
                scan = self.get_scan_by_id(sid)
                if not self._validate(scan):
                    sid = scan_ids[mid - 2]
                    scan = self.get_scan_by_id(sid)

            scan_time = scan.scan_time
            err = abs(scan_time - time)
            if err < best_error:
                best_error = err
                best_match = scan
            if scan_time == time:
                return scan
            elif (hi - lo) == 1:
                return best_match
            elif scan_time > time:
                hi = mid
            else:
                lo = mid
        if hi == 0 and not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")

    def get_scan_by_index(self, index):
        """Retrieve the scan object for the specified scan index.

        This internally calls :meth:`get_scan_by_id` which will
        use its cache.

        Parameters
        ----------
        index: int
            The index to get the scan for

        Returns
        -------
        Scan
        """
        if not self._use_index:
            raise TypeError("This method requires the index. Please pass `use_index=True` during initialization")
        # index_keys = tuple(self.index)
        # id_str = index_keys[index]
        # Use pyteomics index structure to avoid re-traversing the whole index again and again
        id_str = self.index.from_index(index)
        scan = self.get_scan_by_id(id_str)
        if not self._validate(scan):
            warnings.warn("index %d, id=%r does not appear to be a mass spectrum. Most behaviors will fail." % (
                index, id_str), stacklevel=2)
        return scan

    def _yield_from_index(self, scan_source, start):
        raise NotImplementedError()

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True, grouped=True):
        '''Reconstruct an iterator which will start from the scan matching one of ``scan_id``,
        ``rt``, or ``index``. Only one may be provided.

        After invoking this method, the iterator this object wraps will be changed to begin
        yielding scan bunchs (or single scans if ``grouped`` is ``False``).

        This method will trigger several random-access operations, making it prohibitively
        expensive for normally compressed files.

        Arguments
        ---------
        scan_id: str, optional
            Start from the scan with the specified id.
        rt: float, optional
            Start from the scan nearest to specified time (in minutes) in the run. If no
            exact match is found, the nearest scan time will be found, rounded up.
        index: int, optional
            Start from the scan with the specified index.
        require_ms1: bool, optional
            Whether the iterator must start from an MS1 scan. True by default.
        grouped: bool, optional
            whether the iterator should yield scan bunches or single scans. True by default.
        '''
        if scan_id is None:
            if rt is not None:
                scan = self.get_scan_by_time(rt)
            elif index is not None:
                try:
                    scan = self.get_scan_by_index(index)
                except IndexError:
                    if index > len(self.index):
                        index = len(self.index) - 1
                    else:
                        index = 0
                    scan = self.get_scan_by_index(index)

            else:
                raise ValueError("Must provide a scan locator, one of (scan_id, rt, index)")

            scan_id = scan.id
        else:
            scan = self.get_scan_by_id(scan_id)

        # We must start at an MS1 scan, so backtrack until we reach one
        if require_ms1:
            scan = self._locate_ms1_scan(scan)
            scan_id = scan.id

        iterator = self._yield_from_index(self._source, scan_id)
        self.make_iterator(iterator, grouped=grouped)
        return self

    def __repr__(self):
        return "{self.__class__.__name__}({self.source_file!r})".format(self=self) # pylint: disable=missing-format-attribute

    def __reduce__(self):
        return self.__class__, (self.source_file, self._use_index)


@xml._keepstate
def _find_section(source, section):
    value = next(source.iterfind(section))
    return value


@xml._keepstate
def get_tag_attributes(source, tag_name, stop_at=None):
    '''Iteratively parse XML stream in ``source`` until encountering ``tag_name``
    at which point parsing terminates and return the attributes of the matched
    tag.

    Parameters
    ----------
    source: file-like
        A file-like object over an XML document
    tag_name: str
        The name of the XML tag to parse until

    Returns
    -------
    dict
    '''
    g = etree.iterparse(source, ('start', 'end'))
    for event, tag in g:
        if event == 'start':
            if xml._local_name(tag) == tag_name:
                return tag.attrib
            elif stop_at and xml._local_name(tag) == stop_at:
                break
            else:
                continue
        else:
            tag.clear()
    return None


@xml._keepstate
def iterparse_until(source, target_name, quit_name):
    '''Iteratively parse XML stream in ``source``, yielding XML elements
    matching ``target_name``. If at any point a tag matching ``quit_name``
    is encountered, stop parsing.

    Parameters
    ----------
    source: file-like
        A file-like object over an XML document
    tag_name: str
        The name of the XML tag to parse until
    quit_name: str
        The name to stop parsing at.

    Yields
    ------
    lxml.etree.Element
    '''
    g = etree.iterparse(source, ('start', 'end'))
    for event, tag in g:
        if event == 'start':
            if xml._local_name(tag) == quit_name:
                break
            else:
                if xml._local_name(tag) == target_name:
                    yield tag
                else:
                    tag.clear()
