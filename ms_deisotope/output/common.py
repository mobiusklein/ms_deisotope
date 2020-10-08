# pragma: no cover
from ms_deisotope.utils import Base
from ms_deisotope.data_source.common import ScanBunch, ScanIterator


class ScanSerializerBase(object):
    """A common base class for all types which serialize :class:`~.Scan`-like objects to
    disk.

    These types support the context manager protocol, controlling their closing behavior.
    """

    def __init__(self, *args, **kwargs):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def save_scan_bunch(self, bunch, **kwargs):
        """Save a single :class:`~.ScanBunch` to the
        writing stream, recording any information relating
        the precursors and products to each other.

        May call :meth:`save_scan` multiple times.

        Parameters
        ----------
        bunch : :class:`~.ScanBunch`
            The scan bunch to save
        """
        self.save_scan(bunch.precursor, **kwargs)
        for prod in bunch.products:
            self.save_scan(prod, **kwargs)

    def save_scan(self, scan, **kwargs):
        """Save a single scan to the writing stream.

        Parameters
        ----------
        scan : :class:`ScanBase`
            The scan to save
        """
        raise NotImplementedError()

    def save(self, bunch, **kwargs):
        """Save any scan information in `bunch`.

        This method can handle :class:`~.ScanBunch` or :class:`ScanBase`
        instances, dispatching to the appropriate logic.

        Parameters
        ----------
        bunch : :class:`~.ScanBunch` or :class:`ScanBase`
            The scan data to save. May be a collection of related scans or a single scan.

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
        """Finish writing scan data, write any pending metadata
        and close the file stream.

        May call :meth:`complete`.
        """
        self.complete()

    def complete(self):
        """Perform any final operations after all spectra have been
        written, adding any trailing content necessary for the format.

        May be called from :meth:`close`, but before the writing stream
        has actually been closed.
        """
        pass


class ScanDeserializerBase(object):
    def __init__(self, *args, **kwargs):
        super(ScanDeserializerBase, self).__init__(*args, **kwargs)

    def __next__(self):
        return self.next()

    def next(self):
        raise NotImplementedError()

    def get_scan_by_id(self, scan_id):
        raise NotImplementedError()

    def get_scan_by_time(self, rt, require_ms1=False):
        raise NotImplementedError()

    def get_scan_by_index(self, index):
        raise NotImplementedError()


ScanIterator.register(ScanDeserializerBase)


class SampleRun(Base):

    def __init__(self, name, uuid, completed=True, sample_type=None, **kwargs):
        self.name = name
        self.uuid = uuid
        self.sample_type = sample_type
        self.completed = completed
        self.parameters = kwargs
