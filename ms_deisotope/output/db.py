from uuid import uuid4
import os

from sqlalchemy import create_engine, select, func, event
from sqlalchemy.orm import sessionmaker, scoped_session, validates, deferred
from sqlalchemy.orm.session import object_session
from sqlalchemy.engine import Connectable

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import (
    Column, Numeric, Integer, String, ForeignKey, PickleType,
    Boolean)
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.mutable import Mutable, MutableDict

from sqlalchemy import exc

import numpy as np

from ms_peak_picker import FittedPeak as MemoryFittedPeak, PeakIndex, PeakSet
from ms_deisotope import DeconvolutedPeak as MemoryDeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.peak_set import Envelope
from ms_deisotope.averagine import (
    mass_charge_ratio)

from ms_deisotope.data_source.common import (
    ProcessedScan, PrecursorInformation as MemoryPrecursorInformation, ScanBunch)

from ms_deisotope.output.common import ScanSerializerBase, ScanDeserializerBase


def Mass(index=True):
    return Column(Numeric(14, 6, asdecimal=False), index=index)


class MutableList(Mutable, list):

    @classmethod
    def coerce(cls, key, value):
        if not isinstance(value, MutableList):
            if isinstance(value, list):
                return MutableList(value)
            value = Mutable.coerce(key, value)

        return value

    def __setitem__(self, key, value):
        old_value = list.__getitem__(self, key)
        for obj, key in self._parents.items():
            old_value._parents.pop(obj, None)

        list.__setitem__(self, key, value)
        for obj, key in self._parents.items():
            value._parents[obj] = key

        self.changed()

    def __getstate__(self):
        return list(self)

    def __setstate__(self, state):
        self[:] = state


Base = declarative_base()


def find_by_name(session, model_class, name):
    return session.query(model_class).filter(model_class.name == name).first()


def make_unique_name(session, model_class, name):
    marked_name = name
    i = 1
    while find_by_name(session, model_class, marked_name) is not None:
        marked_name = "%s (%d)" % (name, i)
        i += 1
    return marked_name


class HasUniqueName(object):
    name = Column(String(128), default=u"", unique=True)
    uuid = Column(String(64), index=True, unique=True)

    @classmethod
    def make_unique_name(cls, session, name):
        return make_unique_name(session, cls, name)

    @classmethod
    def find_by_name(cls, session, name):
        return find_by_name(session, cls, name)

    @validates("name")
    def ensure_unique_name(self, key, name):
        session = object_session(self)
        if session is not None:
            model_class = self.__class__
            name = make_unique_name(session, model_class, name)
            return name
        else:
            return name


class SampleRun(Base, HasUniqueName):
    __tablename__ = "SampleRun"

    id = Column(Integer, primary_key=True, autoincrement=True)

    ms_scans = relationship(
        "MSScan", backref=backref("sample_run"), lazy='dynamic')
    sample_type = Column(String(128))

    completed = Column(Boolean(), default=False, nullable=False)

    def save_bunch(self, bunch):
        session = object_session(self)
        return serialize_scan_bunch(session, bunch, self.id)

    def __repr__(self):
        return "SampleRun(id=%d, name=%s)" % (self.id, self.name)


class MSScan(Base):
    __tablename__ = "MSScan"

    id = Column(Integer, primary_key=True, autoincrement=True)
    index = Column(Integer, index=True)
    ms_level = Column(Integer)
    scan_time = Column(Numeric(10, 5, asdecimal=False), index=True)
    title = Column(String(512))
    scan_id = Column(String(512), index=True)
    sample_run_id = Column(Integer, ForeignKey(
        SampleRun.id, ondelete='CASCADE'), index=True)

    peak_set = relationship("FittedPeak", backref="scan", lazy="dynamic")
    deconvoluted_peak_set = relationship(
        "DeconvolutedPeak", backref='scan', lazy='dynamic')

    info = deferred(Column(MutableDict.as_mutable(PickleType)))

    def __repr__(self):
        f = "{}({}, {}, {}, {}".format(
            self.__class__.__name__, self.scan_id, self.ms_level, self.scan_time,
            self.deconvoluted_peak_set.count())
        if self.ms_level > 1:
            f = "%s %s" % (f, self.precursor_information)
        f += ")"
        return f

    def convert(self, fitted=True, deconvoluted=True):
        precursor_information = self.precursor_information.convert(
        ) if self.precursor_information is not None else None

        session = object_session(self)
        conn = session.connection()

        if fitted:
            q = conn.execute(select([FittedPeak.__table__]).where(
                FittedPeak.__table__.c.scan_id == self.id)).fetchall()

            peak_set_items = list(
                map(make_memory_fitted_peak, q))

            peak_set = PeakSet(peak_set_items)
            peak_set._index()
            peak_index = PeakIndex(np.array([], dtype=np.float64), np.array(
                [], dtype=np.float64), peak_set)
        else:
            peak_index = PeakIndex(np.array([], dtype=np.float64), np.array(
                [], dtype=np.float64), PeakSet([]))

        if deconvoluted:
            q = conn.execute(select([DeconvolutedPeak.__table__]).where(
                DeconvolutedPeak.__table__.c.scan_id == self.id)).fetchall()

            deconvoluted_peak_set_items = list(
                map(make_memory_deconvoluted_peak, q))

            deconvoluted_peak_set = DeconvolutedPeakSet(
                deconvoluted_peak_set_items)
            deconvoluted_peak_set.reindex()
        else:
            deconvoluted_peak_set = DeconvolutedPeakSet([])

        info = self.info or {}

        scan = ProcessedScan(
            self.scan_id, self.title, precursor_information, int(self.ms_level),
            float(self.scan_time), self.index, peak_index, deconvoluted_peak_set,
            activation=info.get('activation'))
        return scan

    @classmethod
    def _serialize_scan(cls, scan, sample_run_id=None):
        db_scan = cls(
            index=scan.index, ms_level=scan.ms_level,
            scan_time=float(scan.scan_time), title=scan.title,
            scan_id=scan.id, sample_run_id=sample_run_id,
            info={'activation': scan.activation})
        return db_scan

    @classmethod
    def serialize(cls, scan, sample_run_id=None):
        db_scan = cls._serialize_scan(scan, sample_run_id)
        db_scan.peak_set = map(FittedPeak.serialize, scan.peak_set)
        db_scan.deconvoluted_peak_set = map(
            DeconvolutedPeak.serialize, scan.deconvoluted_peak_set)
        return db_scan

    @classmethod
    def serialize_bulk(cls, scan, sample_run_id, session, fitted=True, deconvoluted=True):
        db_scan = cls._serialize_scan(scan, sample_run_id)

        session.add(db_scan)
        session.flush()

        if fitted:
            FittedPeak._serialize_bulk_list(scan.peak_set, db_scan.id, session)
        if deconvoluted:
            DeconvolutedPeak._serialize_bulk_list(
                scan.deconvoluted_peak_set, db_scan.id, session)
        return db_scan


class PrecursorInformation(Base):
    __tablename__ = "PrecursorInformation"

    id = Column(Integer, primary_key=True)
    sample_run_id = Column(Integer, ForeignKey(
        SampleRun.id, ondelete='CASCADE'), index=True)

    precursor_id = Column(Integer, ForeignKey(
        MSScan.id, ondelete='CASCADE'), index=True)
    precursor = relationship(MSScan, backref=backref("product_information"),
                             primaryjoin="PrecursorInformation.precursor_id == MSScan.id",
                             uselist=False)
    product_id = Column(Integer, ForeignKey(
        MSScan.id, ondelete='CASCADE'), index=True)
    product = relationship(MSScan, backref=backref("precursor_information", uselist=False),
                           primaryjoin="PrecursorInformation.product_id == MSScan.id",
                           uselist=False)
    neutral_mass = Mass()
    charge = Column(Integer)
    intensity = Column(Numeric(16, 4, asdecimal=False))
    defaulted = Column(Boolean)
    orphan = Column(Boolean)

    @property
    def extracted_neutral_mass(self):
        return self.neutral_mass

    @property
    def extracted_charge(self):
        return self.charge

    @property
    def extracted_intensity(self):
        return self.intensity

    def __repr__(self):
        return "DBPrecursorInformation({}, {}, {})".format(
            self.precursor.scan_id, self.neutral_mass, self.charge)

    def convert(self, data_source=None):
        return MemoryPrecursorInformation(
            mass_charge_ratio(self.neutral_mass,
                              self.charge), self.intensity, self.charge,
            self.precursor.scan_id, data_source, self.neutral_mass, self.charge,
            self.intensity)

    @classmethod
    def serialize(cls, inst, precursor, product, sample_run_id):
        db_pi = PrecursorInformation(
            precursor_id=precursor.id, product_id=product.id,
            charge=inst.extracted_charge, intensity=inst.extracted_intensity,
            neutral_mass=inst.extracted_neutral_mass, sample_run_id=sample_run_id,
            defaulted=inst.defaulted, orphan=inst.orphan)
        return db_pi


class PeakMixin(object):
    mz = Mass(False)
    intensity = Column(Numeric(16, 4, asdecimal=False))
    full_width_at_half_max = Column(Numeric(7, 6, asdecimal=False))
    signal_to_noise = Column(Numeric(10, 3, asdecimal=False))
    area = Column(Numeric(12, 4, asdecimal=False))

    @classmethod
    def _serialize_bulk_list(cls, peaks, scan_id, session):
        out = cls._prepare_serialize_list(peaks, scan_id)
        session.bulk_save_objects(out)

    @classmethod
    def _prepare_serialize_list(cls, peaks, scan_id):
        out = []
        for peak in peaks:
            db_peak = cls.serialize(peak)
            db_peak.scan_id = scan_id
            out.append(db_peak)
        return out

    def __repr__(self):
        return "DB" + repr(self.convert())


def make_memory_fitted_peak(self):
    return MemoryFittedPeak(
        self.mz, self.intensity, self.signal_to_noise, -1, -1,
        self.full_width_at_half_max, (self.area if self.area is not None else 0.))


class FittedPeak(Base, PeakMixin):
    __tablename__ = "FittedPeak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    scan_id = Column(Integer, ForeignKey(
        MSScan.id, ondelete='CASCADE'), index=True)

    convert = make_memory_fitted_peak

    @classmethod
    def serialize(cls, peak):
        return cls(
            mz=peak.mz, intensity=peak.intensity, signal_to_noise=peak.signal_to_noise,
            full_width_at_half_max=peak.full_width_at_half_max, area=peak.area)


def make_memory_deconvoluted_peak(self):
    return MemoryDeconvolutedPeak(
        self.neutral_mass, self.intensity, self.charge,
        self.signal_to_noise, -1, self.full_width_at_half_max,
        self.a_to_a2_ratio, self.most_abundant_mass, self.average_mass,
        self.score, Envelope(
            self.envelope), self.mz, None, self.chosen_for_msms,
        (self.area if self.area is not None else 0.))


class DeconvolutedPeak(Base, PeakMixin):
    __tablename__ = "DeconvolutedPeak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    neutral_mass = Mass(False)
    average_mass = Mass(False)
    most_abundant_mass = Mass(False)

    charge = Column(Integer)

    score = Column(Numeric(12, 6, asdecimal=False))
    scan_id = Column(Integer, ForeignKey(
        MSScan.id, ondelete='CASCADE'), index=True)
    envelope = Column(MutableList.as_mutable(PickleType))
    a_to_a2_ratio = Column(Numeric(8, 7, asdecimal=False))
    chosen_for_msms = Column(Boolean)

    def __eq__(self, other):
        return (abs(self.neutral_mass - other.neutral_mass) < 1e-5) and (
            abs(self.intensity - other.intensity) < 1e-5)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return hash((self.mz, self.intensity, self.charge))

    convert = make_memory_deconvoluted_peak

    @classmethod
    def serialize(cls, peak):
        return cls(
            mz=peak.mz, intensity=peak.intensity, signal_to_noise=peak.signal_to_noise,
            full_width_at_half_max=peak.full_width_at_half_max,
            neutral_mass=peak.neutral_mass, average_mass=peak.average_mass,
            most_abundant_mass=peak.most_abundant_mass, charge=peak.charge, score=peak.score,
            envelope=list(peak.envelope), a_to_a2_ratio=peak.a_to_a2_ratio,
            chosen_for_msms=peak.chosen_for_msms, area=(peak.area if peak.area is not None else 0.))


def serialize_scan_bunch(session, bunch, sample_run_id=None):
    precursor = bunch.precursor
    db_precursor = MSScan.serialize(precursor, sample_run_id=sample_run_id)
    session.add(db_precursor)
    db_products = [MSScan.serialize(p, sample_run_id=sample_run_id)
                   for p in bunch.products]
    session.add_all(db_products)
    session.flush()
    for scan, db_scan in zip(bunch.products, db_products):
        pi = scan.precursor_information
        db_pi = PrecursorInformation.serialize(
            pi, db_precursor, db_scan, sample_run_id)
        session.add(db_pi)
    session.flush()
    return db_precursor, db_products


def serialize_scan_bunch_bulk(session, bunch, sample_run_id, ms1_fitted=True, msn_fitted=True):
    precursor = bunch.precursor
    db_precursor = MSScan.serialize_bulk(
        precursor, sample_run_id, session, fitted=ms1_fitted)
    db_products = [MSScan.serialize_bulk(p, sample_run_id, session, fitted=msn_fitted)
                   for p in bunch.products]
    for scan, db_scan in zip(bunch.products, db_products):
        pi = scan.precursor_information
        db_pi = PrecursorInformation(
            precursor_id=db_precursor.id, product_id=db_scan.id,
            charge=pi.extracted_charge, intensity=pi.extracted_intensity,
            neutral_mass=pi.extracted_neutral_mass, sample_run_id=sample_run_id)
        session.add(db_pi)
    session.flush()
    return db_precursor, db_products


class ScanSerializationBatcher(object):
    """Encapsulates the process of saving ScanBunch objects to the database
    incrementally, delaying writing peak information until a large collection
    of peaks are waiting to be written in a single burst to minimize overhead.

    This object assumes it can reuse the same Session over and over again.

    Attributes
    ----------
    batch_size : int
        The number of ScanBunch objects to try to save before
        flushing all peaks to disk.
    count : int
        The number of ScanBunch objects saved since the last
        flush.
    deconvoluted_list : list
        Accumulator of db.DeconvolutedPeak objects whose :attr:`scan_id`
        attribute has been set to their parent db.MSScan's primary key
        awaiting being sent to the database.
    fitted_list : list
        Accumulator of db.FittedPeak objects whose :attr:`scan_id`
        attribute has been set to their parent db.MSScan's primary key
        awaiting being sent to the database.
    ms1_fitted : bool
        Whether or not to save FittedPeak objects attatched to scans where
        :attr:`ms_level` == `1`
    msn_fitted : bool
        Whether or not to save FittedPeak objects attatched to scans where
        :attr:`ms_level` > `1`
    sample_run_id : int
        The primary key of the db.SampleRun that all created db.MSScan objects
        should belong to.
    session : Session
        The `sqlalchemy` session object that this object uses to communicate
        with the database
    """
    def __init__(self, session, sample_run_id, batch_size=50, ms1_fitted=True, msn_fitted=True):
        self.session = session
        self.sample_run_id = sample_run_id
        self.batch_size = batch_size
        self.ms1_fitted = ms1_fitted
        self.msn_fitted = msn_fitted
        self.fitted_list = []
        self.deconvoluted_list = []
        self.count = 0

    def add(self, bunch):
        """Save the incoming ScanBunch object partially, sending all contained Scan
        and PrecursorInformation to the database, but only mapping Peak objects to
        their equivalent ORM classes and populating the relevant foreign key attributes
        and accumulating those for later flushing.

        This function invokes :meth:`check_condition` which may cause peaks to be sent
        to the database and have the session flushed.

        Parameters
        ----------
        bunch : ScanBunch
            The set of related scans to be saved
        """
        precursor = bunch.precursor
        db_precursor = MSScan._serialize_scan(precursor, self.sample_run_id)
        db_products = [MSScan._serialize_scan(
            prod, self.sample_run_id) for prod in bunch.products]
        self.session.add(db_precursor)
        self.session.add_all(db_products)
        self.session.flush()
        if self.ms1_fitted:
            self.fitted_list.extend(
                FittedPeak._prepare_serialize_list(precursor.peak_set, db_precursor.id))
        self.deconvoluted_list.extend(
            DeconvolutedPeak._prepare_serialize_list(
                precursor.deconvoluted_peak_set, db_precursor.id))
        for scan, db_scan in zip(bunch.products, db_products):
            pi = scan.precursor_information
            db_pi = PrecursorInformation(
                precursor_id=db_precursor.id, product_id=db_scan.id,
                charge=pi.extracted_charge, intensity=pi.extracted_intensity,
                neutral_mass=pi.extracted_neutral_mass, sample_run_id=self.sample_run_id)
            self.session.add(db_pi)
            if self.msn_fitted:
                self.fitted_list.extend(
                    FittedPeak._prepare_serialize_list(scan.peak_set, db_scan.id))
            self.deconvoluted_list.extend(
                DeconvolutedPeak._prepare_serialize_list(
                    scan.deconvoluted_peak_set, db_scan.id))
        self.count += 1
        self.session.flush()
        self.check_condition()

    def check_condition(self):
        """Checks to see if the batch is done,
        and if the batch is done, send the accumulated
        peaks to the database and clear the accumulators
        """
        if self.count >= self.batch_size:
            self.serialize_peaks()

    def serialize_peaks(self):
        """Sends accumulated peaks to the database and flushes
        :attr:`session`. Also resets the accumulator lists and
        sets :attr:`count` to `0`
        """
        if self.fitted_list:
            self.session.bulk_save_objects(self.fitted_list)
            self.fitted_list = []
        self.session.bulk_save_objects(self.deconvoluted_list)
        self.deconvoluted_list = []
        self.session.flush()
        self.count = 0


def configure_connection(connection, create_tables=True):
    if isinstance(connection, basestring):
        try:
            eng = create_engine(connection)
        except exc.ArgumentError:
            eng = SQLiteConnectionRecipe(connection)()
    elif isinstance(connection, Connectable):
        eng = connection
    elif isinstance(connection, ConnectionRecipe):
        eng = connection()
    elif isinstance(connection, scoped_session):
        eng = connection.get_bind()
    else:
        raise ValueError(
            "Could not determine how to get a database connection from %r" % connection)
    if create_tables:
        Base.metadata.create_all(bind=eng)
    return eng


initialize = configure_connection


class ConnectionRecipe(object):

    def __init__(self, connection_url, connect_args=None, on_connect=None, **engine_args):
        if connect_args is None:
            connect_args = {}
        if on_connect is None:
            def on_connect(connection, connection_record):
                pass

        self.connection_url = connection_url
        self.connect_args = connect_args
        self.on_connect = on_connect
        self.engine_args = engine_args

    def __call__(self):
        connection = create_engine(
            self.connection_url, connect_args=self.connect_args,
            **self.engine_args)
        event.listens_for(connection, 'connect')(self.on_connect)
        return connection


class SQLiteConnectionRecipe(ConnectionRecipe):
    connect_args = {
        'timeout': 180,
    }

    @staticmethod
    def _configure_connection(dbapi_connection, connection_record):
        # disable pysqlite's emitting of the BEGIN statement entirely.
        # also stops it from emitting COMMIT before any DDL.
        iso_level = dbapi_connection.isolation_level
        dbapi_connection.isolation_level = None
        try:
            dbapi_connection.execute("PRAGMA page_size = 5120;")
            dbapi_connection.execute("PRAGMA cache_size = 12000;")
            dbapi_connection.execute("PRAGMA temp_store = 3;")
            if not int(os.environ.get("NOWAL", '0')):
                dbapi_connection.execute("PRAGMA journal_mode = WAL;")
                dbapi_connection.execute("PRAGMA wal_autocheckpoint = 100;")
            # dbapi_connection.execute("PRAGMA foreign_keys = ON;")
            dbapi_connection.execute("PRAGMA journal_size_limit = 1000000;")
            pass
        except Exception as e:
            print(e)
        dbapi_connection.isolation_level = iso_level

    def __init__(self, connection_url, **engine_args):
        super(SQLiteConnectionRecipe, self).__init__(
            self._construct_url(connection_url), self.connect_args, self._configure_connection)

    def _construct_url(self, path):
        if path.startswith("sqlite://"):
            return path
        else:
            return "sqlite:///%s" % path


class DatabaseBoundOperation(object):

    def __init__(self, connection):
        self.engine = self._configure_connection(connection)

        self._original_connection = connection
        self._sessionmaker = scoped_session(
            sessionmaker(bind=self.engine, autoflush=False))
        self._session = None

    def bridge(self, other):
        self.engine = other.engine
        self._sessionmaker = other._sessionmaker
        self._session = other._session

    def _configure_connection(self, connection):
        eng = configure_connection(connection, create_tables=True)
        return eng

    def __eq__(self, other):
        return str(self.engine) == str(other.engine) and\
            (self.__class__ is other.__class__)

    def __hash__(self):
        return hash(str(self.engine))

    def __ne__(self, other):
        return not (self == other)

    def __reduce__(self):
        return self.__class__, (self._original_connection,)

    @property
    def session(self):
        return self._sessionmaker

    def query(self, *args):
        return self.session.query(*args)

    def _analyze_database(self):
        conn = self.engine.connect()
        inner = conn.connection.connection
        cur = inner.cursor()
        cur.execute("analyze;")
        inner.commit()

    def _sqlite_reload_analysis_plan(self):
        conn = self.engine.connect()
        conn.execute("ANALYZE sqlite_master;")

    def _sqlite_checkpoint_wal(self):
        conn = self.engine.connect()
        inner = conn.connection.connection
        cur = inner.cursor()
        cur.execute("PRAGMA wal_checkpoint(SQLITE_CHECKPOINT_RESTART);")

    @property
    def dialect(self):
        return self.engine.dialect

    def is_sqlite(self):
        return self.dialect.name == "sqlite"

    def close(self):
        self.session.remove()
        self.engine.dispose()


class DatabaseScanSerializer(ScanSerializerBase, DatabaseBoundOperation):

    def __init__(self, connection, sample_name=None, overwrite=True, save_fitted=False):
        self.uuid = str(uuid4())

        if sample_name is None:
            sample_name = self.uuid

        self._sample_name = sample_name
        DatabaseBoundOperation.__init__(self, connection)
        self._sample_run = None
        self._sample_run_id = None
        self.save_fitted = save_fitted

        self._overwrite(overwrite)

    def _overwrite(self, overwrite=True):
        try:
            sample_run = self.sample_run
        except exc.IntegrityError:
            self.session.rollback()
            sample_run = self.session.query(SampleRun).filter(
                SampleRun.name == self.sample_name).first()
            if not sample_run.completed or overwrite:
                self.session.delete(sample_run)
                self.session.commit()
                self._construct_sample_run()
            else:
                raise

    def __reduce__(self):
        return self.__class__, (self._original_connection, self.sample_name, False)

    def _construct_sample_run(self):
        sr = SampleRun(uuid=self.uuid, name=self._sample_name)
        self.session.add(sr)
        self.session.commit()
        self._sample_run = sr
        self._sample_run_id = sr.id
        self._sample_name = sr.name

    def complete(self):
        self.sample_run.completed = True
        self.session.add(self.sample_run)
        self.session.commit()

    @property
    def sample_run_id(self):
        if self._sample_run_id is None:
            self._construct_sample_run()

        return self._sample_run_id

    @property
    def sample_run(self):
        if self._sample_run is None:
            self._construct_sample_run()
        return self._sample_run

    @property
    def sample_name(self):
        if self._sample_run is None:
            self._construct_sample_run()
        return self._sample_name

    def save(self, bunch, commit=True, bulk=True):
        if bulk:
            out = serialize_scan_bunch_bulk(
                self.session, bunch, self.sample_run_id,
                ms1_fitted=self.save_fitted, msn_fitted=self.save_fitted)
        else:
            out = serialize_scan_bunch(self.session, bunch, self.sample_run_id)
        if commit:
            self.session.commit()
        return out

    def commit(self):
        self.session.commit()


class BatchingDatabaseScanSerializer(DatabaseScanSerializer):

    def __init__(self, connection, sample_name=None, overwrite=True, save_fitted=False, batch_size=10):
        super(BatchingDatabaseScanSerializer, self).__init__(
            connection, sample_name, overwrite, save_fitted)
        self._batch = None
        self.batch_size = batch_size

    @property
    def batch(self):
        if self._batch is None:
            self._batch = ScanSerializationBatcher(
                self.session, self.sample_run_id, self.batch_size, self.save_fitted, self.save_fitted)
        return self._batch

    def save(self, bunch, commit=True):
        self.batch.add(bunch)

    def commit(self):
        self.batch.serialize_peaks()
        self.session.commit()

    def complete(self):
        self.batch.serialize_peaks()
        super(BatchingDatabaseScanSerializer, self).complete()


def flatten(iterable):
    return [y for x in iterable for y in x]


class DatabaseScanDeserializer(ScanDeserializerBase, DatabaseBoundOperation):

    def __init__(self, connection, sample_name=None, sample_run_id=None):

        DatabaseBoundOperation.__init__(self, connection)

        self._sample_run = None
        self._sample_name = sample_name
        self._sample_run_id = sample_run_id
        self._iterator = None
        self._scan_id_to_retention_time_cache = None

    def _intialize_scan_id_to_retention_time_cache(self):
        self._scan_id_to_retention_time_cache = dict(
            self.session.query(MSScan.scan_id, MSScan.scan_time).filter(
                MSScan.sample_run_id == self.sample_run_id))

    def __reduce__(self):
        return self.__class__, (
            self._original_connection, self.sample_name, self.sample_run_id)

    # Sample Run Bound Handle API

    @property
    def sample_run_id(self):
        if self._sample_run_id is None:
            self._retrieve_sample_run()
        return self._sample_run_id

    @property
    def sample_run(self):
        if self._sample_run is None:
            self._retrieve_sample_run()
        return self._sample_run

    @property
    def sample_name(self):
        if self._sample_name is None:
            self._retrieve_sample_run()
        return self._sample_name

    def _retrieve_sample_run(self):
        session = self.session
        if self._sample_name is not None:
            sr = session.query(SampleRun).filter(
                SampleRun.name == self._sample_name).one()
        elif self._sample_run_id is not None:
            sr = session.query(SampleRun).filter(
                SampleRun.id == self._sample_run_id).one()
        else:
            sr = session.query(SampleRun).first()
        self._sample_run = sr
        self._sample_run_id = sr.id
        self._sample_name = sr.name

    # Scan Generator Public API

    def get_scan_by_id(self, scan_id):
        q = self._get_by_scan_id(scan_id)
        if q is None:
            raise KeyError(scan_id)
        mem = q.convert()
        if mem.precursor_information:
            mem.precursor_information.source = self
        return mem

    def convert_scan_id_to_retention_time(self, scan_id):
        if self._scan_id_to_retention_time_cache is None:
            self._intialize_scan_id_to_retention_time_cache()
        try:
            return self._scan_id_to_retention_time_cache[scan_id]
        except KeyError:
            q = self.session.query(MSScan.scan_time).filter(
                MSScan.scan_id == scan_id, MSScan.sample_run_id == self.sample_run_id).scalar()
            self._scan_id_to_retention_time_cache[scan_id] = q
            return q

    def _select_index(self, require_ms1=True):
        indices_q = self.session.query(MSScan.index).filter(
            MSScan.sample_run_id == self.sample_run_id).order_by(MSScan.index.asc())
        if require_ms1:
            indices_q = indices_q.filter(MSScan.ms_level == 1)
        indices = flatten(indices_q.all())
        return indices

    def _iterate_over_index(self, start=0, require_ms1=True):
        indices = self._select_index(require_ms1)
        try:
            i = indices.index(start)
        except ValueError:
            lo = 0
            hi = len(indices)

            while lo != hi:
                mid = (lo + hi) / 2
                x = indices[mid]
                if x == start:
                    i = mid
                    break
                elif lo == (hi - 1):
                    i = mid
                    break
                elif x > start:
                    hi = mid
                else:
                    lo = mid
        items = indices[i:]
        i = 0
        n = len(items)
        while i < n:
            index = items[i]
            scan = self.session.query(MSScan).filter(
                MSScan.index == index,
                MSScan.sample_run_id == self.sample_run_id).one()
            products = [pi.product for pi in scan.product_information]
            yield ScanBunch(scan.convert(), [p.convert() for p in products])
            i += 1

    def __iter__(self):
        return self

    def next(self):
        if self._iterator is None:
            self._iterator = self._iterate_over_index()
        return next(self._iterator)

    def _get_by_scan_id(self, scan_id):
        q = self.session.query(MSScan).filter(
            MSScan.scan_id == scan_id, MSScan.sample_run_id == self.sample_run_id).first()
        return q

    def _get_scan_by_time(self, rt, require_ms1=False):
        times_q = self.session.query(MSScan.scan_time).filter(
            MSScan.sample_run_id == self.sample_run_id).order_by(MSScan.scan_time.asc())
        if require_ms1:
            times_q = times_q.filter(MSScan.ms_level == 1)
        times = flatten(times_q.all())
        try:
            i = times.index(rt)
        except ValueError:
            lo = 0
            hi = len(times)

            while lo != hi:
                mid = (lo + hi) / 2
                x = times[mid]
                if x == rt:
                    i = mid
                    break
                elif lo == (hi - 1):
                    i = mid
                    break
                elif x > rt:
                    hi = mid
                else:
                    lo = mid
        scan = self.session.query(MSScan).filter(
            MSScan.scan_time == times[i],
            MSScan.sample_run_id == self.sample_run_id).one()
        return scan

    def reset(self):
        self._iterator = None

    def get_scan_by_time(self, rt, require_ms1=False):
        q = self._get_scan_by_time(rt, require_ms1)
        mem = q.convert()
        if mem.precursor_information:
            mem.precursor_information.source = self
        return mem

    def _get_scan_by_index(self, index):
        q = self.session.query(MSScan).filter(
            MSScan.index == index, MSScan.sample_run_id == self.sample_run_id).one()
        return q

    def get_scan_by_index(self, index):
        mem = self._get_scan_by_index(index).convert()
        if mem.precursor_information:
            mem.precursor_information.source = self
        return mem

    def _locate_ms1_scan(self, scan):
        while scan.ms_level != 1:
            scan = self._get_scan_by_index(scan.index - 1)
        return scan

    def start_from_scan(self, scan_id=None, rt=None, index=None, require_ms1=True):
        if scan_id is None:
            if rt is not None:
                scan = self._get_scan_by_time(rt)
            elif index is not None:
                scan = self._get_scan_by_index(index)
        else:
            scan = self._get_by_scan_id(scan_id)

        # We must start at an MS1 scan, so backtrack until we reach one
        if require_ms1:
            scan = self._locate_ms1_scan(scan)
        self._iterator = self._iterate_over_index(scan.index)
        return self

    # LC-MS/MS Database API

    def msms_for(self, neutral_mass, mass_error_tolerance=1e-5, start_time=None, end_time=None):
        m = neutral_mass
        w = neutral_mass * mass_error_tolerance
        q = self.session.query(PrecursorInformation).join(
            PrecursorInformation.precursor).filter(
            PrecursorInformation.neutral_mass.between(m - w, m + w),
            PrecursorInformation.sample_run_id == self.sample_run_id).order_by(
            MSScan.scan_time)
        if start_time is not None:
            q = q.filter(MSScan.scan_time >= start_time)
        if end_time is not None:
            q = q.filter(MSScan.scan_time <= end_time)
        return q

    def ms1_peaks_above(self, threshold=1000):
        accumulate = [
            (x[0], x[1].convert(), x[1].id) for x in self.session.query(MSScan.scan_id, DeconvolutedPeak).join(
                DeconvolutedPeak).filter(
                MSScan.ms_level == 1, MSScan.sample_run_id == self.sample_run_id,
                DeconvolutedPeak.neutral_mass > threshold
            ).order_by(MSScan.index).yield_per(1000)]
        return accumulate

    def precursor_information(self):
        prec_info = self.session.query(PrecursorInformation).filter(
            PrecursorInformation.sample_run_id == self.sample_run_id).all()
        return prec_info
