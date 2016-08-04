from uuid import uuid4

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.orm.session import object_session
from sqlalchemy.engine import Connectable

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import (
    Column, Numeric, Integer, String, ForeignKey, PickleType,
    Boolean)
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.mutable import Mutable

from sqlalchemy import exc

import numpy as np

from ms_peak_picker import FittedPeak as MemoryFittedPeak, PeakIndex, PeakSet
from ms_deisotope import DeconvolutedPeak as MemoryDeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.peak_set import Envelope

from ms_deisotope.data_source.common import (
    ProcessedScan, PrecursorInformation as MemoryPrecursorInformation, ScanBunch)

from .common import ScanSerializerBase, ScanDeserializerBase


def Mass():
    return Column(Numeric(14, 6, asdecimal=False), index=True)


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


class SampleRun(Base):
    __tablename__ = "SampleRun"

    id = Column(Integer, primary_key=True, autoincrement=True)
    uuid = Column(String(32), index=True, unique=True)
    name = Column(String(128), unique=True)

    ms_scans = relationship(
        "MSScan", backref=backref("sample_run"), lazy='dynamic')
    sample_type = Column(String(128))

    completed = Column(Boolean(), default=False, nullable=False)

    def save_bunch(self, bunch):
        session = object_session(self)
        return serialize_scan_bunch(session, bunch, self.id)


class MSScan(Base):
    __tablename__ = "MSScan"

    id = Column(Integer, primary_key=True, autoincrement=True)
    index = Column(Integer, index=True)
    ms_level = Column(Integer)
    scan_time = Column(Numeric(10, 5, asdecimal=False), index=True)
    title = Column(String(128), index=True)
    scan_id = Column(String(128), index=True)
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id, ondelete='CASCADE'), index=True)

    peak_set = relationship("FittedPeak", backref="scan", lazy="dynamic")
    deconvoluted_peak_set = relationship(
        "DeconvolutedPeak", backref='scan', lazy='dynamic')

    def __repr__(self):
        f = "{}({}, {}, {}, {}".format(
            self.__class__.__name__, self.scan_id, self.ms_level, self.scan_time,
            self.deconvoluted_peak_set.count())
        if self.ms_level > 1:
            f = "%s %s" % (f, self.precursor_information)
        f += ")"
        return f

    def convert(self):
        precursor_information = self.precursor_information.convert() if self.precursor_information is not None else None

        peak_set_items = [p.convert() for p in self.peak_set]
        peak_set = PeakSet(peak_set_items)
        peak_set._index()
        peak_index = PeakIndex(np.array([], dtype=np.float64), np.array([], dtype=np.float64), peak_set)

        deconvoluted_peak_set_items = [p.convert() for p in self.deconvoluted_peak_set]
        deconvoluted_peak_set = DeconvolutedPeakSet(deconvoluted_peak_set_items)
        deconvoluted_peak_set._reindex()

        scan = ProcessedScan(
            self.scan_id, self.title, precursor_information, self.ms_level,
            self.scan_time, self.index, peak_index, deconvoluted_peak_set)
        return scan

    @classmethod
    def serialize(self, scan, sample_run_id=None):
        db_scan = MSScan(
            index=scan.index, ms_level=scan.ms_level,
            scan_time=scan.scan_time, title=scan.title,
            scan_id=scan.id, sample_run_id=sample_run_id)
        db_scan.peak_set = map(FittedPeak.serialize, scan.peak_set)
        db_scan.deconvoluted_peak_set = map(DeconvolutedPeak.serialize, scan.deconvoluted_peak_set)
        return db_scan

    @classmethod
    def serialize_bulk(self, scan, sample_run_id, session):
        db_scan = MSScan(
            index=scan.index, ms_level=scan.ms_level,
            scan_time=scan.scan_time, title=scan.title,
            scan_id=scan.id, sample_run_id=sample_run_id)

        session.add(db_scan)
        session.flush()

        FittedPeak._serialize_bulk_list(scan.peak_set, db_scan.id, session)
        DeconvolutedPeak._serialize_bulk_list(scan.deconvoluted_peak_set, db_scan.id, session)
        return db_scan


class PrecursorInformation(Base):
    __tablename__ = "PrecursorInformation"

    id = Column(Integer, primary_key=True)
    sample_run_id = Column(Integer, ForeignKey(SampleRun.id, ondelete='CASCADE'), index=True)

    precursor_id = Column(Integer, ForeignKey(MSScan.id), index=True)
    precursor = relationship(MSScan, backref=backref("product_information"),
                             primaryjoin="PrecursorInformation.precursor_id == MSScan.id",
                             uselist=False)
    product_id = Column(Integer, ForeignKey(MSScan.id), index=True)
    product = relationship(MSScan, backref=backref("precursor_information", uselist=False),
                           primaryjoin="PrecursorInformation.product_id == MSScan.id",
                           uselist=False)
    neutral_mass = Mass()
    charge = Column(Integer)
    intensity = Column(Numeric(12, 4, asdecimal=False))

    def __repr__(self):
        return "PrecursorInformation({}, {}, {})".format(
            self.precursor.scan_id, self.neutral_mass, self.charge)

    def convert(self, data_source=None):
        return MemoryPrecursorInformation(
            0, 0, 0, self.precursor.scan_id, data_source, self.neutral_mass, self.charge, self.intensity)


class PeakMixin(object):
    mz = Column(Numeric(12, 6, asdecimal=False), index=True)
    intensity = Column(Numeric(12, 4, asdecimal=False))
    full_width_at_half_max = Column(Numeric(8, 7, asdecimal=False))
    signal_to_noise = Column(Numeric(8, 7, asdecimal=False))
    scan_peak_index = Column(Integer, index=True)
    area = Column(Numeric(12, 6, asdecimal=False))

    @classmethod
    def _serialize_bulk_list(cls, peaks, scan_id, session):
        out = []
        for peak in peaks:
            db_peak = cls.serialize(peak)
            db_peak.scan_id = scan_id
            out.append(db_peak)
        session.bulk_save_objects(out)


class FittedPeak(Base, PeakMixin):
    __tablename__ = "FittedPeak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    scan_id = Column(Integer, ForeignKey(MSScan.id, ondelete='CASCADE'), index=True)

    def convert(self):
        return MemoryFittedPeak(
            self.mz, self.intensity, self.signal_to_noise, -1, self.scan_peak_index,
            self.full_width_at_half_max, self.area)

    @classmethod
    def serialize(cls, peak):
        return cls(
            mz=peak.mz, intensity=peak.intensity, signal_to_noise=peak.signal_to_noise,
            scan_peak_index=peak.index, full_width_at_half_max=peak.full_width_at_half_max,
            area=peak.area)


class DeconvolutedPeak(Base, PeakMixin):
    __tablename__ = "DeconvolutedPeak"

    id = Column(Integer, primary_key=True, autoincrement=True)
    neutral_mass = Mass()
    average_mass = Mass()
    most_abundant_mass = Mass()

    charge = Column(Integer)

    score = Column(Numeric(12, 6, asdecimal=False))
    scan_id = Column(Integer, ForeignKey(MSScan.id, ondelete='CASCADE'), index=True)
    envelope = Column(MutableList.as_mutable(PickleType))
    a_to_a2_ratio = Column(Numeric(8, 7, asdecimal=False))

    def convert(self):
        return MemoryDeconvolutedPeak(
            self.neutral_mass, self.intensity, self.charge,
            self.signal_to_noise, -1, self.full_width_at_half_max,
            self.a_to_a2_ratio, self.most_abundant_mass, self.average_mass,
            self.score, Envelope(self.envelope), self.mz)

    @classmethod
    def serialize(cls, peak):
        return cls(
            mz=peak.mz, intensity=peak.intensity, signal_to_noise=peak.signal_to_noise,
            scan_peak_index=peak.index.neutral_mass, full_width_at_half_max=peak.full_width_at_half_max,
            neutral_mass=peak.neutral_mass, average_mass=peak.average_mass,
            most_abundant_mass=peak.most_abundant_mass, charge=peak.charge, score=peak.score,
            envelope=list(peak.envelope), a_to_a2_ratio=peak.a_to_a2_ratio)


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
        db_pi = PrecursorInformation(
            precursor_id=db_precursor.id, product_id=db_scan.id,
            charge=pi.extracted_charge, intensity=pi.extracted_intensity,
            neutral_mass=pi.extracted_neutral_mass, sample_run_id=sample_run_id)
        session.add(db_pi)
    session.flush()
    return db_precursor, db_products


def serialize_scan_bunch_bulk(session, bunch, sample_run_id):
    precursor = bunch.precursor
    db_precursor = MSScan.serialize_bulk(precursor, sample_run_id, session)
    db_products = [MSScan.serialize_bulk(p, sample_run_id, session)
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


def initialize(path):
    eng = create_engine("sqlite:///%s" % path)
    Base.metadata.create_all(bind=eng)
    return eng


class DatabaseScanSerializer(ScanSerializerBase):
    def __init__(self, connection, sample_name=None, overwrite=True):
        self.uuid = str(uuid4())

        if sample_name is None:
            sample_name = self.uuid

        self.sample_name = sample_name
        self.engine = self._configure_connection(connection)

        self._sessionmaker = scoped_session(sessionmaker(bind=self.engine, autoflush=False))
        self._session = None
        self._sample_run = None
        self._sample_run_id = None

        self._overwrite(overwrite)

    def _overwrite(self, overwrite=True):
        try:
            sample_run = self.sample_run
        except exc.IntegrityError:
            self.session.rollback()
            sample_run = self.session.query(SampleRun).filter(SampleRun.name == self.sample_name).first()
            if not sample_run.completed or overwrite:
                self.session.delete(sample_run)
                self.session.commit()
                self._construct_sample_run()

    def _construct_sample_run(self):
        sr = SampleRun(uuid=self.uuid, name=self.sample_name)
        self.session.add(sr)
        self.session.commit()
        self._sample_run = sr
        self._sample_run_id = sr.id

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
    def session(self):
        return self._sessionmaker()

    def _configure_connection(self, connection):
        if isinstance(connection, basestring):
            try:
                eng = create_engine(connection)
            except exc.ArgumentError:
                eng = create_engine("sqlite:///%s" % connection)
        elif isinstance(connection, Connectable):
            eng = connection
        Base.metadata.create_all(bind=eng)
        return eng

    def save(self, bunch, commit=True, bulk=True):
        if bulk:
            out = serialize_scan_bunch_bulk(self.session, bunch, self.sample_run_id)
        else:
            out = serialize_scan_bunch(self.session, bunch, self.sample_run_id)
        if commit:
            self.session.commit()
        return out

    def query(self, *args):
        return self.session.query(*args)


def flatten(iterable):
    return [y for x in iterable for y in x]


class DatabaseScanDeserializer(ScanDeserializerBase):

    def __init__(self, connection, sample_name=None):

        self.sample_name = sample_name

        self.engine = self._configure_connection(connection)
        self._sessionmaker = scoped_session(sessionmaker(bind=self.engine, autoflush=False))
        self._session = None
        self._sample_run = None
        self._sample_run_id = None
        self._iterator = None

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
    def session(self):
        return self._sessionmaker()

    def _retrieve_sample_run(self):
        session = self.session
        sr = session.query(SampleRun).filter(SampleRun.name == self.sample_name).one()
        self._sample_run = sr
        self._sample_run_id = sr.id

    def _configure_connection(self, connection):
        if isinstance(connection, basestring):
            try:
                eng = create_engine(connection)
            except exc.ArgumentError:
                eng = create_engine("sqlite:///%s" % connection)
        elif isinstance(connection, Connectable):
            eng = connection
        return eng

    def get_scan_by_id(self, scan_id):
        q = self._get_by_scan_id(scan_id)
        return q.convert()

    def _iterate_over_index(self, start=0, require_ms1=True):
        indices_q = self.session.query(MSScan.index).filter(
            MSScan.sample_run_id == self.sample_run_id)
        if require_ms1:
            indices_q = indices_q.filter(MSScan.ms_level == 1)
        indices = flatten(indices_q.all())
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

        while i < len(items):
            index = items[i]
            scan = self.session.query(MSScan).filter(
                MSScan.index == index,
                MSScan.sample_run_id == self.sample_run_id).one()
            products = [pi.product for pi in scan.product_information]
            yield ScanBunch(scan.convert(), [p.convert() for p in products])
            i += 1

    def __iter__(self):
        return self

    def __next__(self):
        if self._iterator is None:
            self._iterator = self._iterate_over_index()
        return next(self._iterator)

    def _get_by_scan_id(self, scan_id):
        q = self.session.query(MSScan).filter(
            MSScan.scan_id == scan_id, MSScan.sample_run_id == self.sample_run_id).one()
        return q

    def _get_scan_by_time(self, rt, require_ms1=False):
        times_q = self.session.query(MSScan.scan_time).filter(
            MSScan.sample_run_id == self.sample_run_id)
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

    def get_scan_by_time(self, rt, require_ms1=False):
        q = self._get_scan_by_time(rt, require_ms1)
        return q.convert()

    def _get_scan_by_index(self, index):
        q = self.session.query(MSScan).filter(
            MSScan.index == index, MSScan.sample_run_id == self.sample_run_id).one()
        return q

    def get_scan_by_index(self, index):
        return self._get_scan_by_index(index).convert()

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
