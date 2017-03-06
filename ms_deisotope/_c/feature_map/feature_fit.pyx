cimport cython

from ms_deisotope._c.feature_map.lcms_feature cimport LCMSFeature


@cython.freelist(1000000)
cdef class map_coord(object):
    def __init__(self, mz, time):
        self.mz = mz
        self.time = time
    
    def __getitem__(self, i):
        if i == 0:
            return self.mz
        elif i == 1:
            return self.time
        else:
            raise IndexError(i)

    def __iter__(self):
        yield self.mz
        yield self.time

    def __reduce__(self):
        return map_coord, (self.mz, self.time,)


cdef class LCMSFeatureSetFit(object):
    def __init__(self, features, theoretical, score, charge,
                 missing_features=0, data=None):
        self.features = features
        self.theoretical = theoretical
        self.score = score
        self.charge = charge
        self.data = data
        self.missing_features = missing_features
        self.monoisotopic_feature = features[0]

    @property
    def mz(self):
        return self.monoisotopic_feature.mz

    def clone(self):
        return self.__class__(
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.data)

    def __reduce__(self):
        return self.__class__, (
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.data)

    cpdef bint _eq(self, LCMSFeatureSetFit other):
        cdef bint val
        val = (self.score == other.score and
               self.charge == other.charge and
               self.features == other.features and
               self.theoretical == other.theoretical)
        if self.data is not None or other.data is not None:
            val = val & (self.data == other.data)
        return val

    cpdef bint _ne(self, LCMSFeatureSetFit other):
        return not (self == other)

    cpdef bint _lt(self, LCMSFeatureSetFit other):
        return self.score < other.score

    cpdef bint _gt(self, LCMSFeatureSetFit other):
        return self.score > other.score

    def __richcmp__(self, LCMSFeatureSetFit other, int code):
        if other is None:
            if code == 3:
                return True
            else:
                return False

        if code == 0:
            return self._lt(other)
        elif code == 2:
            return self._eq(other)
        elif code == 3:
            return not (self._eq(other))
        elif code == 4:
            return self._gt(other)

    def __hash__(self):
        return hash((self.monoisotopic_feature.mz, self.charge))

    def __iter__(self):
        return iter(self.features)

    def __len__(self):
        return len(self.features)

    @property
    def npeaks(self):
        return len(self)

    def __repr__(self):
        return "LCMSFeatureSetFit(score=%0.5f, charge=%d, size=%d, monoisotopic_mz=%0.5f)" % (
            self.score, self.charge, len(self), self.monoisotopic_feature.mz)

    @property
    def start(self):
        first = self.features[0]
        if first is None:
            raise Exception()
        return map_coord(first.mz, first.start_time)

    @property
    def end(self):
        for last in reversed(self.features):
            if last is None:
                continue
            return map_coord(last.mz, last.end_time)
