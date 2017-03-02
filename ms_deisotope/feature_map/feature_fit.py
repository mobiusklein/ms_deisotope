from ms_deisotope.peak_dependency_network.intervals import SpanningMixin


class LCMSFeatureSetFit(object):
    def __init__(self, features, theoretical, score, charge,
                 missing_features=0, data=None):
        self.features = features
        self.theoretical = theoretical
        self.score = score
        self.charge = charge
        self.data = data
        self.missing_features = missing_features
        self.monoisotopic_feature = features[0]

    def clone(self):
        return self.__class__(
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.data)

    def __reduce__(self):
        return self.__class__, (
            self.features, self.theoretical, self.score, self.charge,
            self.missing_features, self.data)

    def __eq__(self, other):
        val = (self.score == other.score and
               self.charge == other.charge and
               self.features == other.features and
               self.theoretical == other.theoretical)
        if self.data is not None or other.data is not None:
            val = val and (self.data == other.data)
        return val

    def __ne__(self, other):
        return not (self == other)

    def __lt__(self, other):
        return self.score < other.score

    def __gt__(self, other):
        return self.score > other.score

    def __hash__(self):
        return hash((self.monoisotopic_feature.mz, self.charge))

    def __iter__(self):
        yield self.score
        yield self.charge
        yield self.features
        yield self.theoretical

    def __len__(self):
        return len(self.features)

    @property
    def npeaks(self):
        return len(self)

    def __repr__(self):
        return "LCMSFeatureSetFit(score=%0.5f, charge=%d, size=%d, monoisotopic_mz=%0.5f)" % (
            self.score, self.charge, len(self), self.monoisotopic_feature.mz)


class FeatureNode(object):
    def __init__(self, feature, links=None):
        if links is None:
            links = {}
        self.feature = feature
        self.links = links
        self._hash = hash(self.feature)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        try:
            return self.feature == other.feature
        except AttributeError:
            return False

    def __ne__(self, other):
        return not self == other

    def __contains__(self, fit):
        return fit in self.links

    def __repr__(self):
        return "FeatureNode(%s, %s)" % (self.feature, self.links)


class FeatureSetFitNode(SpanningMixin):
    def __init__(self, fit, index=None):
        self.fit = fit
        self.edges = set()
        self.overlap_edges = set()

        self._hash = None
        self.score = fit.score

        self.feature_indices = {(f.mz, f.start_time, f.end_time) for f in self.features if f is not None}

        self.start = (fit.features[0].mz, max(f.start_time for f in fit.features if f is not None))
        self.end = (fit.features[-1].mz, min(f.end_time for f in fit.features if f is not None))

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(self.fit)
        return self._hash

    def __eq__(self, other):
        return self.fit is other.fit

    def __gt__(self, other):
        return self.fit > other.fit

    def __lt__(self, other):
        return self.fit < other.fit

    def __ne__(self, other):
        return self.fit is not other.fit

    def visit(self, other):
        if self.isdisjoint(other):
            self.edges.add(other)
            other.edges.add(self)
        else:
            self.overlap_edges.add(other)
            other.overlap_edges.add(self)

    def isdisjoint(self, other):
        return self.feature_indices.isdisjoint(other.feature_indices)

    def __repr__(self):
        return "FeatureSetFitNode(%r)" % self.fit
