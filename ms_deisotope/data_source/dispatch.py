import os

try:
    from collections import Sequence
except ImportError:
    from collections.abc import Sequence


from .infer_type import MSFileLoader
from .scan import ScanProxyContext


class SubsequenceMappingProxy(object):
    '''Poor man's radix tree.
    '''

    def __init__(self, mapping):
        self.mapping = mapping

    def _resolve(self, key):
        if key in self.mapping:
            return self.mapping[key]
        for k, v in self.mapping.items():
            if not isinstance(k, Sequence):
                continue
            if key in k:
                return v
        best_substr = None
        best_substr_len = -1
        if not isinstance(key, Sequence):
            raise KeyError(key)
        for k, v in self.mapping.items():
            if k in key:
                substr_len = len(k)
                if substr_len > best_substr_len:
                    best_substr = k
                    best_substr_len = substr_len
        if best_substr_len > 0:
            return self.mapping[best_substr]
        raise KeyError(key)

    def __getitem__(self, key):
        return self._resolve(key)

    def __contains__(self, key):
        try:
            found = self[key]
            return found is not None
        except KeyError:
            return False

    def __iter__(self):
        return iter(self.mapping)


class DynamicallyLoadingResolver(object):
    '''Dynamically open data files as they are encountered
    '''

    def __init__(self, root_path=".", mapping=None, opener=None):
        if mapping is None:
            mapping = dict()
        if opener is None:
            opener = self._default_opener
        self.root_path = root_path
        self.mapping = mapping
        self.opener = opener

    def _default_opener(self, path):
        return MSFileLoader(path)

    def _resolve(self, path):
        if path in self.mapping:
            return self.mapping[path]
        _path = os.path.join(self.root_path, path)
        if not os.path.exists(_path):
            raise KeyError(path)
        opened = self.opener(_path)
        self.mapping[path] = opened
        return opened

    def __getitem__(self, key):
        return self._resolve(key)

    def __contains__(self, key):
        try:
            found = self[key]
            return found is not None
        except KeyError:
            return False


class DynamicallyLoadingProxyResolver(DynamicallyLoadingResolver):

    def _default_opener(self, path):
        reader = super(DynamicallyLoadingProxyResolver, self)._default_opener(path)
        return ScanProxyContext(reader)
