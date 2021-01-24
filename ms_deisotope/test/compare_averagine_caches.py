'''A helper to verify the determinism of AveragineCache objects used across multiple processes
behaved identically.

Usage:
    python compare_averagine_caches.py <cache_files*.pkl>
'''
import glob
import sys
import pprint
try:
    import cPickle as pickle
except ImportError:
    import pickle


def load_caches_from_paths(paths):
    caches = []
    for path in paths:
        for p in glob.glob(path):
            with open(p, 'rb') as fh:
                caches.append(pickle.load(fh))
    return caches


def compare_caches(caches):
    print("Comparing Key Sets")
    for cache in caches:
        for cache2 in caches:
            if cache is cache2:
                continue
            if set(cache.backend.keys()) != set(cache2.backend.keys()):
                pprint.pprint(set(cache.backend.keys()) - set(cache2.backend.keys()))
                raise AssertionError("Key Sets don't match")
    print("Comparing Value Sets")
    for cache in caches:
        for cache2 in caches:
            if cache is cache2:
                continue
            for key, value in cache.backend.items():
                assert value == cache2.backend[key]


if __name__ == "__main__":
    args = sys.argv[1:]
    caches = load_caches_from_paths(args)
    compare_caches(caches)
