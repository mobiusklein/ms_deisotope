import ms_deisotope

from ms_deisotope.feature_map import quick_index
from ms_deisotope.test.common import datafile


mzml_path = datafile("small.mzML")


def test_quick_index():
    reader = ms_deisotope.MSFileLoader(mzml_path)
    index, _interval_tree = quick_index.index(reader)
    n_1 = len(index.ms1_ids)
    n_n = len(index.msn_ids)
    assert n_1 == 14
    assert n_n == 34
