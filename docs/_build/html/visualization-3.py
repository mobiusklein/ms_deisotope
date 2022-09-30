import ms_deisotope
from ms_deisotope import plot
from ms_deisotope.test.common import datafile

example_file = datafile("20150710_3um_AGP_001_29_30.mzML.gz")

reader = ms_deisotope.MSFileLoader(example_file)
bunch = next(reader)

bunch.precursor.pick_peaks()
bunch.precursor.deconvolute(
    scorer=ms_deisotope.PenalizedMSDeconVFitter(20., 2.0),
    averagine=ms_deisotope.glycopeptide, use_quick_charge=True)

ax = plot.annotate_scan_single(bunch.precursor, bunch.products[0])
ax.figure.set_figwidth(12)