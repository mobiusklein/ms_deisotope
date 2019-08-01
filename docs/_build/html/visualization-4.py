import ms_deisotope
from ms_deisotope import plot
from ms_deisotope.test.common import datafile

reader = ms_deisotope.MSFileLoader(datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
bunch = next(reader)

bunch.precursor.pick_peaks()
bunch.precursor.deconvolute(
    scorer=ms_deisotope.PenalizedMSDeconVFitter(20., 2.0),
    averagine=ms_deisotope.glycopeptide, use_quick_charge=True)

ax = plot.annotate_scan(bunch.precursor, bunch.products, nperrow=2)
ax.figure.set_figwidth(12)
ax.figure.set_figheight(16)