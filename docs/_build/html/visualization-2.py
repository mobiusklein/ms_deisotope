import ms_deisotope
from ms_deisotope import plot
from ms_deisotope.test.common import datafile

reader = ms_deisotope.MSFileLoader(datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
bunch = next(reader)

bunch.precursor.pick_peaks()
bunch.precursor.deconvolute(
    scorer=ms_deisotope.PenalizedMSDeconVFitter(20., 2.0),
    averagine=ms_deisotope.glycopeptide, use_quick_charge=True)

ax = plot.draw_peaklist(bunch.precursor, color='black')
ax = plot.annotate_isotopic_peaks(bunch.precursor, ax=ax)
ax.set_xlim(1160, 1165)
ax.figure.set_figwidth(12)