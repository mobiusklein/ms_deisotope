A Library for Deisotoping and Charge State Deconvolution For Mass Spectrometry
------------------------------------------------------------------------------

This library combines `brainpy` and `ms_peak_picker` to build a toolkit for
MS and MS/MS data. The goal of these libraries is to provide pieces of the puzzle
for evaluating MS data modularly. The goal of this library is to combine the modules
to streamline processing raw data.


API
---


Averagine
=========

An "Averagine" model is used to describe the composition of an "average amino acid",
which can then be used to approximate the composition and isotopic abundance of a
combination of specific amino acids. Given that often the only solution available is
to guess at the composition of a particular `m/z` because there are too many possible
elemental compositions, this is the only tractable solution.

This library supports arbitrary Averagine formulae, but the Senko Averagine is provided
by default: `{"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}`

.. code:: python

    from ms_deisotope import Averagine
    from ms_deisotope import utils

    peptide_averagine = Averagine({"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417})
    
    utils.draw_peaklist(peptide_averagine.isotopic_cluster(1266.321, charge=1))


