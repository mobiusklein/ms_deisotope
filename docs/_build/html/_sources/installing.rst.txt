

Installing ms_deisotope
-----------------------

:mod:`ms_deisotope` was written to be an fusion between a simple peak picking
library, `ms_peak_picker <https://github.com/mobiusklein/ms_peak_picker>`_ and
a simple isotopic pattern generation library `brainpy <https://github.com/mobiusklein/brainpy>`_.

It later grew a data access layer written directly on top of `pyteomics <https://pyteomics.readthedocs.io/en/latest/>`_
and `lxml <https://lxml.de/>`_.

It has two optional dependencies: `Matplotlib <https://matplotlib.org/>`_ for drawing spectra, and
`psims <https://github.com/mobiusklein/psims>`_ for writing mzML and updating its compiled controlled
vocabullary.


Building From Source
====================

:mod:`ms_deisotope` is written for direct use from Python, but large
chunks of its internals have optional but highly recommended C implementations
that significantly improve the speed of the provided algorithms. To build the 
library from source with the C extensions included, you'll need:

    1. A C compiler matching your Python version
    2. A recent version of `Cython <http://cython.org/>`_ >= 0.27.0
    3. NumPy and SciPy
    4. `ms_peak_picker <https://github.com/mobiusklein/ms_peak_picker>`_ built with C extensions
    5. `brainpy <https://github.com/mobiusklein/brainpy>`_ built with C extensions

.. note::
    :mod:`ms_deisotope`'s Python implementations and C implementations should
    behave in the same way. If you think they are not, please
    `open a bug report <https://github.com/mobiusklein/ms_deisotope/issues/new>`_ and
    let us know!

Once you have all of these dependencies installed, running ``python setup.py install`` should take
care of the rest. If you want to profile the C components, the ``--include-diagonistics`` will enable
profiling hooks.


