Formats
~~~~~~~

mzML and mzXML
==============

:mod:`ms_deisotope` supports reading from :mod:`mzML <ms_deisotope.data_source.mzml>` [Martens2011]_
and :mod:`mzXML <ms_deisotope.data_source.mzxml>` [Pedrioli2004a]_ files on all platforms, and can
provide fast random access to uncompressed files. By default, these file types are assumed to be
collections of MS1 and MSn scans, and will iterate over scan bunches. To iterate over single scans,
call :meth:`make_iterator` or another iterator creation method with the keyword argument ``grouped=False``.


MGF
===

:mod:`ms_deisotope` supports reading :mod:`MGF <ms_deisotope.data_source.mgf>` files on all platforms,
and can provide fast random access to uncompressed files. As this format does not store MS1 scans, only
single MSn scans are produced by iteration, as if ``grouped=False`` were passed to :meth:`make_iterator`.


Vendor Readers
==============

If the appropriate COM DLL has been registered with Windows, :mod:`ms_deisotope` is able to read a subset
of the commonly used vendor data file formats. These depend upon the :mod:`comtypes` library, which can be
installed from the Python Package Index.

Thermo Fisher RAW
*****************

If the `MSFileReader <https://thermo.flexnetoperations.com/control/thmo/index>`_ package that provides the
``XRawfile2_<arch>.dll`` library has been installed and registered with Windows, :mod:`ms_deisotope` can open
:mod:`Thermo Fisher RAW <ms_deisotope.data_source.thermo_raw>` files. The implementation of the MSFileReader
bindings are derived from the work of `François Allain <https://github.com/frallain/MSFileReader-Python-bindings>`_,
distributed under the MIT license, which have been included in this codebase.


Agilent .d
**********

If Agilent's MassSpecDataReader.dll and its supporting libraries have been installed, and had their
type libraries built and registered (files whose names match the DLL name, but whose extension is ``.tlb``)
as described in the vendor's installation instructions, :mod:`ms_deisotope` can open
:mod:`Agilent .d <ms_deisotope.data_source.agilent_d>` directories directly. Because there is no standard registered
installation location, you must explicitly tell :mod:`comtypes` where to look for the DLLs before you
can use this feature. See :func:`ms_deisotope.data_source.agilent_d.register_dll_dir` for more details.



References
~~~~~~~~~~

.. [Martens2011]
    Martens, L., Chambers, M., Sturm, M., Kessner, D., Levander, F., Shofstahl, J.,
    … Deutsch, E. W. (2011). mzML--a community standard for mass spectrometry data.
    Molecular & Cellular Proteomics : MCP, 10(1), R110.000133.
    https://doi.org/10.1074/mcp.R110.000133
.. [Pedrioli2004a]
    Pedrioli, P. G. A., Eng, J. K., Hubley, R., Vogelzang, M., Deutsch, E. W., Raught, B.,
    … Aebersold, R. (2004). A common open representation of mass spectrometry data and its
    application to proteomics research. Nature Biotechnology, 22(11), 1459–1466.
    https://doi.org/10.1038/nbt1031