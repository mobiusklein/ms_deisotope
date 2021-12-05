Formats
~~~~~~~

There are many different types of mass spectrometry data files, each requiring different methods
to read them. When working with a path or file-like object, you can use
:func:`~ms_deisotope.data_source.infer_format.MSFileLoader` to open most formats using the same
uniform API. Otherwise, you may directly use the file reading classes in the appropriate data source
submodule.

mzML and mzXML
==============

:mod:`ms_deisotope` supports reading from :mod:`mzML <ms_deisotope.data_source.mzml>` [Martens2011]_
and :mod:`mzXML <ms_deisotope.data_source.mzxml>` [Pedrioli2004a]_ files on all platforms, and can
provide fast random access to uncompressed files. By default, these file types are assumed to be
collections of MS1 and MSn scans, and will iterate over scan bunches. To iterate over single scans,
call :meth:`~.ScanIterator.make_iterator` or another iterator creation method with the keyword argument ``grouped=False``.

mzMLb
=====
With :mod:`h5py`, :mod:`ms_deisotope` supports reading :mod:`mzMLb <ms_deisotope.data_source.mzmlb>` [Bhamber2021]_.
The reader provides identical behaviors to those for mzML, while providing full random access even while
compressed conventionally.

MGF
===

:mod:`ms_deisotope` supports reading :mod:`MGF <ms_deisotope.data_source.mgf>` files on all platforms,
and can provide fast random access to uncompressed files. As this format does not store MS1 scans, only
single MSn scans are produced by iteration, as if ``grouped=False`` were passed to :meth:`~.ScanIterator.make_iterator`.


Vendor Readers
==============

If the appropriate COM DLL has been registered with Windows, :mod:`ms_deisotope` is able to read a subset
of the commonly used vendor data file formats. These depend upon the platform and external libraries. Some
of these external libraries may be automatically detected while others will need to be formally registered
before use.

.. note::
    By using an instrument vendor's library to read their proprietary file format, you are agreeing to
    the requisite license and terms associated with that library.


Thermo Fisher RAW
*****************

If the .NET run time is available and :mod:`pythonnet` is installed, Thermo's ``RawFileReader`` library will be used
to open :mod:`Thermo Fisher RAW <ms_deisotope.data_source.thermo_raw_net>`.

With :mod:`comtypes`, if the `MSFileReader <https://thermo.flexnetoperations.com/control/thmo/index>`_ package
that provides the ``XRawfile2_<arch>.dll`` library has been installed and registered with Windows, :mod:`ms_deisotope`
can open :mod:`Thermo Fisher RAW <ms_deisotope.data_source.thermo_raw>` files. The implementation of the MSFileReader
bindings are derived from the work of `François Allain <https://github.com/frallain/MSFileReader-Python-bindings>`_,
distributed under the MIT license, which have been included in this codebase.


Agilent .d
**********

With :mod:`comtypes`, if Agilent's MassSpecDataReader.dll and its supporting libraries have been installed, and
had their type libraries built and registered (files whose names match the DLL name, but whose extension is ``.tlb``)
as described in the vendor's installation instructions, :mod:`ms_deisotope` can open
:mod:`Agilent .d <ms_deisotope.data_source.agilent_d>` directories directly. Because there is no standard registered
installation location, you must explicitly tell :mod:`comtypes` where to look for the DLLs before you
can use this feature. See :func:`ms_deisotope.data_source.agilent_d.register_dll_dir` for more details.


Waters
******

On Windows, if the MassLynx SDK is installed either on the system path or are registered with :mod:`ms_deisotope`'s
configuration file, :mod:`ms_deisotope` can open :mod:`Waters .RAW <ms_deisotope.data_source.masslynx>` directories directly.
To programmatically register the SDK, see :func:`ms_deisotope.data_source._vendor.masslynx.libload.register_dll` for more
details.


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
.. [Bhamber2021] Bhamber, R. S., Jankevics, A., Deutsch, E. W., Jones, A. R., & Dowsey, A. W. (2021).
    MzMLb: A Future-Proof Raw Mass Spectrometry Data Format Based on Standards-Compliant
    mzML and Optimized for Speed and Storage Requirements. Journal of Proteome Research,
    20(1), 172–183. https://doi.org/10.1021/acs.jproteome.0c00192