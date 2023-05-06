# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][Keep a Changelog] and this project adheres to [Semantic Versioning][Semantic Versioning].

## [Unreleased]

### Added
1. More type annotations

### Changed
1. Arrays used to reconstruct deconvoluted peak sets are not guaranteed to be copied.

### Removed

### Fixed
1. `WrappedScan` objects now properly propagate overridden attributes when copied.
2. More docstring formatting was fixed and missing arguments were described

### Security


## [v0.0.46] - 2023-01-22


### Added
1. Documentation for the `ms_deisotope.peak_set` module.
2. Added `FunctionScorer` to `ms_deisotope.scoring` to allow user-defined scoring functions to be used without
   sub-classing. This incurs some extra overhead on every call, so performance may be worse than a built-in
   sub-class.

### Changed
1. `ms_deisotope.deconvolution.prepare_peaklist` now coerces a tuple of parallel sequences or
   `numpy` arrays of m/z and intensities properly, so the caller does not need to zip them into
   a list of `(mz, intensity)` pairs.

### Fixed
1. `ms_deisotope.deconvolute_peaks` was incorrectly calling `ms_deisotope.deconvolution.prepare_peaklist` too late, so no useful
   coercion was done for the user. `prepare_peaklist` is now correctly called **before** passing the
   provided peak list to the deconvolution machinery, removing the need for users to call it explicitly.


## [v0.0.45] - 2023-01-14

### Added
- More type annotations for `ProcessedScan`-yielding readers.
- C accelerators for `pyteomics`'s unit-scalars and CV-params to speed up parsing and reduce memory usage.
- Created a low level C API for stripped-down peak lists that does not use the GIL, limited use-cases currently.
- Loads of type annotations.
- Added a `-X` option to pass raw key-value pairs to the mzML converter in `ms-index mzml` to make experimental
  configuration easier.

## [v0.0.44] - 2022-10-09

### Added

### Changed
1. The minimum NumPy version on Python 3.10 was increased to 1.23.2

### Deprecated

### Removed

### Fixed
1. Fix software ID in mzML writing.

### Security


## [v0.0.43] - 2022-09-29

### Added
1. Add support for `os.PathLike` objects to `ms_deisotope.MSFileLoader`.

### Changed
1. mzML and mzXML now properly capture the relationship between instrument configuration and software.

### Deprecated

### Removed

### Fixed
1. Fixed division in `_binsearch` used for `precursor-clustering`.

### Security


## [v0.0.42] - 22-06-15

### Added
1. Added `close` parameter to `MzMLSerializer` and `MzMLbSerializer` to reflect the change made upstream in
   `psims`. These govern the automatic closing behavior of the buffer after completing the data document. If
   you *need* the buffer to remain open after completing file writing, pass `close=False` to the writer.

### Changed

### Deprecated

### Removed

### Fixed

### Security


## [v0.0.40] - 22-06-14

### Added
- Make `CVPrimitiveWrapperMeta` hashable.

### Changed
- Updated automatic release workflow
- Allow `ProcessingMethod` to repeat terms

### Deprecated

### Removed

### Fixed

### Security


## [v0.0.39] - 22-05-29

### Added

### Changed
- When using `ms_deisotope.test.common.datafile`, if the file isn't found in the
  package data directory, it will be downloaded from GitHub. Downloaded files will
  be stored in the `test_data/` directory under the path returned by
  `ms_deisotope.config.get_config_dir`.
- Use internal C extensions to speed up unit primitives.
- When `isal` is installed, use `isal.isal_zlib` in place of the stock `zlib`.


### Deprecated

### Removed

### Fixed

### Security



## [v0.0.38] - 22-04-28

### Added
- Added `IonMobilityAware3DMzMLbSerializer` to `ms_deisotope.output.mzmlb`.
- Added IBT10 reporter ions to `TMTReporterExtractor`
- Added `extract-reporter-ions` to `ms-index` to extract reporter ions to CSV.

### Changed
- Raised minimum `ms_peak_picker` version to `v0.1.38`
- Changed the default `dx` of `Scan.average` and `Scan.average_with` from 0.01 to 0.002 to reduce distortion
  when centroiding averaged spectra.
- Change logging to be consistent across modules.

### Deprecated

### Removed
- Removed deprecated `ScanProcessor.get_precursor_peak_for_product_scans`

### Fixed

### Security


## [v0.0.36] - 2022-01-15

### Added
1. Added `get_writer` to `ms_deisotope.output` which will guess which writer (and compression) to use from a file path
   or name of a file-like object.
2. Added an option to skip envelope decoding during `deserialize_deconvoluted_peak_set` as a theoretical optimization.

### Changed

### Deprecated

### Removed

### Fixed
1. Restore unpickling of `ThermoRawLoader`s

### Security


## [v0.0.35] - 2021-12-29

### Added
1. Added `toggle_peak_loading` to `PeakSetDeserializingMixin` to toggle loading and parsing of deconvoluted peak sets.
2. Refactored `MSFileLoader` into a callable class, `ms_deisotope.data_source.infer_type.FormatGuesser`. There is now a
   `ms_deisotope.output.ProcessedMSFileLoader` callable that can be used to determine the correct opener for processed
   mzML, MGF, and mzMLb files.

### Changed

### Deprecated

### Removed

### Fixed
1. Added a missing parameter `use_index` to `ProcessedMGFLoader` that made positional
   argument passing fail when unpickling instances of this subclass.
2. Repaired `peak_collection` fused functions to optimize `DeconvolutedPeakSetIndexed` too.
3. Fixed `TheoreticalIsotopicPattern.total` using an uninitialized accumulator.

### Security


## [v0.0.33] - 2021-12-05

### Added
1. Added an example script for splitting a FAIMS file into separate mzML files.
2. Added an example script for ion mobility and charge state deconvolving Waters HDMSe data.
3. Added Waters MassLynx SDK bindings for reading Waters .RAW folders directly on Windows
   when the SDK C library is present.
4. Added a new command `ms-index maintenance register-waters-masslynx` to register the location
   of the Waters MassLynx SDK if it is not installed in the expected location.
5. Improved documentation of `mzMLb` reading and writing.

### Changed

### Deprecated

### Removed

### Fixed

### Security


## [v0.0.32]

### Added

### Changed
1. Added an extra `ms_level` index to Thermo readers.
2. PROXI server now uses appropriate JSON types in its USI response.
3. `ProcessedMGFLoader`'s iterator behavior will bind generated `ProcessedScan` objects, keeping
   parity with `ProcessedMzMLLoader`.
4. The `mzml`, `mzmlb` and `mgf` commands for `ms-index` can now handle reading mzML files from
   `STDIN`.

### Deprecated

### Removed

### Fixed

### Security


## [v0.0.31]

### Added

### Changed

### Deprecated

### Removed

### Fixed
1. `PeakSetMethods` facet `__call__` properly detects peak modes.

### Security


## [v0.0.30] - 2021-9-2

### Added
1. Added `ms_deisotope.clustering.alignment` to provide a home for `SpectrumAlignment` and `SpectrumAlignmentGraph`.

### Changed

### Deprecated

### Removed

### Fixed
1. `MzMLSerializer` correctly interprets centroid-only but not peak-picked `Scan` objects.
2. Using `_InterleavedGroupedScanIteratorImpl` on Py2 no longer uses a `deque` as it is missing `insert`.
3. `LCMSFeatureSetFit` now has `n_points` in the pure Python version.
4. `SpectrumAlignment` now works properly.
5. Fix error in `MzMLSerializer` when trying to serialize a spectrum that has been peak-picked but where
   the peak set was empty.

### Security


## [v0.0.29] - 2021-7-8

### Fixed
1. Fix specifying the precursor m/z, charge and polarity of a spectrum from CSV or via `data_source.memory.make_scan`.


## [v0.0.28] - 2021-7-7

### Added
1. Added string parsing and formatting to `IDFormat` for those terms which define a nativeID format.
2. Added an `id_format` property to `FileInformation` and `ScanFileMetadataBase` which retrieves the
   nativeID format or formats for a given file.

### Changed
1. The default behavior of `_InterleavedGroupedScanIteratorImpl`, the implementation of grouped iterators,
   when producing a new `ScanBunch` that any product scans whose precursor ID that has been in the product map
   for more than `self.ms1_buffering` productions will be added to the produced `ScanBunch` to prevent loss
   of information if iteration is interrupted.
2. `quick_index.index` now tries much harder to start from an MS1 scan if it can.


### Deprecated

### Removed

### Fixed

### Security

## [v0.0.27] - 2021-6-21

### Added
1. Added TMT11 to `ms_deisotope.qc.signature`.
2. Added `ms_deisotope.output.mzml.IonMobilityAware3DMzMLSerializer` to write processed 3D IMS-MS spectra
    from mobility frames with feature maps.
3. Added `ms_deisotope.output.mzml.ProcessedGeneric3DIonMobilityFrameSource` to read processed feature maps
   out from 3D IMS-MS spectra.
4. Made `Generic3DIonMobilityFrameSource` wrapper iterable and more sequence-like.
5. Added `ms-index ms1-spectrum-diagnostics` to the CLI to collect relatively low level MS1 spectrum metadata.

### Changed
1. Made `LCMSFeatureProcessor` consider fewer combinations of feature sets, lowering the upper bound on the
   combinations. Such cases that required this should be quite rare.
2. `IonMobilityFrame` is closer to being a first-class object instead of an immutable data wrapper.

### Deprecated

### Removed

### Fixed
1. Iterating over MS3-containing datasets in grouped mode will now properly group MS3 spectra with their MS2
   spectra. Applies for higher exponentiated MSn as well.
2. `TheoreticalIsotopicPattern.incremental_truncation` consistently respects its truncation threshold.

### Security


## [v0.0.26] - 2021-5-30

### Added
1. Add `ion_mobility_type` property to `ScanBase` to allow checking ion mobility type on the scan object itself.
2. Add a prototype `FAIMSFilter` and `FAIMSDemultiplexingIterator` to `data_source.query`.
3. Add a new examples set to the repository showing other ways to use the library's features.

### Changed
1. Ensure that any lingering MSn scans are flushed with the final MS1 scan when an interleaved scan iteration
   strategy is wrapping up.
2. Made `LCMSFeatureProcessor` substantially faster and more memory efficient during dependence graph solving
   by introducing a Cython implementation. The solver is still vulnerable to high density noise clusters slowing
   it down, but will no longer completely OOM when these are wide enough.

### Deprecated

### Removed

### Fixed

### Security


## [v0.0.25] - 2021-5-02

### Added
1. Added `has_array` method to `RawDataArrays` and `RawDataArrays3D` to semantically query whether an
    array collection has an array of a particular type (like a flavor of ion mobility) without iteratively
    probing.
2. Added `ms_deisotope.data_source.scan.mobility_frame.Generic3DIonMobilityFrameSource` which can serve a
   ion mobility-aware data structure when there is a profile-mode m/z, intensity and ion mobility array
   in a single spectrum like the one produced by MSConvert's combine-ion-mobility-spectra option.

### Changed
1. `Scan.plot` will automatically call `pick_peaks` if `not self.is_profile and self.peak_set is None` to avoid
   unexpectedly ending up without a plot in this common scenario.

### Deprecated

### Removed

### Fixed
1. Fixed Docker container build process to use latest released versions of libraries.

### Security



## [v0.0.24] - 2021-4-10

### Added

### Changed

### Deprecated

### Removed

### Fixed
1. Made `feature_fit.map_coord` helper structure properly comparable.

### Security


## [v0.0.23rc] - 2021-4-9

### Added


### Changed
1. Re-aliased `Processed*Deserializer` to `Processed*Loader` where `*` is MzML, MGF, or MzMLb.

### Deprecated

### Removed

### Fixed
1. The Waters SDK no longer uses absolute local imports.
2. The Waters driver registration code no longer introduces `None` into the list of paths and breaks on Py3.

### Security


## [v0.0.22] - 2021-4-4

### Added
1. `ms_deisotope.data_source` now exports `scan_from_csv` and `make_scan` helper methods.
2. `ms_deisotope.data_source.mzmlb.MzMLbLoader` is now available for reading mzMLb HDF5 files when `pyteomics.mzmlb` is available.
3. `ms_deisotope.output.MzMLbSerializer` is now available for writing mzMLb HDF5 files when `psims.mzmlb` is available. This is further
   exposed through `ms-index mzmlb` for CLI conversion. Expect this feature to undergo further evolution as the extended indices used
   for other features may also be stored in mzMLb as extra datasets.

### Changed
1. `Scan` objects explicitly are not hashable.
2. Changed the `delimiter` argument of `ms_deisotope.data_source.text.scan_from_csv` to be a regular expression to handle
   arbitrary whitespace delimiters, and added an optional `skiprow` argument to allow you to skip headers in the all-too-common
   text spectrum exports that spectrum viewers provide.
3. Removed the resampling API from `RunningWeightedAverage`, improving memory efficiency.

### Deprecated

### Removed

### Fixed
1. Made random-access gzip compressor interface `idzip` compatible with Py3 buffered IO streams.
2. Added an additional flag `-D` to indicate to `ms-index metadata-index` that the input file is an
   a deisotoped and charge state deconvolved mzML to extract additional fields.
3. Fixed an interface error in `ms_deisotope._c.feature_map.processor`.

### Security


## [v0.0.21]

### Added

### Changed
1. All grouped scan iterators now use an interleaved iteration strategy to handle interleaving of MS1 scans that are
   not the precursor of a subsequent series of MSn scans. When using `start_from_scan`, product scans which follow the
   first MS1 but not actually produced from that first MS1 scan will be included in the product scan list for that
   `ScanBunch`.
2. `guess_type` now includes the input object in the error message when it fails to locate a loader type to reduce the
   amount of blind guessing at the problem.

### Deprecated

### Removed

### Fixed
1. Retrieving the precursor scan is now safer with `thermo_raw_net.ThermoRawLoader`
2. `ms-index spectrum-clustering` now remembers source files when scans are loaded in memory
3. Use a more thorough precursor scan determination algorithm for `ThermoRawLoader` implementations
4. Fixed `ms-deisotope`'s sequential scan numbering when there is scan interleaving

### Security


## [v0.0.19]  - 2020-11-24

### Added

### Changed
1. Make `peak_set.merge` able to merge more than 2 `DeconvolutedPeakSet` instances

### Deprecated

### Removed

### Fixed
1. Always infer the precursor scan relationship in ThermoRaw readers using the "Master Scan Number" header

### Security

## [v0.0.18]  - 2020-10-09
### Added
0. Added `CHANGELOG.md`
1. Added the option for `ProcessedMzMLLoader` to detect and handle deconvoluted peak lists from other tools and to load
   them without crashing while looking for `ms_deisotope`-specific extra data. This behavior can be customized by overriding
   the `deserialize_external_deconvoluted_peak_set` method in a derived class.
2. Added a `time` attribute to `RandomAccessScanSource`-based readers which lets you use index and slice notation with time
   values (`float`s) using the standard minute unit common across `ms_deisotope`.
3. Made `ms_deisotope.data_source.mgf.MGFLoader`, `ms_deisotope.peak_set.decharge`, and `ms_deisotope.processor.process` part
   of the top-level module's API.
4. Added `RawDataArray.size` property to report the length of the m/z and intensity array.
5. Added `default_precursor_ion_selection_window` (`-D`) to the `ms-deisotope` command line tool.
6. Add drift time aware variant of `DeconvolutedPeak`, `IonMobilityDeconvolutedPeak` which has a `drift_time` attribute.
7. Added more documentation on the `truncate_after` parameter for deconvolution.

### Changed
1. The `_reindex` method of `DeconvolutedPeakSet` has been renamed `reindex` to reflect that it should be part of the type's
   public API since user code may need to call it if they apply some new transformation to the list of peaks.
2. The `composition_list` argument to `HybridAveragineCompositionListPeakDependenceGraphDeconvoluter` is now optional.
3. The `incremental_truncation` option for `deconvolute_peaks` and `ExhaustivePeakSearchDeconvoluterBase`-based strategies
   now apply truncation to *all* fits, not just those passing the initial full-width fit. This required more invasive changes
   to the implementations of `AveragineDeconvoluterBase` and `MultiAveragineDeconvoluterBase` but is now more consistent with the
   original intent behind `incremental_truncation`.
4. `Scan.clear` now takes a `full` parameter which will discard the `peak_set`, `deconvoluted_peak_set`, and `product_scans` attributes'
   data. This more aggressively frees memory.
5. `Scan.average` and `Scan.average_with` now skip scans with empty signal arrays.
6. The `ScanProcessor.deconvolute_precursor_scan` method now explicitly requests `product_scans` be passed as an argument. It also extracts
   co-isolating precursor ions from the isolation window even if a deconvoluted peak is not found for the precursor ion.
7. `MzMLSerializer` now propagates all `PrecursorInformation.annotations` as parameters.
8. `ScanProcessor` now reports coisolation even when the precursor peak is not found (though it still omits it if an unacceptable solution is reported).

### Fixed
1. When using the `CompositionList`-based deconvoluters with a mass shift, the theoretical isotopic pattern will now have the correct
   monoisotopic m/z.
2. `IntervalTreeNode.overlaps` is now consistent with `SpanningMixin.overlaps`.
3. `ms-deisotope` CLI now properly handles the case when the MSn average option `-an` is not passed.
4. `ms-deisotope` CLI now properly builds the extended scan index when processing a file with only MSn spectra

---

<!-- Links -->
[Keep a Changelog]: https://keepachangelog.com/
[Semantic Versioning]: https://semver.org/

<!-- Versions -->
[Unreleased]: https://github.com/mobiusklein/ms_deisotope/compare/v0.0.46...HEAD
[Released]: https://github.com/mobiusklein/ms_deisotope/releases
[v0.0.46]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.46
[v0.0.45]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.45
[v0.0.44]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.44
[v0.0.43]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.43
[v0.0.42]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.42
[v0.0.40]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.40
[v0.0.39]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.39
[v0.0.38]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.38
[v0.0.37]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.37
[v0.0.36]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.36
[v0.0.35]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.35
[v0.0.33]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.33
[v0.0.32]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.32
[v0.0.31]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.31
[v0.0.30]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.30
[v0.0.29]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.29
[v0.0.28]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.28
[v0.0.27]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.27
[v0.0.26]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.26
[v0.0.25]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.25
[v0.0.24]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.24
[v0.0.23]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.23rc
[v0.0.22]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.22
[v0.0.21]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.21
[v0.0.19]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.19
[v0.0.18]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.18
[v0.0.16]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.16