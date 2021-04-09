# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][Keep a Changelog] and this project adheres to [Semantic Versioning][Semantic Versioning].

## [Unreleased]

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


## [Released]

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
[Unreleased]: https://github.com/mobiusklein/ms_deisotope/compare/v0.0.22...HEAD
[Released]: https://github.com/mobiusklein/ms_deisotope/releases
[v0.0.22]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.22
[v0.0.21]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.21
[v0.0.19]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.19
[v0.0.18]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.18
[v0.0.16]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.16