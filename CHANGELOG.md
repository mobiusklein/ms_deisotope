# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][Keep a Changelog] and this project adheres to [Semantic Versioning][Semantic Versioning].

## [Unreleased]

### Added

### Changed

### Deprecated

### Removed

### Fixed
1. Retrieving the precursor scan is now safer with `thermo_raw_net.ThermoRawLoader`
2. `ms-index spectrum-clustering` now remembers source files when scans are loaded in memory

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
[Unreleased]: https://github.com/mobiusklein/ms_deisotope/compare/v0.0.19...HEAD
[Released]: https://github.com/mobiusklein/ms_deisotope/releases
[0.0.16]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.16
[0.0.18]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.18
[0.0.18]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.19