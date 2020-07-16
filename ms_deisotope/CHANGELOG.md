# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog][Keep a Changelog] and this project adheres to [Semantic Versioning][Semantic Versioning].

## [Unreleased]
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

### Fixed
1. When using the `CompositionList`-based deconvoluters with a mass shift, the theoretical isotopic pattern will now have the correct
   monoisotopic m/z.

---

## [Released]

---

<!-- Links -->
[Keep a Changelog]: https://keepachangelog.com/
[Semantic Versioning]: https://semver.org/

<!-- Versions -->
[Unreleased]: https://github.com/mobiusklein/ms_deisotope/compare/v0.0.16...HEAD
[Released]: https://github.com/mobiusklein/ms_deisotope/releases
[0.0.16]: https://github.com/mobiusklein/ms_deisotope/releases/v0.0.16