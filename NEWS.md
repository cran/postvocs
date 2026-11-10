# postvocs 0.2.6

## New Features
- `filter_by_frequency()` gains a new argument `remove_blank_compounds` (default `FALSE`).
  - When `FALSE`, blank samples are excluded from frequency calculations, but compounds detected in blanks are retained. A new flag column `In_Blank` (`"*"` if detected in any blank sample, otherwise empty) is added to `abundance`, `round1`, `round2`, and `summary_table`. Blank sample columns are retained in `abundance` (placed before the frequency columns).
  - When `TRUE`, compounds detected in any blank sample are removed from all samples. `summary_table` keeps the `Removed_in_Blank` column (`"-"` for removed compounds) and computes their frequencies (overall and within each group) from the original data, excluding blank samples from the denominator. No `In_Blank` column is added in this mode.

## Enhancements
- `summary_table` in `filter_by_frequency()` no longer contains the `Removed_in_Blank` column when `remove_blank_compounds = FALSE`.
- Blank and internal standard samples are always excluded from frequency calculations, regardless of `remove_blank_compounds`.

## Documentation
- Updated `filter_by_frequency()` documentation to describe the new `remove_blank_compounds` argument, the `In_Blank` flag column, and the revised `summary_table` layout.

# postvocs 0.2.5

## Breaking Changes
- None.

## New Features
- None.

## Enhancements
- `annotate_compounds()` example now uses a user-provided library (`lib_source = "user"`) instead of webchem, eliminating network dependency and ensuring stable CRAN checks.
- Added pre-annotated example data (`inst/extdata/annotated_example.rds`) for use in `filter_by_frequency()` and `save_postvocs_results()` examples.
- `filter_by_frequency()` and `save_postvocs_results()` examples now load pre-annotated data, avoiding repeated `annotate_compounds()` calls.
- Added defensive checks for `prop_df` in `annotate_compounds()` to handle edge cases when `pc_prop()` returns unexpected results.

## Bug Fixes
- Fixed macOS CRAN check failure caused by `webchem` service unavailability in `annotate_compounds()` example.
- Removed commented-out code lines from all examples per CRAN policy.

## Documentation
- Updated examples to be self-contained and executable without network access.
- Added `inst/extdata/annotated_example.rds` as example data for downstream functions.

# postvocs 0.2.4

## Bug Fixes
- Removed all commented-out code lines from examples in `batch_process_gcms()`, `extract_peak_areas()`, and `process_gcms_txt()` per CRAN review. The previously commented examples are now uncommented and executable.
# postvocs 0.2.3

## Bug Fixes
- Replaced `\dontrun{}` with `\donttest{}` in all examples per CRAN policy.
- Removed commented-out code lines from examples.
- Removed default output paths from functions to comply with CRAN policies.

# postvocs 0.2.2

## Bug Fixes
- Removed redundant LICENSE file and simplified License field to `GPL (>= 3)` as requested by CRAN reviewer.

# postvocs 0.2.1

## Bug Fixes
- Fixed FileID matching in `batch_process_gcms()` to support exact filename mapping.

# postvocs 0.2.0

## Breaking Changes
- The `extract_cas_abundance()` alias has been removed. Please use `extract_peak_areas()` instead.
- `build_cas_abundance_matrix()` has been renamed to `build_cas_abundance()` with simplified interface.

## New Features
- Added `extract_peak_areas()` as a more flexible replacement for `extract_cas_abundance()`:
  - Supports data.frame input, file paths, and previous function outputs.
  - Can automatically detect and process all samples from a folder.
- Added `build_cas_abundance()` as a replacement for `build_cas_abundance_matrix()`:
  - Directly accepts output from `extract_peak_areas()`.
  - No automatic file saving; returns a data.frame for user control.
- Added `save_postvocs_results()` for unified saving of analysis results:
  - Supports CSV and XLSX formats.
  - Automatically adds single quote prefix to CAS columns to prevent Excel date conversion.
  - For XLSX output, multiple tables are saved as separate worksheets in one file.
- Added `filter_by_frequency()` for two-step occurrence frequency screening:
  - Removes compounds detected in blank samples.
  - Filters by total frequency and treatment-specific frequency.
  - Returns detailed summary tables with flags for each compound.

## Enhancements
- `annotate_compounds()` now supports Excel input and automatically cleans CAS numbers.
- `filter_by_frequency()` now accepts `annotate_compounds()` output directly.
- `save_postvocs_results()` now checks whether a CAS column already has a leading single quote before adding one, avoiding duplication.
- Improved documentation and examples for all functions.
- CAS numbers are automatically cleaned (spaces removed) during parsing.

## Bug Fixes
- Improved handling of column name detection in `extract_peak_areas()`.
- Enhanced error messages for better user guidance.

## Documentation
- Added comprehensive README.md with workflow examples.
- All functions now have complete roxygen2 documentation.
- Added example data files (`inst/extdata/txt/` with 13 sample .txt files and `SampleID.xlsx`) so all function examples are now self-contained and directly runnable.

# postvocs 0.1.0

- Initial development version.
- Basic package structure with core functions for GC-MS data processing.
- Internal development only; not released on CRAN.
