## Resubmission

This is a resubmission in response to CRAN feedback.

The following changes were made:
* revised the `Description` field in `DESCRIPTION` to cite the methodological reference in CRAN-recommended form with authors, year, and DOI
* removed unnecessary `\donttest{}` wrappers where feasible
* reduced example sizes so runnable examples are checked automatically
* refined examples and documentation accordingly

## Test environments
* local macOS (Apple silicon M3 Max), R 4.5.2
* win-builder (R-devel)
* R-hub: linux, macos-arm64, windows

## R CMD check results
0 errors | 0 warnings | 0 notes

* win-builder shows only the standard "New submission" NOTE.
* A possible misspelling note for `debiased` may appear in `DESCRIPTION`; this is intended technical terminology.

## Downstream dependencies
There are currently no downstream dependencies.
