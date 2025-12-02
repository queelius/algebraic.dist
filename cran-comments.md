## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Alexander Towell <lex@metafunctor.com>'
  New submission

## Notes on print.dist method

This package defines a `print.dist` method for probability distribution objects.
The `stats` package also has a `print.dist` method for distance matrix objects
(from `stats::dist()`). Our implementation checks the object's class and defers
to `stats:::print.dist` via `NextMethod()` when the object is not one of our
probability distribution types, ensuring proper interoperability.

## Test environments

* Local: Ubuntu 24.04.3 LTS, R 4.3.3 (338 tests passing)
* GitHub Actions (R-CMD-check):
  - ubuntu-latest (R oldrel-1): success
  - ubuntu-latest (R release): success
  - ubuntu-latest (R devel): success
  - windows-latest (R release): success
  - macos-latest (R release): success
* R-hub:
  - linux (R-devel): success
  - windows (R-devel): success
  - macos-arm64 (R-devel): success
* win-builder:
  - R-devel: pending
  - R-release: pending
  - R-oldrelease: pending

## Downstream dependencies

This is a new package with no downstream dependencies.
