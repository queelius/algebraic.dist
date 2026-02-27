## R CMD check results

0 errors | 0 warnings | 1 note

The single NOTE is:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Alexander Towell <lex@metafunctor.com>'

  This is a version update (0.1.0 -> 0.9.1) of an existing CRAN package.

## Changes in this version

Major feature release adding:
- 12 new distribution types (gamma, Weibull, beta, chi-squared, uniform,
  log-normal, Poisson, mixture, plus expression and realized distributions)
- Distribution algebra with automatic simplification (20+ rules)
- MVN closed-form conditioning and affine transforms
- Mixture distributions with Bayesian conditional updates
- Limiting distributions (CLT, LLN, delta method, normal approximation)
- Full documentation: all 67+ exports have @return and @examples
- 3 vignettes: Examples, Distribution Algebra, Multivariate/Mixture

## Test environments

* Local: Ubuntu 24.04.3 LTS, R 4.3.3
* GitHub Actions (R-CMD-check):
  - ubuntu-latest (R oldrel-1): success
  - ubuntu-latest (R release): success
  - ubuntu-latest (R devel): success
  - windows-latest (R release): success
  - macos-latest (R release): success

## Downstream dependencies

`algebraic.mle` (on CRAN) depends on this package. We have verified
the update does not break `algebraic.mle`.
