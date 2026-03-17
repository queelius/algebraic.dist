## R CMD check results

0 errors | 0 warnings | 0 notes

## Major release (v0.9.1 → v1.0.0)

API stabilization release. Key changes:

- Fixed mvn CDF spurious warning (pmvnorm error tolerance check)
- Added `nparams` default method for all parametric distributions
- Added `dim.edist` for expression distributions
- Added `hazard`/`surv` fallback methods for all continuous distributions
- Standardized closure signatures (removed parameter overrides, added `...`)

## Coordinated submission

This is part of a coordinated 6-package submission. All packages are
maintained by me. Updated versions being submitted simultaneously:

- algebraic.dist 1.0.0 (this package)
- algebraic.mle 2.0.2
- likelihood.model 1.0.0
- compositional.mle 2.0.0
- flexhaz 0.5.1
- maskedcauses 0.9.3

## Test environments

* local Ubuntu 24.04, R 4.3.3
