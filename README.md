
- [R package: `algebraic.dist`](#r-package-algebraicdist)
  - [Installation](#installation)
  - [Overview](#overview)
    - [Distribution types](#distribution-types)
    - [Automatic simplification](#automatic-simplification)
    - [Multivariate operations](#multivariate-operations)
    - [Limiting distributions](#limiting-distributions)
  - [Quick example](#quick-example)
  - [Documentation](#documentation)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.dist`

<!-- badges: start -->

[![GPL-3
License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/queelius/algebraic.dist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/queelius/algebraic.dist/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/algebraic.dist)](https://CRAN.R-project.org/package=algebraic.dist)
<!-- badges: end -->

<!-- summary-start -->

An algebra over probability distributions: compose, transform, and
automatically simplify distribution expressions in R.
<!-- summary-end -->

<!-- tags-start -->

**Tags**: probability distributions, distribution algebra, automatic
simplification, multivariate normal, mixture models, CLT, delta method,
Monte Carlo, R, statistics <!-- tags-end -->

## Installation

Install the stable release from CRAN:

``` r
install.packages("algebraic.dist")
```

Or install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("queelius/algebraic.dist")
```

## Overview

`algebraic.dist` lets you build and manipulate probability distributions
as first-class R objects. Algebraic operations (`+`, `-`, `*`, `/`, `^`,
`exp`, `log`, `min`, `max`, …) on distribution objects automatically
simplify to closed-form distributions when mathematical identities
apply, and fall back to lazy Monte Carlo expressions (`edist`)
otherwise.

### Distribution types

Normal, exponential, gamma, Weibull, chi-squared, uniform, beta,
log-normal, Poisson, multivariate normal (MVN), mixture, and empirical
distributions.

### Automatic simplification

Over 20 built-in rules, including:

- `Normal + Normal` → `Normal`
- `Gamma + Gamma` (same rate) → `Gamma`
- `exp(Normal)` → `LogNormal`
- `Normal(0,1)^2` → `ChiSq(1)`
- `min(Exp, ..., Exp)` → `Exp`
- `c * Uniform(a,b)` → `Uniform`

When no rule matches, the result is a lazy `edist` that samples from its
components and evaluates the expression on demand.

### Multivariate operations

- **MVN conditioning**: closed-form Schur complement via `conditional()`
- **Affine transforms**: `affine_transform(x, A, b)` for exact linear
  maps
- **Mixture conditioning**: Bayesian weight updates via Bayes’ rule
- **Marginals**: exact for MVN and mixture distributions

### Limiting distributions

- `clt()` — Central Limit Theorem
- `lln()` — Law of Large Numbers
- `delta_clt()` — delta method for transformed means
- `normal_approx()` — moment-matching normal approximation

## Quick example

``` r
library(algebraic.dist)

# Sum of normals simplifies to a normal
x <- normal(1, 4)
y <- normal(2, 5)
z <- x + y
z
#> Normal distribution (mu = 3, var = 9)
```

``` r
# exp of a normal simplifies to log-normal
w <- exp(normal(0, 1))
w
#> Log-normal distribution (meanlog = 0, sdlog = 1)
```

``` r
# Gamma addition with matching rates
g <- gamma_dist(3, 2) + gamma_dist(4, 2)
g
#> Gamma distribution (shape = 7, rate = 2)
```

``` r
# CLT: the standardized sample mean converges to N(0, 1)
clt(normal(5, 4))
#> Normal distribution (mu = 0, var = 4)
```

## Documentation

The full documentation is available at
<https://queelius.github.io/algebraic.dist/>.

Vignettes:

- [Getting
  started](https://queelius.github.io/algebraic.dist/articles/example.html)
  — core distribution objects, sampling, and basic operations
- [Distribution
  algebra](https://queelius.github.io/algebraic.dist/articles/algebra.html)
  — simplification rules, limiting distributions, and the CLT/LLN/delta
  method
- [Multivariate
  operations](https://queelius.github.io/algebraic.dist/articles/multivariate.html)
  — MVN conditioning, affine transforms, and Gaussian mixture models
