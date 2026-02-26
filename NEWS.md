# algebraic.dist 0.2.0

## New Distributions

* `gamma_dist(shape, rate)` — Gamma distribution with hazard/survival functions
* `weibull_dist(shape, scale)` — Weibull distribution with closed-form hazard
* `chi_squared(df)` — Chi-squared distribution with hazard/survival functions
* `uniform_dist(min, max)` — Uniform distribution on [min, max]
* `beta_dist(shape1, shape2)` — Beta distribution on (0, 1)
* `lognormal(meanlog, sdlog)` — Log-normal distribution with hazard/survival
* `poisson_dist(lambda)` — Poisson distribution with exact expectation via truncated summation
* `mixture(components, weights)` — Mixture distributions with law of total variance

## New Operators and Algebra

* `*.dist` — Scalar multiplication (`c * dist`, `dist * c`, `dist * dist`)
* `^.dist` — Power operator (`dist ^ n`)
* `Math.dist` — Group generic for `exp()`, `log()`, `sqrt()`, `abs()`, etc.
* `Summary.dist` — Group generic for `sum()`, `prod()`, `min()`, `max()`
* Extended `+.dist` and `-.dist` for numeric location shifts

## Simplification Rules

* `c * Normal(mu, v)` simplifies to `Normal(c*mu, c^2*v)`
* `c * Gamma(a, r)` simplifies to `Gamma(a, r/c)` for c > 0
* `c * Exponential(r)` simplifies to `Gamma(1, r/c)` for c > 0
* `Normal(mu, v) + c` simplifies to `Normal(mu+c, v)`
* `Gamma(a1, r) + Gamma(a2, r)` simplifies to `Gamma(a1+a2, r)` (same rate)
* `Exp(r) + Exp(r)` simplifies to `Gamma(2, r)` (same rate)
* `ChiSq(d1) + ChiSq(d2)` simplifies to `ChiSq(d1+d2)`
* `Poisson(l1) + Poisson(l2)` simplifies to `Poisson(l1+l2)`
* `Normal(0,1)^2` simplifies to `ChiSq(1)`
* `exp(Normal(mu, v))` simplifies to `LogNormal(mu, sqrt(v))`
* `log(LogNormal(ml, sl))` simplifies to `Normal(ml, sl^2)`
* `min(Exp(r1), ..., Exp(rk))` simplifies to `Exp(sum(r))`

## New Infrastructure

* `realize()` generic — materialize any distribution to `empirical_dist` by sampling
* Auto-fallback methods for `edist`: `cdf`, `density`, `sup`, `conditional`, `rmap`, `inv_cdf`
* `countable_set` R6 class for countably infinite support (Poisson)
* `inv_cdf.empirical_dist` — quantile function for empirical distributions

## Improvements

* Informative error messages in all constructors (replaced `stopifnot`)
* `format()` methods for all distribution types
* Standardized `print()` methods delegating to `format()`
* Fixed `vcov.exponential` — was returning `rate` instead of `1/rate^2`
* Fixed `sampler.edist` crash when n=1
* `conditional.empirical_dist` gives informative error on zero matches
* Zero-variance guard in `expectation_data()` CI computation

---

# algebraic.dist 0.1.0

* Initial CRAN release.

## Features

* Core distribution types: `normal`, `mvn`, `exponential`, `empirical_dist`
* Expression distributions (`edist`) for lazy composition of distributions
* Algebraic operations (`+`, `-`) on distributions with automatic simplification
* Support classes: `finite_set`, `interval` for representing distribution domains
* Generic methods: `sampler`, `mean`, `vcov`, `density`, `cdf`, `params`
* Monte Carlo estimation for `expectation`, `conditional`, and `rmap` operations
