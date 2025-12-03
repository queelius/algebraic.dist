# algebraic.dist 0.1.0

* Initial CRAN release.

## Features

* Core distribution types: `normal`, `mvn`, `exponential`, `empirical_dist`
* Expression distributions (`edist`) for lazy composition of distributions
* Algebraic operations (`+`, `-`) on distributions with automatic simplification
* Support classes: `finite_set`, `interval` for representing distribution domains
* Generic methods: `sampler`, `mean`, `vcov`, `density`, `cdf`, `params`
* Monte Carlo estimation for `expectation`, `conditional`, and `rmap` operations
