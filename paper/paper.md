---
title: 'algebraic.dist: An Algebra over Probability Distributions in R'
tags:
  - R
  - probability distributions
  - distribution algebra
  - Monte Carlo methods
  - statistical computing
authors:
  - name: Alexander Towell
    orcid: 0000-0001-6443-9897
    affiliation: 1
affiliations:
  - name: Southern Illinois University Edwardsville
    index: 1
date: 24 February 2026
bibliography: paper.bib
---

# Summary

`algebraic.dist` is an R package that treats probability distributions as
first-class algebraic objects. Distributions---normal, exponential, multivariate
normal, and empirical---are represented as S3 objects with a unified interface
for sampling, density evaluation, cumulative distribution functions, moments,
expectations, conditioning, and functional mappings. Algebraic operators (`+`,
`-`) applied to distribution objects produce *expression distributions*
(`edist`), lazy objects that capture both the algebraic expression and its
component distributions. When a mathematical identity applies---for example, the
sum of independent normal random variables is itself normal---the expression is
automatically simplified to the closed-form result. When no closed form exists,
all operations fall back to Monte Carlo estimation with Central Limit Theorem
based confidence intervals. The package is available on CRAN and serves as the
foundation layer of an ecosystem of R packages for maximum likelihood estimation
and reliability engineering [@algebraic.mle; @likelihood.model].

# Statement of Need

Researchers who build statistical models by composing simpler distributional
components---such as sums, differences, or transformations of random
variables---frequently need to propagate distributional properties through those
compositions. In reliability engineering, for instance, the lifetime of a series
system is the minimum of component lifetimes; in Bayesian updating, posterior
distributions arise as transformations of prior and likelihood. R provides
excellent infrastructure for individual distributions through the `d`/`p`/`q`/`r`
function families in the `stats` package [@R], but no built-in mechanism for
composing distributions as objects and reasoning about the result.

`algebraic.dist` targets researchers and practitioners in statistics, reliability
engineering, and simulation science who need to compose distributions
programmatically, propagate moments and densities through algebraic operations,
and obtain either exact closed-form results or Monte Carlo approximations with
quantified uncertainty. The package fills the gap between low-level distribution
functions and high-level probabilistic programming by providing a lightweight
algebraic layer where distributions are values that can be combined, simplified,
and queried.

# State of the Field

Several R packages represent probability distributions as objects. The `distr`
family of packages [@Ruckdeschel2006; @Ruckdeschel2014] provides a comprehensive
S4 class hierarchy with arithmetic operations on distributions, automatic
generation of density, CDF, and quantile functions via FFT-based convolution, and
extensions for expectations and statistical models. Its breadth is a strength but
also a source of complexity: the S4 dispatch system imposes a steep learning
curve, and the heavy dependency chain (over ten interrelated packages) can be
prohibitive for lightweight applications.

`distr6` [@Sonabend2021] reimplements much of the `distr` philosophy using R6
classes, adding over 50 built-in distributions, decorators for exotic methods such
as hazard functions and anti-derivatives, and support for truncation, mixture, and
product distributions. Its focus is on providing a complete object-oriented
interface; algebraic composition and automatic simplification to closed forms are
not central design goals.

`distributions3` [@Hayes2024] takes a tidyverse-compatible approach, representing
distributions as S3 objects with generics `pdf()`, `cdf()`, `quantile()`, and
`random()` designed for pedagogical clarity. The `distributional` package
[@OHaraWild2024] stores distributions in a vectorised `vctrs` format for
integration with tibbles and tidy workflows. The `distionary` package
[@Coia2026], part of the probaverse ecosystem, provides distribution building
blocks with automatic property derivation from partial specifications.

The `actuar` package [@Dutang2008] provides compound and aggregate distributions
for actuarial science, including convolutions and mixture operations, but targets
the actuarial loss-modeling domain rather than general-purpose algebraic
composition with symbolic simplification.

`algebraic.dist` occupies a distinct position in this landscape. Rather than
maximizing the catalogue of supported distributions or optimizing for tidy data
integration, it focuses on the *algebraic structure*: composing distributions via
operations, automatically simplifying to closed forms when identities apply, and
falling back to Monte Carlo estimation otherwise. The expression distribution
(`edist`) abstraction---which captures both the symbolic expression and its
component distributions as a single lazy object---is unique to this package and
enables extensible simplification rules. Additionally, `algebraic.dist` serves as
foundation infrastructure for a chain of downstream packages in reliability
engineering and maximum likelihood estimation, a role that shaped its minimal
dependency footprint (only `stats`, `mvtnorm`, and `R6`).

# Software Design

The package is organized around two class systems. Distribution objects use S3,
providing lightweight polymorphism with a flat class hierarchy: `dist` at the
root, with `univariate_dist` and `multivariate_dist` as dimensional subtypes and
`continuous_dist` and `discrete_dist` as measure-theoretic subtypes. Concrete
types---`normal`, `exponential`, `mvn`, and `empirical_dist`---inherit
from the appropriate combination. Support objects (`interval` and `finite_set`)
use R6 for encapsulated state, representing the domain of each distribution with
operations for membership testing, infimum, and supremum.

A key design decision is the *sampler pattern*: `sampler(x)` returns a
*function* rather than a sample. This closure captures the distribution's
parameters and can be called repeatedly with different sample sizes, enabling lazy
evaluation and composition. Expression distributions exploit this pattern
directly: `sampler.edist` constructs a closure that samples from each component
distribution and evaluates the captured R expression against those samples.

Algebraic operations on distributions are implemented as S3 methods on the `+`
and `-` operators for the `dist` class. Each operation first constructs an `edist`
wrapping the expression and its operands, then immediately calls `simplify()`.
The `simplify.edist` method dispatches on the operator and operand types: if both
operands are `normal`, the sum is reduced to a `normal` with the appropriate mean
and variance. If no simplification rule matches, the `edist` is returned
unchanged and subsequent operations (sampling, moments) use Monte Carlo
estimation. This design makes the simplification system extensible---new rules
for additional distribution families can be added without modifying existing code.

When closed-form methods are unavailable, the `expectation` generic falls back to
Monte Carlo estimation with configurable sample size and CLT-based confidence
intervals. For univariate continuous distributions, numerical integration via
`stats::integrate` is used instead, exploiting the known support bounds from the
distribution's `sup()` method.

# Research Impact Statement

`algebraic.dist` is published on CRAN (version 0.1.0) and available through
r-universe at <https://queelius.r-universe.dev>. It serves as the foundation
layer for a suite of packages addressing maximum likelihood estimation
(`algebraic.mle`, also on CRAN), Fisherian likelihood modeling
(`likelihood.model`), composable MLE solvers (`compositional.mle`), and a
reliability engineering chain for masked failure data analysis (`flexhaz`,
`serieshaz`, `maskedhaz`, `maskedcauses`). These downstream packages depend on the
distribution abstraction, sampler pattern, and support algebra provided by
`algebraic.dist`. As the ecosystem matures, the package is expected to serve as
shared infrastructure for simulation studies and Monte Carlo experiments in
reliability engineering and survival analysis research.

# AI Usage Disclosure

Claude (Anthropic) was used as a drafting assistant for this manuscript. The
author provided the package design, code, and technical content; the AI assisted
with prose structuring and JOSS formatting. All claims about the software, its
design rationale, and its relationship to other packages were reviewed and
verified by the author. No AI-generated code is included in the `algebraic.dist`
package itself.

# Acknowledgements

The author thanks the CRAN maintainers for their review of the package, and the
R community for the `mvtnorm` [@Genz2009], `R6` [@Chang2021], and `stats`
packages on which `algebraic.dist` depends.

# References
