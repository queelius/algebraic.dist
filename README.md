
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.dist`

<!-- badges: start -->

[![GPL-3
License](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

<!-- summary-start -->
An algebra over distributions (random elements).
<!-- summary-end -->

**Tags**:
<!-- tags-start -->
multivariate distributions, multivariate normal distribution, multivariate empirical distribution, data generating process, R, data-science, statistics, inference, likelihood-models, probability-theory
<!-- tags-end -->

**Table Of Contents**

- [R package: `algebraic.dist`](#r-package-algebraicdist)
   - [Installation](#installation)
   - [About](#about)

## GitHub Pages Documentation

The GitHub documentation can be viewed [here](https://queelius.github.io/algebraic.dist/).

See the vignette [algebraic.dist: Example](https://queelius.github.io/algebraic.dist/articles/example.html)
for a quick introduction to the package.

## Installation

You can install the development version of `algebraic.dist` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("queelius/algebraic.dist")
```

## About

The R package `algebraic.dist` provides an algebra over distributions.
It's not fully-formed yet, but I plan on using it for a lot of my future work.
For instance, I'll move a lot of the code in `algebraic.mle` and 
`likelihood.model` to this package.

After that, I want to experiment with using the `algebraic.dist` to do the
following:

- Compose distributions such that operations over distributions generate
other known distributions.

  There are a lot of well-known compositions, such as
the exponential distribution being the minimum of independent exponential distributions, or the sum of independent normal
distributions being a normal distribution, but there is a very large space of
possible compositions that are not as well-known or well-studied that I want to
explore.

- Let people use an R expression to lazily compose functions of distributions.
Simplifying a distribution expression will generate a most simple R expression
that represents the same distribution.

Sometimes, this may result in a simple close-form distribution, like a 
multivariate normal distribution, but in other cases it may result in a
(hopefully simpler) expression that composes multiple distributions and
operations over them.

- With these R expressions that represent distributions, we can define more
operations, like taking the limiting distribution of a sequence of
distributions, say $\lim_{n \to \infty} \frac{1}{n} \sum_{i=1}^n X_i$, which is
of normal by the central limit theorem.

- Deduce various properties of these distributions, such as their moments,
variances, etc. Sometimes, this may require numerical integration or Monte
Carlo methods, but if the expression simplifies to a known distribution, then
we can use the known properties of that distribution.

I have a lot of this code in place in C++, but I want to
re-implement it in R so that it's more accessible to others. I may also
implement some of the more interesting compositions in C++ and expose them to R
via Rcpp, but I'm not sure yet. I use a lot of templates and metaprogramming in
C++, and I'm not sure how well that will translate to Rcpp.

