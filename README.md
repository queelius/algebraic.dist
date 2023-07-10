
  - [R package: `algebraic.dist`](#r-package-algebraicdist)
      - [Installation](#installation)
      - [TODO](#todo)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.dist`

<!-- badges: start -->

<!-- badges: end -->

An algebra over distributions (random elements).

## Installation

You can install the development version of `algebraic.dist` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("queelius/algebraic.dist")
```

See the vignette [algebraic.dist:
Example](http://queelius.github.io/algebraic.dist/docs/articles/example.html)
for a quick introduction to the package.

## TODO

1.  Flesh out the vignettes to show how to use the package.

2.  Reincorporate the `edist` (expresion distribution), which is just a
    an R expression combined with an environment that has a dictionary
    for each of the variables in the expression where if a variable is a
    `dist`-type object, then we treat that variable like a random
    variable. When we, say, sample from the `edist` object, we’ll sample
    from the `dist`-type objects in the environment and then evaluate
    the expression with the sampled values.

3.  For taking the PDF, CDF, expectation, and so on, for now we just use
    the sampler and then put it into a `empirical_dist` object, which
    permits those kinds of operations, but there are clearly a lot of
    interesting things we can do, e.g., `min(x1, x2)` when `x1` and `x2`
    are exponential distributions is also an exponential distribution
    with a failure rate `lambda1 + lambda2`. There are a host of
    relations between various named distributions, and there are even
    more cases where the distribution, no matter how complex, has an
    asymptotically normal distribution.
