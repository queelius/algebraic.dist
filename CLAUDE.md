# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`algebraic.dist` is an R package that provides an algebra over distributions (random elements). It enables operations on distribution objects, composition of distributions, and lazy evaluation of distribution expressions.

### Key Dependencies

- **mvtnorm**: For multivariate normal distribution functions (dmvnorm, pmvnorm, rmvnorm)
- **R6**: For support classes (finite_set, interval)
- **stats**: For standard distribution functions and statistical methods

## Key Architecture

### Class System Design

This package uses **S3 for distributions** and **R6 for support objects**:
- Distribution classes (`dist`, `normal`, `mvn`, `edist`, etc.) use S3 for lightweight polymorphism
- Support classes (`finite_set`, `interval`) use R6 for encapsulated state and methods

### Core Class Hierarchy

The package uses S3 classes with the following inheritance structure:

- **`dist`**: Base class for all distributions
  - **`univariate_dist`**: Single-dimensional distributions
    - `normal`: Univariate normal distribution
    - `exponential`: Exponential distribution
  - **`multivariate_dist`**: Multi-dimensional distributions
    - `mvn`: Multivariate normal distribution
    - `empirical_dist`: Empirical/bootstrap distribution from data
  - **`edist`**: Expression distributions (lazy composition)
    - Wraps expressions like `x + y` where x, y are distributions
    - Enables algebraic operations: `+`, `-` on dist objects
    - Supports deferred evaluation and potential simplification

### Distribution Types

- **`continuous_dist`**: Continuous distributions (normal, mvn, exponential)
- **`discrete_dist`**: Discrete distributions (empirical_dist)

### Support Classes (R6)

Support objects represent the set of values a distribution can realize. All supports must implement:
- `has(x)`: Check if value(s) are in the support
- `infimum()`: Get the infimum (minimum)
- `supremum()`: Get the supremum (maximum)
- `dim()`: Get the dimension

Implementations:
- **`finite_set`**: Finite discrete sets (used for empirical distributions)
- **`interval`**: Continuous intervals with open/closed bounds (used for continuous distributions)

## Key Generic Methods

All distribution objects implement these core generics (defined in `R/generic_dist.R`):

- `sampler(x)`: Returns a function that generates n samples
- `mean(x)`, `vcov(x)`: Statistical moments
- `params(x)`: Distribution parameters
- `density(x)`, `cdf(x)`: Probability functions
- `expectation(x, g)`: Expected value of function g under distribution x
- `marginal(x, indices)`: Marginal distributions
- `conditional(x, P)`: Conditional distributions given predicate P
- `rmap(x, g)`: Apply function g to distribution x
- `sup(x)`: Support of the distribution
- `simplify(x)`: Simplify distribution expressions (stub for future work)

## Development Commands

### Quick Start Workflow

```r
# After making changes to R code:
devtools::document()        # Update documentation
devtools::install()         # Install locally
devtools::check()           # Run R CMD check

# If you modified README.Rmd:
rmarkdown::render("README.Rmd")

# If you modified vignettes or want to update pkgdown site:
pkgdown::build_site()
```

### Building and Documentation

```r
# Generate documentation from roxygen2 comments
devtools::document()

# Build package tarball
devtools::build()

# Install locally (automatically runs document() first)
devtools::install()

# Check package (R CMD check)
devtools::check()

# Load package for interactive development
devtools::load_all()
```

### Testing

**IMPORTANT**: This package currently has **no tests directory**.

When creating the test infrastructure:

```r
# Create tests directory structure
usethis::use_testthat()

# Create a test file for a specific module
usethis::use_test("empirical_dist")

# Run all tests
devtools::test()

# Run test coverage
covr::package_coverage()

# Run a single test file
testthat::test_file("tests/testthat/test-empirical_dist.R")
```

Per user's global instructions, always:
1. Run unit and integration tests after changes
2. Think up new tests for new features
3. Run test coverage to identify gaps

### Vignettes and Website

```r
# Build vignettes
devtools::build_vignettes()

# Build pkgdown site (GitHub Pages)
pkgdown::build_site()
```

### README

README.md is generated from README.Rmd:

```r
# Rebuild README
rmarkdown::render("README.Rmd")
```

## Important Implementation Notes

### Data Structure Conventions

**Critical**: Throughout this package, data matrices follow a consistent convention:
- **Rows = observations (samples)**
- **Columns = dimensions (variables)**

This applies to:
- `empirical_dist` constructor input
- Sampler outputs for multivariate distributions
- `obs(x)` return values
- All internal data storage

Example:
```r
# 100 samples from a 3-dimensional distribution
samples <- sampler(my_mvn)(100)  # Returns 100×3 matrix
dim(samples)  # [1] 100   3
```

### Expression Distributions (`edist`)

Located in `R/edist.R` and `R/algebra.R`:
- `edist` objects store an expression `e` and a list of distributions `vars`
- When sampled, they evaluate the expression against samples from component distributions
- Created via algebraic operations: `x + y` where x, y are `dist` objects
- Class names are dynamically generated to include the expression and component types (e.g., `"x_+_y_normal_exponential"`)
- **Simplification is planned but not yet implemented** (see `R/algebra.R:40-49`)
- `simplify.edist()` currently just returns the object unchanged

### Empirical Distributions

Located in `R/empirical_dist.R`:
- Constructed from data matrices where **rows = observations, columns = dimensions**
- Methods like `conditional` and `rmap` filter/transform the underlying data directly
- Support is represented as a `finite_set` of observed values
- Can be univariate or multivariate depending on input data structure
- Inherits from both `discrete_dist` and either `univariate_dist` or `multivariate_dist`

### Monte Carlo Approximations

Default fallback for generic `dist` operations:
- Generic `dist` objects use MC estimation for `conditional()` and `rmap()`
- **Default sample size**: `n = 10000` (see `R/dist.R`)
- `expectation()` accepts a `control` list parameter:
  - `control$n`: Number of MC samples (default: 10000)
  - `control$compute_stats`: Whether to compute CIs (default: FALSE)
  - `control$alpha`: Significance level for CI (default: 0.05)
- When `compute_stats=TRUE`, returns list with `value`, `ci`, and `n`
- Uses CLT for confidence intervals in `expectation_data()` (see `R/utils.R`)

### Sampler Pattern

All distributions implement `sampler(x)` which returns a **function** (not a sample):
```r
# sampler returns a function
samp_fn <- sampler(my_dist)

# Call the returned function to generate n samples
samples <- samp_fn(n = 100)

# For univariate: returns vector of length n
# For multivariate: returns n×d matrix (rows = observations)
```

This pattern enables:
- Lazy evaluation (sampler is created once, called multiple times)
- Consistent interface across all distribution types
- Composition in `edist` (see `sampler.edist()` in `R/edist.R:85-119`)

### Future Direction

Per README.Rmd, planned features include:

- Automatic simplification of composed distributions (e.g., sum of normals → normal)
- Lazy expression trees that simplify to closed forms when possible
- Limiting distributions (see `limit.edist()` stub in `R/algebra.R:60-65`)
- Automatic deduction of moments and properties
- Moving code from `algebraic.mle` and `likelihood.model` packages to this package
