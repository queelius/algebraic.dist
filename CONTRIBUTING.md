# Contributing to algebraic.dist

Thank you for considering contributing to `algebraic.dist`.

## How to Contribute

### Reporting Bugs

Open an issue at <https://github.com/queelius/algebraic.dist/issues> with:

- A minimal reproducible example
- Your R version (`sessionInfo()`)
- The package version (`packageVersion("algebraic.dist")`)

### Suggesting Features

Open an issue describing the feature, its use case, and how it fits the package's algebraic approach to distributions.

### Pull Requests

1. Fork the repository and create a feature branch from `main`.
2. Add tests for any new functionality (we use `testthat` edition 3).
3. Run `devtools::check()` and ensure there are no errors, warnings, or notes.
4. Submit a pull request describing your changes.

### Code Style

- Follow existing conventions: S3 for distribution classes, R6 for support objects.
- Document all exported functions with roxygen2.
- Use the sampler-returns-function pattern for new distribution types.

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating you agree to abide by its terms.
