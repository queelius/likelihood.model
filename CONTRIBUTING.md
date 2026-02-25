# Contributing to likelihood.model

Thank you for your interest in contributing to likelihood.model!

## Reporting Bugs

Open an issue at <https://github.com/queelius/likelihood.model/issues> with:

- A minimal reproducible example
- Output of `sessionInfo()`
- What you expected vs. what happened

## Suggesting Features

Open an issue describing the use case and how it relates to the likelihood model
framework. If proposing a new distribution or observation type, explain whether
it fits as a `likelihood_name()` wrapper, a `likelihood_contr_model` contribution,
or a standalone model class.

## Pull Requests

1. Fork the repository and create a feature branch
2. Follow the existing code style (S3 generics, roxygen2 documentation)
3. Add tests for new functionality in `tests/test.R`
4. Run `devtools::check()` and ensure no errors or warnings
5. Update documentation with `devtools::document()` if you changed any roxygen comments
6. Submit the PR with a clear description of the change

## Development Setup

```r
# Install dependencies
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::install_deps(dependencies = TRUE)

# Load for development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Code of Conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
