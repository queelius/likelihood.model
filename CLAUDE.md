# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`likelihood.model` is an R package providing the core framework for likelihood-based statistical inference in the Fisherian tradition. It defines the `likelihood_model` concept (generic interface) with methods `loglik`, `score`, `hess_loglik`, `fim`, and provides infrastructure for MLE fitting, Fisherian inference, likelihood ratio tests, and bootstrap sampling.

The contribution-based model builder (`likelihood_contr_model`) and the named distribution wrapper (`likelihood_name`) have been extracted to the separate `likelihood.contr` package.

The package is designed to interoperate with the `algebraic.mle` package for maximum likelihood estimation and the `algebraic.dist` package for distributions.

## Development Commands

### Building and Testing
```r
# Install dependencies (core + suggested)
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))

# Install companion packages (needed for full functionality)
devtools::install_github("queelius/algebraic.mle")
devtools::install_github("queelius/algebraic.dist")

# Load package during development (fastest workflow)
devtools::load_all()

# Run tests
testthat::test_file("tests/test.R")

# Run test coverage analysis
covr::package_coverage()

# Build documentation from roxygen2 comments
devtools::document()

# Check package (CRAN standards)
devtools::check()

# Build package tarball
devtools::build()

# Install from local source
devtools::install()
```

### Documentation
```r
# Generate documentation from roxygen2 comments
devtools::document()

# Build vignettes
devtools::build_vignettes()

# Preview pkgdown site
pkgdown::build_site()
```

## Architecture

The package is organized in layers, reflected in the file naming convention.

### Core Concept: Likelihood Model

The package defines `likelihood_model` as a **concept** (not a concrete class) that objects must satisfy. To be a likelihood model, an object must:
- Have class `"likelihood_model"` in its inheritance chain
- Implement `loglik(model, ...)` method that returns a log-likelihood function
- Optionally implement `score(model, ...)` and `hess_loglik(model, ...)` (defaults use numerical differentiation via `numDeriv`)

Generic functions that operate on likelihood models:
- `fit.likelihood_model()`: MLE solver using `optim`, returns `fisher_mle` objects
- `sampler.likelihood_model()`: Bootstrap sampling distribution via `boot` package
- `lrt()`: Likelihood ratio test between nested models
- `fim()`: Fisher information matrix
- `assumptions()`: Retrieve model assumptions

### Example Implementation

**`exponential_lifetime`** (in `R/example-exponential.R`):
- Exponential distribution with optional right-censoring support
- Demonstrates closed-form MLE (`fit()` override bypasses `optim` entirely)
- Provides analytical score, hessian, and FIM
- Includes `rdata()` method for Monte Carlo validation
- Sufficient statistics: d (exact count), T (total time), MLE: lambda_hat = d/T

### Key Implementation Details

- **Method dispatch**: Uses S3 generics. Models add `"likelihood_model"` to class vector to inherit default methods
- **Numerical defaults**: `score.likelihood_model()` and `hess_loglik.likelihood_model()` use Richardson extrapolation with `r=6` via `numDeriv`
- **Optimization**: `fit.likelihood_model()` supports multiple `optim` methods (Nelder-Mead, BFGS, SANN, etc.). SANN useful for multimodal likelihoods but should be followed by local search
- **Bootstrap**: `sampler.likelihood_model()` uses `boot::boot()` with multicore parallelization support

### File Structure

**Core layer** (concept + infrastructure):
- `R/core-generics.R`: Concept definition, generics (`loglik`, `score`, `hess_loglik`, `fim`, etc.), default methods, `print.likelihood_model`
- `R/core-fit.R`: `fit` re-export, `fit.likelihood_model` (optim-based solver), `sampler.likelihood_model` (bootstrap)
- `R/core-fisher_mle.R`: `fisher_mle` and `fisher_boot` result objects with S3 methods (`coef`, `vcov`, `confint`, `se`, `aic`, `bic`, `summary`)
- `R/core-fisherian.R`: Fisherian inference (`support`, `relative_likelihood`, `likelihood_interval`, `profile_loglik`, `evidence`)
- `R/core-lrt.R`: Likelihood ratio test (`lrt()`)

**Example implementation**:
- `R/example-exponential.R`: `exponential_lifetime` with closed-form MLE and right-censoring

**Tests and vignettes**:
- `tests/test.R`: ~236 unit tests covering core framework and exponential_lifetime
- `vignettes/getting-started.Rmd`: Quick introduction with examples
- `vignettes/exponential-lifetime.Rmd`: Detailed exponential model coverage
- `vignettes/algebraic-mle-integration.Rmd`: Three-package ecosystem demonstration

## Important Notes

### Testing and Coverage
- `tests/test.R` contains ~236 tests covering all core infrastructure and the exponential_lifetime model
- When adding new features, write complete tests and run `covr::package_coverage()` to ensure coverage
- Focus on testing edge cases: censored data, empty data frames, parameter validation
- See vignettes for realistic examples to base tests on

### Adding New Specialized Models
Create a new class (like `exponential_lifetime`):
- Define S3 methods: `loglik.<class>`, `score.<class>`, `hess_loglik.<class>`
- Add `"likelihood_model"` to class vector to inherit default behavior
- Optionally override `fit.<class>` for closed-form MLEs (see `exponential_lifetime`)
- Optionally implement `fim.<class>` for analytical FIM
- Optionally implement `rdata.<class>` for Monte Carlo validation

For contribution-based models and standard distribution wrappers, use the `likelihood.contr` package.

### Performance Considerations
- **Numerical differentiation is expensive**: Default `score` and `hess_loglik` use `numDeriv` with Richardson extrapolation (`r=6`)
- Provide analytical derivatives whenever possible (see `example-exponential.R`)
- Always validate inputs: check for NULL, positive values where required, sufficient sample size

### Package Dependencies
- **Required**: `algebraic.mle` (>= 2.0.0), `generics`, `stats`, `numDeriv`, `boot`
- **Suggested**: `mvtnorm`, `algebraic.dist`, `testthat` (>= 3.0.0), `knitr`, `rmarkdown`
- **Roxygen2**: Version 7.3.3 for documentation generation

## Common Workflows

### Using the Exponential Lifetime Model
```r
# Create model
model <- exponential_lifetime("t")

# Fit (closed-form, no initial guess needed)
mle <- fit(model)(df)
coef(mle)
se(mle)
confint(mle)

# With right-censoring
model_c <- exponential_lifetime("t", censor_col = "status")
mle_c <- fit(model_c)(df_censored)
```

### Fisherian Inference
```r
# Support, relative likelihood, likelihood intervals
s <- support(mle, theta, df, model)
rl <- relative_likelihood(mle, theta, df, model)
li <- likelihood_interval(mle, df, model, k = 8)
prof <- profile_loglik(mle, df, model, param = 1)
ev <- evidence(model, df, theta1, theta2)
```

### Likelihood Ratio Test
```r
lrt_result <- lrt(
  null = null_model,
  alt = alt_model,
  data = df,
  null_par = params(null_result),
  alt_par = params(alt_result),
  dof = length(params(alt_result)) - length(params(null_result))
)
lrt_result$p.value
```
