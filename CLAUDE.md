# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`likelihood.model` is an R package for specifying and using likelihood models for statistical inference. It defines a generic framework where likelihood models implement key methods (`loglik`, `score`, `hess_loglik`) and provides two main implementations: `likelihood_contr_model` (for likelihood contributions with different observation types) and `likelihood_name_model` (for standard R distributions).

The package is designed to interoperate with the `algebraic.mle` package for maximum likelihood estimation and the `algebraic.dist` package for distributions.

## Development Commands

### Building and Testing
```r
# Install dependencies (core + suggested)
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "dplyr", "tibble"))

# Install companion packages (needed for full functionality)
devtools::install_github("queelius/algebraic.mle")
devtools::install_github("queelius/algebraic.dist")

# Load package during development (fastest workflow)
devtools::load_all()

# Run tests
testthat::test_file("tests/test.R")

# Run test coverage analysis (important - tests are incomplete)
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

### Core Concept: Likelihood Model

The package defines `likelihood_model` as a **concept** (not a concrete class) that objects must satisfy. To be a likelihood model, an object must:
- Have class `"likelihood_model"` in its inheritance chain
- Implement `loglik(model, ...)` method that returns a log-likelihood function
- Optionally implement `score(model, ...)` and `hess_loglik(model, ...)` (defaults use numerical differentiation via `numDeriv`)

Generic functions that operate on likelihood models:
- `fit.likelihood_model()`: MLE solver using `optim`, returns `algebraic.mle::mle_numerical` objects
- `sampler.likelihood_model()`: Bootstrap sampling distribution via `boot` package
- `lrt()`: Likelihood ratio test between nested models
- `fim()`: Fisher information matrix
- `assumptions()`: Retrieve model assumptions

### Main Implementations

**`likelihood_contr_model`** (R6 class in `R/likelihood_contr_model.R`):
- **Uses R6 classes** (not S3) for object-oriented programming with mutable state
- Handles heterogeneous observation types (exact, left-censored, right-censored, interval-censored)
- Uses `obs_type` function to classify rows in data frame
- Splits data by type and sums likelihood contributions (assumes i.i.d.)
- Dynamically dispatches to type-specific functions:
  - Looks for `loglik.<type>`, `score.<type>`, `hess_loglik.<type>` functions
  - Falls back to numerical differentiation if analytical methods not found
- Stores dispatchers in lists: `$logliks`, `$scores`, `$hess_logliks`
- **Important**: R6 methods are called with `$` (e.g., `model$loglik(df, par)`) not S3 generics
- S3 wrapper methods exist: `loglik.likelihood_contr_model()` wraps `model$loglik()`

**`likelihood_name_model`** (in `R/likelihood_name.R`):
- Wraps standard R distributions (e.g., "norm", "weibull", "exp")
- Uses R naming convention: `d<name>` (pdf), `p<name>` (cdf)
- Handles censoring via `censor_col`: "left", "right", "interval", or "exact"/NA
- `prepare_args()` utility matches parameters to distribution function arguments

**`likelihood_exact_weibull`** (in `R/likelihood_exact_weibull.R`):
- Example of a specialized model with analytical derivatives
- Demonstrates providing explicit `loglik`, `score`, and `hess_loglik` methods
- Can be used standalone or as a contribution in `likelihood_contr_model`

### Key Implementation Details

- **Method dispatch**: Uses S3 generics. Models add `"likelihood_model"` to class vector to inherit default methods
- **Numerical defaults**: `score.likelihood_model()` and `hess_loglik.likelihood_model()` use Richardson extrapolation with `r=6` via `numDeriv`
- **Optimization**: `fit.likelihood_model()` supports multiple `optim` methods (Nelder-Mead, BFGS, SANN, etc.). SANN useful for multimodal likelihoods but should be followed by local search
- **Bootstrap**: `sampler.likelihood_model()` uses `boot::boot()` with multicore parallelization support

### File Structure

- `R/likelihood_model.R`: Core generic framework and default implementations
- `R/likelihood_contr_model.R`: R6 class for contribution-based models
- `R/likelihood_name.R`: Named distribution wrapper
- `R/likelihood_exact_weibull.R`: Example specialized model with analytical derivatives
- `R/utils.R`: Helper functions (`prepare_args`, `get_params`, `method_exists`)
- `R/tests.R`: Likelihood ratio test (`lrt()`)
- `tests/test.R`: Basic unit tests (**incomplete** - has placeholder tests with `...`)
- `vignettes/likelihood-contributions.Rmd`: Detailed examples with series systems
- `vignettes/likelihood-name-model.Rmd`: Examples using named distributions

## Important Notes

### Testing and Coverage
- **Tests are incomplete**: `tests/test.R` contains placeholder tests with `...` that need implementation
- When adding new features, write complete tests and run `covr::package_coverage()` to ensure coverage
- Focus on testing edge cases: censored data, empty data frames, parameter validation
- See vignettes for realistic examples to base tests on

### Adding New Distributions
When adding new distribution support, choose one of these approaches:
1. **For standard R distributions**: Use `likelihood_name("dist_name", ob_col, censor_col)`
   - Works with any distribution that has `d<name>` and `p<name>` functions
   - Automatically handles censoring types
2. **For custom likelihood contributions**: Create `loglik.<type>`, `score.<type>`, `hess_loglik.<type>` functions
   - Used by `likelihood_contr_model` via dynamic dispatch
   - Provide analytical derivatives when possible (much faster than numerical)
3. **For specialized models**: Create a new class (like `likelihood_exact_weibull`)
   - Define S3 methods: `loglik.<class>`, `score.<class>`, `hess_loglik.<class>`
   - Add `"likelihood_model"` to class vector to inherit default behavior

### Performance Considerations
- **Numerical differentiation is expensive**: Default `score` and `hess_loglik` use `numDeriv` with Richardson extrapolation (`r=6`)
- Provide analytical derivatives whenever possible (see `likelihood_exact_weibull.R` for example)
- For Weibull: analytical derivatives ~10-100x faster than numerical
- Parameters can be named or unnamed; `prepare_args()` handles both cases
- Always validate inputs: check for NULL, positive values where required, sufficient sample size

### Package Dependencies
- **Required**: `mvtnorm`, `generics`, `MASS`, `stats`, `numDeriv`, `boot`, `dplyr`
- **Companion packages**: `algebraic.mle` (MLE objects), `algebraic.dist` (distributions)
- **Roxygen2**: Version 7.3.1 for documentation generation

## Common Workflows

### Using a Standard Distribution (Simple Case)
```r
# Create model for normal distribution
model <- likelihood_name("norm", ob_col = "x", censor_col = NULL)

# Fit the model
solver <- fit(model)
result <- solver(df, par = c(0, 1))  # initial guess: mean=0, sd=1

# Get parameter estimates
params(result)
```

### Handling Censored Data
```r
# Data with censoring column ("exact", "left", "right", "interval")
df <- data.frame(
  x = c(1.2, 3.4, 2.1),
  censor = c("exact", "right", "left")
)

model <- likelihood_name("weibull", ob_col = "x", censor_col = "censor")
solver <- fit(model)
result <- solver(df, par = c(shape = 2, scale = 1))
```

### Using Likelihood Contributions (Advanced Case)
```r
# Define observation type function
obs_type <- function(df) {
  # Returns "exact", "censored", etc. for each row
  df$type
}

# Define type-specific log-likelihood functions
loglik.exact <- function(df, par) {
  # Implementation for exact observations
}

loglik.censored <- function(df, par) {
  # Implementation for censored observations
}

# Create contribution model
model <- likelihood_contr_model$new(
  obs_type = obs_type,
  assumptions = c("iid", "exponential components")
)

# Fit as usual
solver <- fit(model)
result <- solver(df, par = initial_params)
```

### Likelihood Ratio Test
```r
# Nested models (null is subset of alternative)
null_result <- fit(null_model)(df, par_null)
alt_result <- fit(alt_model)(df, par_alt)

# Perform LRT
lrt_result <- lrt(
  null = null_model,
  alt = alt_model,
  data = df,
  null_par = params(null_result),
  alt_par = params(alt_result),
  dof = length(params(alt_result)) - length(params(null_result))
)

# Check p-value
lrt_result$p.value
```
