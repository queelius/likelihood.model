# likelihood.model 1.0.0

## Breaking Changes
* Removed `likelihood_contr_model` (R6 contribution-based model builder). This functionality has been extracted to the new `likelihood.contr` package.
* Removed `likelihood_name()` (named distribution wrapper). Use `likelihood.contr::contr_name()` instead.
* Removed `weibull_uncensored` / `likelihood_exact_weibull` (Weibull example model). Use `likelihood.contr` for Weibull likelihood contributions.
* Removed `observed_info()` generic (was deprecated in v0.10.0). Use `-hess_loglik(model)(df, par)` or `observed_fim()` on MLE results.
* Dropped `R6` from Imports (no longer needed).
* `mvtnorm` moved to Suggests (only needed for multivariate `sampler()`).
* Dropped `dplyr` and `tibble` from Suggests.

## Notes
* The core framework (`loglik`, `score`, `hess_loglik`, `fim`, `fit`, Fisherian inference, `lrt`, `fisher_mle`, `fisher_boot`) is unchanged.
* `exponential_lifetime` remains as the reference implementation.
* For contribution-based models with heterogeneous observation types, use the `likelihood.contr` package.

# likelihood.model 0.10.0

## Breaking Changes
* `observed_info()` is deprecated; use `-hess_loglik(model)(df, par)` or `observed_fim()` on MLE results
* `likelihood_exact_weibull()` is deprecated in favor of `weibull_uncensored()`
* Unrecognized censoring values in `loglik.likelihood_name_model()` now error instead of warning (silently returning 0 produced incorrect MLEs during optimization)
* Removed redundant R6 instance methods `$loglik()`, `$score()`, `$hess_loglik()` from `likelihood_contr_model`; use the S3 generics `loglik(model)`, `score(model)`, `hess_loglik(model)` instead
* `mvtnorm` moved from Imports to Suggests (only needed for multivariate `sampler()`)

## Bug Fixes
* Fixed `fim.likelihood_model()` Monte Carlo estimator: replaced biased outer-product-of-scores with negative expected Hessian approach, which is unbiased at any parameter value
* Fixed `fim.likelihood_model()` not forwarding `...` to `hess_loglik`, and added NaN/Inf guard to prevent silent corruption of FIM estimates
* Fixed `...` not forwarded in `likelihood_contr_model` numerical fallback for `score` and `hess_loglik`
* Fixed `parent.frame()` in `likelihood_contr_model$new()` capturing R6 internals instead of the caller's environment
* Fixed `sampler.likelihood_model()` `...` scoping bug in bootstrap statistic function
* Added optimizer-friendly boundary guards to `weibull_uncensored`: `loglik` returns penalty, `score`/`hess_loglik` return zeros for invalid parameters (k <= 0 or l <= 0)
* Added input validation in `loglik.likelihood_name_model()`: errors with a clear message when `ob_col` is missing from the data frame

## New Features
* New `rdata.weibull_uncensored()` method for generating Weibull data from the model's DGP
* `bias.fisher_mle()` now supports optional Monte Carlo bias estimation when `theta` and `model` are provided; returns zeros (asymptotic result) otherwise
* `mse.fisher_mle()` now forwards `model` and `n_sim` to `bias()` for MC-based MSE
* `likelihood_contr_model$new()` gains `envir` parameter controlling where `loglik.<type>`, `score.<type>`, and `hess_loglik.<type>` are looked up
* `fim.likelihood_model()` now adds parameter names to the output matrix when `theta` is named, and includes `tryCatch` with informative warnings on partial failures
* `evidence()` is now an S3 generic, extensible by model classes
* `lrt()` returns a proper `lrt_result` class with a `print` method
* New `loglik.likelihood_model` and `rdata.likelihood_model` default methods provide clear error messages when a model class doesn't implement the required concept methods

# likelihood.model 0.9.1

## New Features
* Added interval censoring support to `likelihood_name()` via `ob_col_upper` parameter
* Input validation in `likelihood_name()` constructor: checks that distribution functions exist
* Parameter count validation in `prepare_args_list()` with clear error messages
* Warning on unrecognized censoring values (e.g., typos like "righ")
* Auto-switch from Nelder-Mead to BFGS for 1-dimensional optimization
* New Fisherian inference functions: `support()`, `relative_likelihood()`, `likelihood_interval()`, `profile_loglik()`, `evidence()`
* New `exponential_lifetime` model with closed-form MLE and right-censoring support
* New `fisher_mle` and `fisher_boot` result objects with full S3 method suite

## Bug Fixes
* Fixed `rowSums` error in `likelihood_contr_model` score for single-parameter models
* Fixed `print.likelihood_model` crash when `ob_col` is NULL (e.g., for contr models)

## algebraic.mle Integration
* New re-exported generics from `algebraic.mle`: `params()`, `nparams()`, `observed_fim()`, `obs()`, `mse()`
* Methods for `fisher_mle` objects provide a unified interface across the likelihood ecosystem

## Documentation
* New `algebraic-mle-integration` vignette demonstrating the three-package ecosystem (`likelihood.model` + `algebraic.mle` + `algebraic.dist`)
* Expanded `likelihood-name-model` vignette: left-censoring, interval censoring, mixed types, Fisherian inference, multiple distributions
* Expanded `likelihood-contributions` vignette: Case 2 inference, single-parameter example, best practices, interpretive prose for computed results
* New `getting-started` vignette with comprehensive introduction
* New `exponential-lifetime` vignette

# likelihood.model 0.9.0

* Initial public release with `likelihood_name()`, `likelihood_contr_model`, and `weibull_uncensored`
