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

## Documentation
* Expanded `likelihood-name-model` vignette: left-censoring, interval censoring, mixed types, Fisherian inference, multiple distributions
* Expanded `likelihood-contributions` vignette: Case 2 inference, single-parameter example, best practices
* New `getting-started` vignette with comprehensive introduction
* New `exponential-lifetime` vignette

# likelihood.model 0.9.0

* Initial public release with `likelihood_name()`, `likelihood_contr_model`, and `weibull_uncensored`
