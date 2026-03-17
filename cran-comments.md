## Resubmission

This is v0.10.0, an update to the previously submitted v0.9.1 (2026-02-07).

Key changes since v0.9.1:
- Fixed biased FIM Monte Carlo estimator (outer-product-of-scores replaced with negative expected Hessian)
- Added `rdata.weibull_uncensored()` for data generation
- Improved `bias.fisher_mle()` with optional MC estimation
- Added `envir` parameter to `likelihood_contr_model` for explicit lookup control
- Input validation and `...` forwarding fixes

## R CMD check results

0 errors | 0 warnings | 1 note

* CRAN incoming feasibility: "Possibly misspelled words: Fisherian, Royall"
  - Fisherian: standard adjective referring to R.A. Fisher's school of inference
  - Royall: surname of Richard Royall, author of the cited reference

## Test environments
* local: Ubuntu 24.04.3 LTS, R 4.3.3 (2024-02-29)

## Notes
* Package depends on `algebraic.mle` (available on CRAN) and suggests `algebraic.dist` (GitHub only, used in one vignette behind `requireNamespace` guard).

## Downstream dependencies
None.
