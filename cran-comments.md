## Resubmission

This is a resubmission addressing feedback from Konstanze Lauseker (2026-02-03):

1. **DESCRIPTION**: Removed single quotes around non-package/software terms.
2. **DESCRIPTION**: Added `()` after all function names (loglik(), score(), etc.).
3. **DESCRIPTION**: Added reference in proper CRAN format: Royall (1997) <ISBN:978-0412044113>.
4. **Vignettes**: All five vignettes now save and restore `options()` via `old_opts <- options(...); ...; options(old_opts)`.
5. **Rd files**: Added missing `\value` tags to `hess_loglik.likelihood_model` and `obs.fisher_mle`.

Previous NOTEs from incoming pre-tests have also been addressed:
- Escaped braces in `likelihood_interval.Rd` LaTeX formula.
- Wrapped `<name>` in `\code{}` in `likelihood_name.Rd` to prevent HTML parsing.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments
* local: Ubuntu 24.04.3 LTS, R 4.3.3 (2024-02-29)

## Notes
* This is a first CRAN submission (resubmission after reviewer feedback).
* Package depends on `algebraic.mle` (available on CRAN) and suggests `algebraic.dist` (GitHub only, used in one vignette).

## Downstream dependencies
None.
