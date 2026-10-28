## likelihood.model 1.0.2

This patch corrects profile likelihood intervals that inferred positivity of
nuisance parameters from the signs of fitted values. Optimizer and root-finder
scales now use parameter uncertainty, with regression controls for a correlated
quadratic and normal observations at scales 1, 1e-3 and 1e-6.

## Checks

Ubuntu 24.04.3 LTS, R 4.3.3, x86_64 Linux, 2026-10-28:
R CMD check --as-cran: 0 errors, 0 warnings, 2 environment notes.

* Current time could not be verified.
* HTML/manual math validation skipped because tidy and V8 are unavailable.

Examples, regression tests, the existing full test script, vignette rebuilding
and PDF manual passed. The final repeat check disabled remote incoming checks;
the preceding full remote incoming check completed successfully. R-devel
win-builder validation is requested separately.
