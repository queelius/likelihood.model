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

## Final constrained-domain correction

Both profiling APIs now share a centered objective, verify finite converged
nuisance fits, and use derivative-free fallback when finite differences cross
a genuine model domain. The expanded regression script passes 22 assertions;
current coverage is 98.88%. The final archive passed local --as-cran with
0 errors, 0 warnings and the same 2 environment notes.

Current R-devel (2026-10-26 r90766), --no-manual --as-cran: Status: OK.
https://github.com/queelius/likelihood.model/actions/runs/40645641026/job/122673698013
