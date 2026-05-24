## R CMD check results

0 errors | 0 warnings | 1 note

* One NOTE: "unable to verify current time" (system clock check, benign).

## Patch (v1.0.0 -> v1.0.1)

Bug fix in `fim.likelihood_model`. The default Monte Carlo FIM
previously forwarded the inner closure's `...` to both `rdata` and
`hess_loglik`. The docstring already documented the correct contract
("Additional arguments passed to rdata"); the implementation now
matches.

Two motivations:

1. The likelihood is a function of data alone. Kwargs that describe
   the data-generating process (censoring time, masking probability,
   observation functor, etc.) belong to the sampling layer and must
   not appear in the likelihood evaluator. The previous behavior
   violated the standard sampling/inference separation.

2. R partial-argument-matching turned the leak into a deterministic
   crash: a DGP kwarg whose name shared a prefix with a `hess_loglik`
   formal (e.g., `p` colliding with the universal `par` formal) was
   silently misrouted at the call site, corrupting the parameter
   vector. Downstream packages with masked-data models could not call
   `fim()` with `p`-style kwargs at all.

The existing test for the leaky behavior was rewritten in two parts:
one pinning that `...` does NOT reach `hess_loglik`, the other
pinning that `...` does reach `rdata`. No public API changes.

## Test environments

* local Ubuntu 24.04, R 4.3.3
* win-builder (R-devel and R-release)
