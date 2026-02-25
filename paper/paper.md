---
title: 'likelihood.model: A Generic Framework for Likelihood-Based Inference in R'
tags:
  - R
  - likelihood inference
  - maximum likelihood estimation
  - censored data
  - survival analysis
authors:
  - name: Alexander Towell
    orcid: 0000-0001-6443-9897
    affiliation: 1
affiliations:
  - name: Department of Computer Science, Southern Illinois University Edwardsville
    index: 1
date: 24 February 2026
doi: 10.5281/zenodo.18463604
bibliography: paper.bib
---

# Summary

`likelihood.model` is an R [@R] package that provides a generic framework
for specifying likelihood models and performing inference in the Fisherian
tradition following @Royall1997. The package defines a `likelihood_model`
concept through S3 generics---`loglik()`, `score()`, `hess_loglik()`, and
`fim()`---that any model class can implement. When only a log-likelihood
is provided, the score and Hessian are computed automatically via
Richardson extrapolation [@numDeriv]. Users who supply analytical
derivatives gain 10--100x speedups without changing the inference
interface.

The package offers three model-building strategies at different levels of
generality. `likelihood_name()` wraps any standard R distribution (normal,
Weibull, exponential, gamma, etc.) into a likelihood model with automatic
support for exact, left-censored, right-censored, and interval-censored
observations. `likelihood_contr_model`, an R6 class, handles heterogeneous
observation types by dispatching to type-specific log-likelihood
contributions and summing them under an independence assumption---a
pattern common in reliability engineering and survival analysis where a
single data set mixes exact failure times with various forms of censoring.
Specialized model classes such as `weibull_uncensored` and
`exponential_lifetime` demonstrate how to provide analytical derivatives
or closed-form MLEs that bypass numerical optimization entirely.

All models produce `fisher_mle` result objects with a full S3 method
suite (`coef`, `vcov`, `confint`, `se`, `aic`, `bic`, `summary`) and
Fisherian inference functions: `support()`, `relative_likelihood()`,
`likelihood_interval()`, `profile_loglik()`, and `evidence()`. Bootstrap
sampling distributions and likelihood ratio tests for nested models are
also supported.

# Statement of Need

Researchers in survival analysis, reliability engineering, and
applied statistics who work with censored or heterogeneous data need a
framework that separates *model specification* from *inference machinery*.
Existing R tools for maximum likelihood estimation either focus on
optimization mechanics (requiring users to manage log-likelihoods,
gradients, and Hessians manually) or are tightly coupled to specific model
families. `likelihood.model` fills this gap by defining a minimal
interface---implement `loglik()` and inherit everything else---that lets
researchers focus on the statistical model while the package handles
derivative computation, optimization, bootstrap resampling, and Fisherian
inference automatically.

The package is designed as infrastructure for an ecosystem of
domain-specific packages. It sits between `algebraic.mle` [@algebraic.mle],
which provides MLE result algebra (delta method transformations, marginal
distributions), and downstream packages for reliability analysis such as
`flexhaz`, `maskedcauses`, and `maskedhaz`, which use the
`likelihood_contr_model` class to implement likelihood contributions for
series systems with masked failure causes.

# State of the Field

Several R packages address maximum likelihood estimation, but each
occupies a different niche. The `stats4::mle` function and its extension
`bbmle` [@bbmle] provide formula-based MLE fitting with profile likelihood
and model comparison, but require users to write the negative
log-likelihood as a single function and do not provide a generic interface
for score functions, Hessians, or Fisherian inference tools such as
likelihood intervals and evidence functions. `maxLik` [@maxLik] offers a
unified wrapper around multiple optimization algorithms (including
Newton-Raphson) with automatic gradient and Hessian computation, but its
focus is on the optimization step rather than on defining a reusable
likelihood model concept that downstream packages can extend.

For censored data specifically, `flexsurv` [@flexsurv] provides flexible
parametric survival models with user-definable hazard functions and
supports multiple censoring types, but it is oriented toward regression
modeling rather than providing a generic likelihood interface.
`fitdistrplus` [@fitdistrplus] fits univariate distributions to censored
data via several methods (MLE, moment matching, quantile matching) but
does not expose score functions, Fisher information, or Fisherian
inference tools.

`likelihood.model` differs from these packages in three ways. First, it
defines a composable *likelihood model concept* with S3 generics that
downstream packages implement, rather than a monolithic fitting function.
Second, its contribution-based architecture (`likelihood_contr_model`)
supports heterogeneous observation types---mixed exact, left-censored,
right-censored, and interval-censored data within a single model---with
dynamic dispatch to type-specific likelihood functions. Third, it provides
a complete suite of Fisherian inference tools (support functions, relative
likelihood, likelihood intervals, profile likelihood, evidence) as
first-class citizens rather than post-hoc additions.

# Software Design

The package is organized in three layers. The *core layer* defines the
`likelihood_model` concept as a set of S3 generics (`loglik`, `score`,
`hess_loglik`, `fim`) with default methods that fall back to numerical
differentiation via `numDeriv` [@numDeriv]. Any object with
`"likelihood_model"` in its class vector inherits a complete inference
toolkit: `fit()` for MLE via `optim`, `sampler()` for bootstrap sampling
distributions via `boot` [@boot], `lrt()` for likelihood ratio tests, and
the Fisherian inference functions. The `fisher_mle` result class inherits
from `algebraic.mle::mle`, ensuring compatibility with the `algebraic.mle`
interface for delta method transformations (`rmap`) and marginal
distributions (`marginal`).

The *model builder layer* provides two ready-to-use implementations.
`likelihood_name()` constructs a model from any R distribution name by
convention---if `dnorm` and `pnorm` exist, `likelihood_name("norm", ...)`
produces a complete likelihood model with censoring support.
`likelihood_contr_model` uses an R6 class with a closure-based caching
strategy: the S3 wrappers (`loglik(model)`, `score(model)`) return
functions that cache the data frame split and resolved dispatchers,
eliminating repeated overhead during the hundreds of function evaluations
that `optim` performs.

The *example layer* provides reference implementations---`weibull_uncensored`
with analytical score and Hessian, and `exponential_lifetime` with a
closed-form MLE that bypasses optimization---demonstrating how to
specialize the generic interface for performance or analytic convenience.

# Research Impact Statement

`likelihood.model` serves as the inference backbone for an ecosystem of
R packages targeting reliability engineering with masked failure data.
The `flexhaz` package uses `likelihood_model` to define distributions
via arbitrary hazard functions; `maskedcauses` and `maskedhaz` use
`likelihood_contr_model` to implement MLE for series systems where the
component cause of failure may be masked or censored; and
`compositional.mle` composes MLE solvers built on `likelihood.model`
via sequential chaining and parallel racing. All packages in this
ecosystem are available on the author's
[r-universe](https://queelius.r-universe.dev) and are undergoing CRAN
submission. The package has been submitted to CRAN as version 0.9.1, has
322 unit tests with approximately 97% code coverage, and includes five
vignettes with worked examples spanning basic MLE, censored data,
likelihood contributions for series systems, Fisherian inference, and
ecosystem integration.

# AI Usage Disclosure

The `likelihood.model` package was designed and implemented by the author
without AI assistance. Claude (Anthropic) was used to help write unit
tests. The JOSS paper draft was prepared with AI assistance and reviewed
by the author for accuracy and completeness.

# Acknowledgements

The author thanks Richard Royall, whose book *Statistical Evidence: A
Likelihood Paradigm* provided the conceptual foundation for the Fisherian
inference functions in this package.

# References
