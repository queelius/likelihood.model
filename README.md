# likelihood.model

Likelihood-based statistical inference in the Fisherian tradition for R.

## Why likelihood?

Classical statistics often jumps from data to p-values or posterior distributions. The likelihood approach takes a different path: the likelihood function *is* the complete summary of what the data say about the parameters. No priors, no long-run frequency guarantees -- just the evidence in the data.

This package makes it easy to:

- **Specify** likelihood models for standard distributions or custom observation types
- **Estimate** parameters via maximum likelihood, with analytical or numerical derivatives
- **Quantify evidence** using the likelihood function directly -- support functions, relative likelihood, and likelihood intervals
- **Handle censored data** (left, right, interval) as a natural part of the likelihood framework

## Installation

```r
# From CRAN
install.packages("likelihood.model")

# Development version from GitHub
devtools::install_github("queelius/likelihood.model")
```

## Quick start

### Fit a distribution

```r
library(likelihood.model)

# Any R distribution works: "norm", "weibull", "exp", "gamma", ...
model <- likelihood_name("weibull", ob_col = "time", censor_col = "status")
mle <- fit(model)(df, par = c(shape = 1, scale = 1))

coef(mle)      # parameter estimates
se(mle)        # standard errors
confint(mle)   # confidence intervals
```

### Handle censored data properly

Ignoring censoring biases your estimates. `likelihood.model` handles it correctly:

```r
df <- data.frame(
  time   = c(1.2, 3.4, 5.0, 2.1, 5.0),
  status = c("exact", "exact", "right", "exact", "right")
)

model <- likelihood_name("weibull", ob_col = "time", censor_col = "status")
mle <- fit(model)(df, par = c(shape = 1, scale = 1))
```

Right-censored observations contribute through the survival function rather than the density -- the likelihood framework handles this seamlessly.

### Fisherian inference

Instead of confidence intervals that rely on repeated-sampling arguments, likelihood intervals identify parameter values well-supported by the data:

```r
# Likelihood interval: {theta : L(theta)/L(theta_hat) >= 1/8}
li <- likelihood_interval(mle, df, model, k = 8)

# How well does the data support theta1 over theta2?
evidence(model, df, theta1 = c(2, 3), theta2 = c(1, 5))
```

## Three ways to build a model

| Approach | Use when | Example |
|----------|----------|---------|
| `likelihood_name()` | Standard R distribution | `likelihood_name("norm", "x")` |
| Specialized model | You need analytical derivatives or closed-form MLE | `weibull_uncensored("x")`, `exponential_lifetime("t")` |
| `likelihood_contr_model` | Heterogeneous observation types (mixed exact/censored) | Custom likelihood contributions |

### `likelihood_name()` -- the quick path

Wraps any R distribution that has `d<name>` and `p<name>` functions. Supports exact, left-censored, right-censored, and interval-censored observations automatically.

```r
model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
```

### Specialized models -- the fast path

When you have analytical derivatives, estimation is 10-100x faster than numerical differentiation:

```r
# Weibull with hand-derived score and Hessian
model <- weibull_uncensored("x")

# Exponential with closed-form MLE (bypasses optim entirely)
model <- exponential_lifetime("t", censor_col = "status")
mle <- fit(model)(df)  # no initial guess needed
```

### `likelihood_contr_model` -- the flexible path

For complex models with different observation types, define type-specific log-likelihood contributions:

```r
model <- likelihood_contr_model$new(
  obs_type = function(df) df$type,
  logliks = list(exact = loglik_exact, censored = loglik_censored)
)
```

## The likelihood calculus

Every likelihood model provides a consistent set of operations:

```r
loglik(model)       # log-likelihood function
score(model)        # score function (gradient)
hess_loglik(model)  # Hessian matrix
fim(model)          # Fisher information matrix
```

If you provide only `loglik()`, the package computes `score()` and `hess_loglik()` via numerical differentiation (Richardson extrapolation). Provide analytical versions for better speed and accuracy.

## Inference toolkit

### Maximum likelihood estimation

```r
mle <- fit(model)(df, par = c(1, 1))
coef(mle)           # point estimates
vcov(mle)           # variance-covariance matrix
confint(mle)        # Wald confidence intervals
summary(mle)        # full summary
```

### Bootstrap

```r
boot <- sampler(model, df = df, par = c(1, 1))(n = 1000)
se(boot)            # bootstrap standard errors
confint(boot)       # bootstrap confidence intervals (BCa)
bias(boot)          # bootstrap bias estimate
```

### Model comparison

```r
lrt(null_model, alt_model, df,            # likelihood ratio test
    null_par = p0, alt_par = p1, dof = 1)
aic(mle)                                   # Akaike information criterion
bic(mle)                                   # Bayesian information criterion
```

### Fisherian likelihood inference

```r
support(mle, theta, df, model)             # S(theta) = logL(theta) - logL(theta_hat)
relative_likelihood(mle, theta, df, model) # R(theta) = L(theta) / L(theta_hat)
likelihood_interval(mle, df, model, k = 8) # {theta : R(theta) >= 1/k}
profile_loglik(mle, df, model, param = 1)  # profile out nuisance parameters
evidence(model, df, theta1, theta2)        # log-likelihood ratio
```

## Ecosystem

`likelihood.model` interoperates with two companion packages:

- **[algebraic.mle](https://github.com/queelius/algebraic.mle)** -- MLE result objects with a unified interface (`params`, `rmap`, `marginal`)
- **[algebraic.dist](https://github.com/queelius/algebraic.dist)** -- probability distribution objects for simulation and analysis

## Learn more

- [Getting started](https://queelius.github.io/likelihood.model/articles/getting-started.html) -- worked examples from basic MLE to Fisherian inference
- [Named distributions](https://queelius.github.io/likelihood.model/articles/likelihood-name-model.html) -- censoring, interval data, model comparison
- [Likelihood contributions](https://queelius.github.io/likelihood.model/articles/likelihood-contributions.html) -- custom models for series systems and mixed observation types
- [Exponential lifetime](https://queelius.github.io/likelihood.model/articles/exponential-lifetime.html) -- closed-form MLE, analytical derivatives, Monte Carlo validation
- [algebraic.mle integration](https://queelius.github.io/likelihood.model/articles/algebraic-mle-integration.html) -- the three-package ecosystem
- [Reference](https://queelius.github.io/likelihood.model/reference/) -- full API documentation
