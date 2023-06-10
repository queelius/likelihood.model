Likelihood contribution model
================

# `likelihood.model`

`likelihood.model` is an R package that allows you to…

## Installation

You can install the development version of `likelihood.model` from
[GitHub](https://github.com/queelius/likelihood.model) with:

``` r
if (!require(devtools)) {
    install.packages("devtools")
}
devtools::install_github("queelius/likelihood.model")
```

Now, we load the requisite libraries:

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(algebraic.mle)
library(likelihood.model)
```

    ## 
    ## Attaching package: 'likelihood.model'

    ## The following objects are masked from 'package:algebraic.mle':
    ## 
    ##     fim, score

## Exponentially distributed system lifetime

Suppose we have a population of systems which have exponentially
distributed lifetimes, in particular, the `i`-th system has lifetime
`T[i] ~ EXP(1)`. We are interested in the failure rate `1` of the
systems.

Let’s draw a sample of size
![n=25](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D25
"n=25").

``` r
n <- 100
rate <- 1
set.seed(123)
library(tibble)

df <- tibble(t = rexp(n, rate))
head(df)
```

    ## # A tibble: 6 × 1
    ##        t
    ##    <dbl>
    ## 1 0.843 
    ## 2 0.577 
    ## 3 1.33  
    ## 4 0.0316
    ## 5 0.0562
    ## 6 0.317

Analytically, we know that the MLE of the population failure rate is
given by `rate.hat = mean(df$t)`, and that the sampling distribution of
the MLE is given by `rate.hat ~ Normal(rate, rate^2 / n)`.

``` r
(rate.hat <- 1/mean(df$t))
```

    ## [1] 0.9562801

### Likelihood model

Let’s derive the likelihood model for this problem anyway.

``` r
loglik.exp <- function(row, par) {
    dexp(row$t, par, log = TRUE)
}

score.exp <- function(row, par) {
    1/par - row$t
}

hess_loglik.exp <- function(row, par) {
    -1/par^2
}

model.exp <- likelihood_contr_model$new(
    obs_type = function(row) "exp"
)
```

The MLE may be found with:

``` r
exp.mle <- mle_numerical(
    optim(
        rate,
        fn = loglik(model.exp),
        gr = score(model.exp),
        df = df, hessian = TRUE, lower = 0, upper = 5,
        method = "Brent", control = list(fnscale = -1)))
summary(exp.mle)
```

We see that the two MLEs, the analytical one and the numerical one, are
the same.

### Right-censoring

Suppose that we only observe the systems in the sample for
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1
"1") unit of time. Any system that fails before time
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1
"1") is observed, but any system that fails after time
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1
"1") is not observed. We say that the system is *right-censored* because
we only know that it survived at least until time
![1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1
"1").

We could have also specified the likelihood model by explicitly
specifying the log-likelihood, score, and Hessian functions:

``` r
model.exp.rt_censored <- likelihood_contr_model$new(
    obs_type = function(row) {
        row$obs_type
    },
    logliks = list(
        exact = loglik.exp,
        right_censored = function(row, par) {
            pexp(1, par, lower.tail = FALSE, log.p = TRUE)
        }
    )
)

df.rt_censored <- df %>% mutate(
    obs_type = if_else(t <= 1, "exact", "right_censored"))

ll.exp.rt_censored <- loglik(model.exp.rt_censored)
exp.mle.rt_censored <- algebraic.mle::mle_numerical(
    optim(rate, fn = ll.exp.rt_censored, df = df.rt_censored,
    hessian = TRUE, lower = 0, upper = 5,
    method = "Brent", control = list(fnscale = -1)))
summary(exp.mle.rt_censored)
```

We can add more types of observations to the model, for example, we can
add left-censored observations, where we only know that the system
failed before time
![.1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;.1
".1"), interval censoring, where we only know that the system failed
between times
![.1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;.1
".1") and
![.9](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;.9
".9"), exact times when the system failed between
![.9](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;.9
".9") and
![1.25](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1.25
"1.25"), and right-censored observations, where we only know that the
system survived at least until time
![1.25](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1.25
"1.25").

``` r
model.exp.censored <- likelihood_contr_model$new(
    obs_type = function(row) {
        row$obs_type
    },
    logliks = list(
        exact = loglik.exp,
        right_censored = function(row, par) {
            pexp(1.25, par, lower.tail = FALSE, log.p = TRUE)
        },
        left_censored = function(row, par) {
            pexp(.1, par, log.p = TRUE)
        },
        interval_censored = function(row, par) {
            pexp(.9, par, log.p = TRUE) - pexp(.1, par, log.p = TRUE)
        }
    )
)

df.censored <- df %>% mutate(
    obs_type = case_when(
        t <= .1 ~ "left_censored",
        t > .1 & t <= .9 ~ "interval_censored",
        t > .9 & t <= 1.25 ~ "exact",
        t > 1.25 ~ "right_censored"
    )
)

ll.exp.censored <- loglik(model.exp.censored)
s.exp.censored <- score(model.exp.censored)
exp.mle.censored <- algebraic.mle::mle_numerical(
    optim(rate, fn = ll.exp.censored, df = df.censored,
    hessian = TRUE, method="L-BFGS-B", lower = 0, control = list(fnscale = -1, trace = 1, REPORT = 1)))
summary(exp.mle.censored)
s.exp.censored(df.censored, rate)
```
