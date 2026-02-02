## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)
options(digits = 2)

## ----install-likelihood-model, eval = FALSE-----------------------------------
# if (!require(devtools)) {
#     install.packages("devtools")
# }
# devtools::install_github("queelius/likelihood.model")

## ----load-packages, message=FALSE, warning=FALSE, echo = FALSE----------------
library(dplyr)
library(likelihood.model)
library(generics)

## ----data-gen-exp-series------------------------------------------------------
n <- 150
rates <- c(1.1, 1.2, 1.3)
set.seed(1235)

df <- data.frame(
    t1 = rexp(n, rates[1]),
    t2 = rexp(n, rates[2]),
    t3 = rexp(n, rates[3])
)
df$t <- apply(df, 1, min)

# map each observation to the corresponding index (component) that was minimum
df$k <- apply(df, 1, which.min)

for (i in 1:n) {
    for (p in 1:3) {
        if (p == df$k[i]) next
        df[i, p] <- NA
    }
}
head(df)

## ----observed-likelihood, cache=TRUE------------------------------------------
score.observed <- function(df, rates, ...) {
  rep(-sum(df$t), length(rates)) + (df |> dplyr::count(k))$n / rates
}

loglik.observed <- function(df, rates, ...) {
  sum(log(rates[df$k])) - sum(rates) * sum(df$t)
}

model.observed <- likelihood_contr_model$new(
  obs_type = function(df) { "observed" }
)

## ----mle-optim-observed, cache=TRUE-------------------------------------------
optim(rates, fn = loglik(model.observed), df=df, control=list(fnscale=-1))

## ----mle-observed-likelihood-model, cache=TRUE--------------------------------
rates.hat <- fit(model.observed)(df, par = rates)
summary(rates.hat)

## ----hess-loglike-------------------------------------------------------------
hess_loglik.observed <- function(df, rates, ...) {
  p <- length(rates)
  H <- matrix(0, p, p)
  counts <- df |> dplyr::count(k)
  for (j in 1:p) {
    H[j, j] <- -counts$n[j] / rates[j]^2
  }
  H
}
# Observed information from MLE (negative Hessian)
cat("Observed information matrix (from MLE):\n")
print(-rates.hat$hessian)

# Computed directly via hess_loglik
cat("\nHessian at MLE:\n")
print(hess_loglik(model.observed)(df, coef(rates.hat)))

## ----bootstrap, cache=TRUE----------------------------------------------------
model.samp <- sampler(model.observed, df = df, par = rates)
boot_result <- model.samp(n = 500)

cat("Bootstrap MLE:\n")
print(coef(boot_result))

cat("\nBootstrap standard errors:\n")
print(se(boot_result))

cat("\nBootstrap bias:\n")
print(bias(boot_result))

cat("\nBootstrap covariance:\n")
print(vcov(boot_result))

## -----------------------------------------------------------------------------
# if df$t > .5, then the system was right-censored
df.censored <- df |> mutate(
    censored = t > .5,
    t = ifelse(censored, .5, t)
)

df.censored[df.censored[, "censored"], "k"] <- NA
df.censored[df.censored[, "censored"], paste0("t", 1:3)] <- rep(NA, 3)
head(df.censored)

## -----------------------------------------------------------------------------
loglik.censored <- function(df, rates) {
    -sum(rates) * sum(df$t)
}
score.censored <- function(df, rates) {
    rep(-sum(df$t), length(rates))
}
hess_loglik_censored <- function(df, rates) {
    p <- length(rates)
    matrix(0, p, p)
}

model.censored <- likelihood_contr_model$new(
    obs_type = function(df) {
        ifelse(df$censored, "censored", "observed")
    }
)

mle.censored <- fit(model.censored)(df.censored, par = rates)
summary(mle.censored)

## ----case2-inference----------------------------------------------------------
cat("Uncensored estimates:\n")
cat("  Estimates:", coef(rates.hat), "\n")
cat("  Std errors:", se(rates.hat), "\n")

cat("\nCensored estimates:\n")
cat("  Estimates:", coef(mle.censored), "\n")
cat("  Std errors:", se(mle.censored), "\n")

cat("\n95% Confidence Intervals (uncensored):\n")
print(confint(rates.hat))

cat("\n95% Confidence Intervals (censored):\n")
print(confint(mle.censored))

cat("\nTrue parameters:", rates, "\n")

## ----case2-fisherian----------------------------------------------------------
# Support at MLE (should be 0)
s_at_mle <- support(mle.censored, coef(mle.censored), df.censored,
                    model.censored)
cat("Support at MLE:", s_at_mle, "\n")

# Evidence for MLE vs true parameters
ev <- evidence(model.censored, df.censored, coef(mle.censored), rates)
cat("Evidence (MLE vs true):", ev, "\n")

# Likelihood interval for rate 1 (first component)
li <- likelihood_interval(mle.censored, df.censored, model.censored,
                          k = 8, param = 1)
cat("\n1/8 likelihood interval for lambda_1:\n")
print(li)

## ----single-param-------------------------------------------------------------
# Define contribution types for a single-rate exponential
loglik.exp_exact <- function(df, par, ...) {
  lambda <- par[1]
  sum(log(lambda) - lambda * df$t)
}

loglik.exp_right <- function(df, par, ...) {
  lambda <- par[1]
  -lambda * sum(df$t)
}

score.exp_exact <- function(df, par, ...) {
  lambda <- par[1]
  nrow(df) / lambda - sum(df$t)
}

score.exp_right <- function(df, par, ...) {
  lambda <- par[1]
  -sum(df$t)
}

# Simulate right-censored exponential data
set.seed(42)
lambda_true <- 2.0
n_single <- 200
censor_time_single <- 0.4

raw_lifetimes <- rexp(n_single, lambda_true)
df_single <- data.frame(
  t = pmin(raw_lifetimes, censor_time_single),
  type = ifelse(raw_lifetimes > censor_time_single, "exp_right", "exp_exact")
)

cat("Exact observations:", sum(df_single$type == "exp_exact"), "\n")
cat("Right-censored:", sum(df_single$type == "exp_right"), "\n")

model_single <- likelihood_contr_model$new(
  obs_type = function(df) df$type,
  assumptions = c("exponential lifetimes", "non-informative right censoring")
)

# Fit and verify
mle_single <- fit(model_single)(df_single, par = c(lambda = 1))
cat("\nMLE:", coef(mle_single), "(true:", lambda_true, ")\n")

# Closed-form MLE for comparison: d / T
d <- sum(df_single$type == "exp_exact")
total_T <- sum(df_single$t)
cat("Closed-form MLE:", d / total_T, "\n")
cat("95% CI:", confint(mle_single)[1, ], "\n")

## ----analytical-example, eval=FALSE-------------------------------------------
# # Slow (default): numerical differentiation of loglik.my_type
# loglik.my_type <- function(df, par, ...) { ... }
# 
# # Fast: provide analytical derivatives
# score.my_type <- function(df, par, ...) { ... }
# hess_loglik.my_type <- function(df, par, ...) { ... }

## ----cache-tip, eval=FALSE----------------------------------------------------
# # Good: cached S3 wrapper
# ll <- loglik(model)
# result <- optim(par, fn = ll, df = df, control = list(fnscale = -1))
# 
# # Also good: using fit() which does this internally
# result <- fit(model)(df, par = par)
# 
# # Slower: calling R6 method directly in a loop (no caching)
# # model$loglik(df, par) -- re-splits df every call

