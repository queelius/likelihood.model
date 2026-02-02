## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 6,
    fig.height = 4
)
options(digits = 4)

## ----load, message=FALSE, warning=FALSE---------------------------------------
library(likelihood.model)
library(algebraic.mle)

## ----basic-fit----------------------------------------------------------------
set.seed(42)
true_shape <- 2.5
true_scale <- 3.0
n <- 200
df <- data.frame(x = rweibull(n, shape = true_shape, scale = true_scale))

model <- weibull_uncensored("x")
mle_result <- fit(model)(df, par = c(1, 1))

# algebraic.mle generics — these work because fisher_mle inherits from mle
params(mle_result)
nparams(mle_result)

## ----observed-fim-------------------------------------------------------------
observed_fim(mle_result)

## ----fim-pd-------------------------------------------------------------------
eigen(observed_fim(mle_result))$values

## ----vcov-vs-fim--------------------------------------------------------------
vcov(mle_result)
solve(observed_fim(mle_result))

## ----rmap---------------------------------------------------------------------
# Transform Weibull (shape, scale) -> mean lifetime
mean_life_mle <- rmap(mle_result, function(p) {
  c(mean_lifetime = p[2] * gamma(1 + 1 / p[1]))
})

params(mean_life_mle)
se(mean_life_mle)

## ----rmap-compare-------------------------------------------------------------
true_mean <- true_scale * gamma(1 + 1 / true_shape)
cat("True mean lifetime:", true_mean, "\n")
cat("Estimated mean lifetime:", params(mean_life_mle), "\n")
cat("95% CI:", confint(mean_life_mle), "\n")

## ----rmap-multi---------------------------------------------------------------
# Derive mean, variance, and median of the Weibull distribution
derived_mle <- rmap(mle_result, function(p) {
  k <- p[1]; lam <- p[2]
  c(
    mean   = lam * gamma(1 + 1/k),
    var    = lam^2 * (gamma(1 + 2/k) - gamma(1 + 1/k)^2),
    median = lam * log(2)^(1/k)
  )
})

params(derived_mle)
se(derived_mle)

## ----marginal-----------------------------------------------------------------
# Marginal for shape parameter
shape_mle <- marginal(mle_result, 1)
params(shape_mle)
se(shape_mle)
confint(shape_mle)

# Marginal for scale parameter
scale_mle <- marginal(mle_result, 2)
params(scale_mle)
se(scale_mle)
confint(scale_mle)

## ----expectation--------------------------------------------------------------
set.seed(123)
# E[shape^2] under the asymptotic distribution
e_shape_sq <- expectation(mle_result, function(p) p[1]^2,
                          control = list(n = 10000L))
cat("E[shape^2]:", e_shape_sq, "\n")
cat("shape^2 at MLE:", params(mle_result)[1]^2, "\n")

# Probability that shape > 2 (under asymptotic distribution)
pr_shape_gt_2 <- expectation(mle_result, function(p) as.numeric(p[1] > 2),
                             control = list(n = 10000L))
cat("P(shape > 2):", pr_shape_gt_2, "\n")

## ----mse----------------------------------------------------------------------
mse(mle_result)
all.equal(mse(mle_result), vcov(mle_result))

## ----bootstrap, cache=TRUE----------------------------------------------------
set.seed(42)
boot_sampler <- sampler(model, df = df, par = c(1, 1))
boot_result <- boot_sampler(n = 200)

# Same algebraic.mle generics work
params(boot_result)
nparams(boot_result)
se(boot_result)
bias(boot_result)

## ----boot-compare-------------------------------------------------------------
cat("Asymptotic 95% CI:\n")
confint(mle_result)

cat("\nBootstrap percentile 95% CI:\n")
confint(boot_result, type = "perc")

## ----dist-check, include=FALSE------------------------------------------------
has_dist <- requireNamespace("algebraic.dist", quietly = TRUE)

## ----dist-section, eval=has_dist----------------------------------------------
library(algebraic.dist)

## ----dist-compare, eval=has_dist----------------------------------------------
set.seed(42)
lambda_true <- 2.0
n_exp <- 200
df_exp <- data.frame(t = rexp(n_exp, rate = lambda_true))

model_exp <- exponential_lifetime("t")
mle_exp <- fit(model_exp)(df_exp)

cat("MLE:", params(mle_exp), "\n")
cat("SE:", se(mle_exp), "\n")

# Theoretical asymptotic distribution of the MLE
asymp_var <- lambda_true^2 / n_exp
asymp_dist <- normal(mu = lambda_true, var = asymp_var)
cat("\nTheoretical asymptotic distribution:\n")
cat("  Mean:", params(asymp_dist)[1], "\n")
cat("  Variance:", params(asymp_dist)[2], "\n")

# Compare: sample from the MLE's estimated distribution
mle_sampler <- sampler(mle_exp)
set.seed(1)
mle_samples <- mle_sampler(5000)

# vs. sample from the theoretical distribution
dist_sampler <- sampler(asymp_dist)
set.seed(1)
dist_samples <- dist_sampler(5000)

cat("\nMLE sampler mean:", mean(mle_samples), "\n")
cat("Theoretical sampler mean:", mean(dist_samples), "\n")
cat("MLE sampler sd:", sd(mle_samples), "\n")
cat("Theoretical sampler sd:", sd(dist_samples), "\n")

