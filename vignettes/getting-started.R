## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 6,
    fig.height = 4
)
options(digits = 4)

## ----install, eval=FALSE------------------------------------------------------
# # Install from GitHub
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("queelius/likelihood.model")

## ----load, message=FALSE, warning=FALSE---------------------------------------
library(likelihood.model)

## ----weibull-example----------------------------------------------------------
# Generate synthetic survival data
set.seed(42)
n <- 150
true_shape <- 2.5
true_scale <- 3.0
df <- data.frame(x = rweibull(n, shape = true_shape, scale = true_scale))

# Create the likelihood model
model <- weibull_uncensored("x")

# View model assumptions
assumptions(model)

# Fit the MLE
solver <- fit(model)
mle <- solver(df, par = c(1, 1))  # initial guess: shape=1, scale=1

# View results
summary(mle)

## ----weibull-results----------------------------------------------------------
# Parameter estimates
cat("Estimated parameters:\n")
print(coef(mle))

# Confidence intervals (Wald-based)
cat("\n95% Confidence Intervals:\n")
print(confint(mle))

# Standard errors
cat("\nStandard Errors:\n")
print(se(mle))

# Log-likelihood value
cat("\nLog-likelihood:", loglik_val(mle), "\n")

# AIC for model selection
cat("AIC:", aic(mle), "\n")

## ----censoring-example--------------------------------------------------------
# Generate normal data with right-censoring
set.seed(123)
n <- 200
true_mean <- 10
true_sd <- 2
censor_time <- 11

# Generate latent (true) values
x_latent <- rnorm(n, true_mean, true_sd)

# Apply censoring
censored <- x_latent > censor_time
df_cens <- data.frame(
  x = ifelse(censored, censor_time, x_latent),
  censor = ifelse(censored, "right", "exact")
)

cat("Censoring rate:", mean(censored) * 100, "%\n")

## ----censoring-fit------------------------------------------------------------
# Create model that handles censoring
model_cens <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

# Fit the model
solver_cens <- fit(model_cens)
mle_cens <- suppressWarnings(solver_cens(df_cens, par = c(0, 1)))

cat("MLE (accounting for censoring):\n")
cat("  Mean:", coef(mle_cens)[1], "(true:", true_mean, ")\n")
cat("  SD:  ", coef(mle_cens)[2], "(true:", true_sd, ")\n")

## ----censoring-comparison-----------------------------------------------------
# Naive estimates (ignoring censoring)
naive_mean <- mean(df_cens$x)
naive_sd <- sd(df_cens$x)

cat("\nComparison of estimates:\n")
cat("                     Mean      SD\n")
cat(sprintf("True values:      %7.3f   %5.3f\n", true_mean, true_sd))
cat(sprintf("MLE (correct):    %7.3f   %5.3f\n", coef(mle_cens)[1], coef(mle_cens)[2]))
cat(sprintf("Naive (biased):   %7.3f   %5.3f\n", naive_mean, naive_sd))
cat("\nNote: Naive estimates are biased downward due to censoring!\n")

## ----bootstrap-example, cache=TRUE--------------------------------------------
# Using the Weibull example from before
model <- weibull_uncensored("x")
df <- data.frame(x = rweibull(150, shape = 2.5, scale = 3.0))

# Create bootstrap sampler
boot_sampler <- sampler(model, df = df, par = c(1, 1))

# Generate bootstrap samples (100 replicates for speed)
boot_result <- boot_sampler(n = 100)

# Compare asymptotic vs bootstrap standard errors
mle <- fit(model)(df, par = c(1, 1))
asymp_se <- se(mle)
boot_se <- se(boot_result)  # Uses vcov from bootstrap replicates

cat("Standard Error Comparison:\n")
cat("           Asymptotic  Bootstrap\n")
cat(sprintf("Shape:     %10.4f  %9.4f\n", asymp_se[1], boot_se[1]))
cat(sprintf("Scale:     %10.4f  %9.4f\n", asymp_se[2], boot_se[2]))

# Bootstrap bias estimate
cat("\nBootstrap Bias Estimate:\n")
print(bias(boot_result))

## ----score-check--------------------------------------------------------------
model <- weibull_uncensored("x")
df <- data.frame(x = rweibull(100, 2, 1.5))

# Fit MLE
mle <- fit(model)(df, par = c(1, 1))

# Evaluate score at MLE
score_func <- score(model)
score_at_mle <- score_func(df, coef(mle))

cat("Score at MLE (should be near zero):\n")
print(score_at_mle)

# The score is also stored in the MLE object
cat("\nScore stored in MLE object:\n")
print(score_val(mle))

## ----fisherian-example--------------------------------------------------------
# Fit a model
model <- weibull_uncensored("x")
df <- data.frame(x = rweibull(100, 2.0, 1.5))
mle <- fit(model)(df, par = c(1, 1))

# Support function: log relative likelihood
# S(theta) = logL(theta) - logL(theta_hat)
# Support at MLE is always 0
s_at_mle <- support(mle, coef(mle), df, model)
cat("Support at MLE:", s_at_mle, "\n")

# Support at alternative parameter values (negative = less supported)
s_alt <- support(mle, c(1.5, 1.0), df, model)
cat("Support at (1.5, 1.0):", s_alt, "\n")

# Relative likelihood = exp(support)
rl <- relative_likelihood(mle, c(1.5, 1.0), df, model)
cat("Relative likelihood at (1.5, 1.0):", rl, "\n")

## ----likelihood-intervals-----------------------------------------------------
# Compute 1/8 likelihood interval (roughly equivalent to 95% CI)
# This is the set of theta where R(theta) >= 1/8
li <- likelihood_interval(mle, df, model, k = 8)
print(li)

# Compare with Wald confidence interval
cat("\nWald 95% CI:\n")
print(confint(mle))

## ----exponential-example------------------------------------------------------
# Generate exponential survival data
set.seed(99)
df_exp <- data.frame(t = rexp(200, rate = 3.0))

# Create model and fit -- no initial guess needed!
model_exp <- exponential_lifetime("t")
mle_exp <- fit(model_exp)(df_exp)

cat("Closed-form MLE (no optim):\n")
cat("  lambda_hat:", coef(mle_exp), "(true: 3.0)\n")
cat("  SE:", se(mle_exp), "\n")

# The score at the MLE is exactly zero (by construction)
cat("  Score at MLE:", score_val(mle_exp), "\n")

## ----exponential-censored-----------------------------------------------------
# Generate data with right-censoring at time 0.5
set.seed(99)
true_lambda <- 3.0
x <- rexp(200, rate = true_lambda)
censored <- x > 0.5
df_cens_exp <- data.frame(
  t = pmin(x, 0.5),
  status = ifelse(censored, "right", "exact")
)

cat("Censoring rate:", mean(censored) * 100, "%\n")

model_cens_exp <- exponential_lifetime("t", censor_col = "status")
mle_cens_exp <- fit(model_cens_exp)(df_cens_exp)

cat("MLE (censored):", coef(mle_cens_exp), "(true:", true_lambda, ")\n")
cat("95% CI:", confint(mle_cens_exp)[1, ], "\n")

## ----exponential-crossval-----------------------------------------------------
# Compare log-likelihood values at the MLE
lambda_hat <- coef(mle_exp)
ll_analytical <- loglik(model_exp)(df_exp, lambda_hat)

# likelihood_name("exp") uses dexp(x, rate), so pass unnamed parameter
df_exp_name <- data.frame(t = df_exp$t, censor = rep("exact", nrow(df_exp)))
model_generic <- likelihood_name("exp", ob_col = "t", censor_col = "censor")
ll_generic <- loglik(model_generic)(df_exp_name, unname(lambda_hat))

cat("Analytical loglik:", ll_analytical, "\n")
cat("Generic loglik:   ", ll_generic, "\n")
cat("Match:", isTRUE(all.equal(unname(ll_analytical), ll_generic, tolerance = 1e-6)), "\n")

## ----session------------------------------------------------------------------
sessionInfo()

