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

## ----uncensored---------------------------------------------------------------
set.seed(42)
true_lambda <- 2.5
df <- data.frame(t = rexp(200, rate = true_lambda))

model <- exponential_lifetime("t")
assumptions(model)

## ----fit-uncensored-----------------------------------------------------------
mle <- fit(model)(df)

cat("MLE:", coef(mle), "(true:", true_lambda, ")\n")
cat("SE:", se(mle), "\n")
cat("95% CI:", confint(mle)[1, ], "\n")
cat("Score at MLE:", score_val(mle), "(exactly zero by construction)\n")

## ----verify-mle---------------------------------------------------------------
n <- nrow(df)
total_time <- sum(df$t)
manual_mle <- n / total_time

cat("fit() result:   ", coef(mle), "\n")
cat("Manual n/T:     ", manual_mle, "\n")
cat("Match:", all.equal(unname(coef(mle)), manual_mle), "\n")

## ----censored-data------------------------------------------------------------
set.seed(42)
true_lambda <- 2.0
censor_time <- 0.5

# Generate latent failure times
x <- rexp(300, rate = true_lambda)
censored <- x > censor_time

df_cens <- data.frame(
  t = pmin(x, censor_time),
  status = ifelse(censored, "right", "exact")
)

cat("Sample size:", nrow(df_cens), "\n")
cat("Censoring rate:", round(mean(censored) * 100, 1), "%\n")
cat("Exact observations (d):", sum(!censored), "\n")
cat("Total time (T):", round(sum(df_cens$t), 2), "\n")

## ----fit-censored-------------------------------------------------------------
model_cens <- exponential_lifetime("t", censor_col = "status")
assumptions(model_cens)

mle_cens <- fit(model_cens)(df_cens)

cat("\nMLE:", coef(mle_cens), "(true:", true_lambda, ")\n")
cat("SE:", se(mle_cens), "\n")
cat("95% CI:", confint(mle_cens)[1, ], "\n")

## ----censoring-bias-----------------------------------------------------------
# Wrong: treat all observations as exact
model_wrong <- exponential_lifetime("t")
mle_wrong <- fit(model_wrong)(df_cens)

cat("Comparison:\n")
cat("  True lambda:          ", true_lambda, "\n")
cat("  MLE (correct):        ", coef(mle_cens), "\n")
cat("  MLE (ignoring censor):", coef(mle_wrong), "\n")
cat("\nIgnoring censoring overestimates the failure rate!\n")

## ----derivatives--------------------------------------------------------------
set.seed(99)
df_small <- data.frame(t = rexp(50, rate = 3))
model_small <- exponential_lifetime("t")
lambda_test <- 2.5

# Analytical score
score_fn <- score(model_small)
analytical_score <- score_fn(df_small, lambda_test)

# Numerical score (via numDeriv)
ll_fn <- loglik(model_small)
numerical_score <- numDeriv::grad(
  function(p) ll_fn(df_small, p), lambda_test
)

cat("Score at lambda =", lambda_test, ":\n")
cat("  Analytical:", analytical_score, "\n")
cat("  Numerical: ", numerical_score, "\n")
cat("  Match:", all.equal(unname(analytical_score), numerical_score), "\n")

# Analytical Hessian
hess_fn <- hess_loglik(model_small)
analytical_hess <- hess_fn(df_small, lambda_test)

numerical_hess <- numDeriv::hessian(
  function(p) ll_fn(df_small, p), lambda_test
)

cat("\nHessian at lambda =", lambda_test, ":\n")
cat("  Analytical:", analytical_hess[1, 1], "\n")
cat("  Numerical: ", numerical_hess[1, 1], "\n")
cat("  Match:", all.equal(analytical_hess[1, 1], numerical_hess[1, 1]), "\n")

## ----fim----------------------------------------------------------------------
fim_fn <- fim(model_small)
n_obs <- nrow(df_small)
lambda_hat <- coef(fit(model_small)(df_small))

fim_analytical <- fim_fn(lambda_hat, n_obs)

cat("FIM at MLE (lambda =", lambda_hat, "):\n")
cat("  Analytical n/lambda^2:", n_obs / lambda_hat^2, "\n")
cat("  fim() result:         ", fim_analytical[1, 1], "\n")
cat("  Match:", all.equal(n_obs / unname(lambda_hat)^2, fim_analytical[1, 1]), "\n")

## ----fisherian----------------------------------------------------------------
set.seed(42)
df_fish <- data.frame(t = rexp(100, rate = 2.0))
model_fish <- exponential_lifetime("t")
mle_fish <- fit(model_fish)(df_fish)

# Support function: S(lambda) = logL(lambda) - logL(lambda_hat)
# Support at MLE is always 0
s_mle <- support(mle_fish, coef(mle_fish), df_fish, model_fish)
cat("Support at MLE:", s_mle, "\n")

# At a different value, support is negative
s_alt <- support(mle_fish, c(lambda = 1.5), df_fish, model_fish)
cat("Support at lambda=1.5:", s_alt, "\n")

# Relative likelihood: R(lambda) = L(lambda)/L(lambda_hat) = exp(S)
rl_alt <- relative_likelihood(mle_fish, c(lambda = 1.5), df_fish, model_fish)
cat("Relative likelihood at lambda=1.5:", rl_alt, "\n")

## ----likelihood-interval------------------------------------------------------
li <- likelihood_interval(mle_fish, df_fish, model_fish, k = 8)
print(li)

cat("\nWald 95% CI for comparison:\n")
print(confint(mle_fish))

## ----profile, fig.width=6, fig.height=4---------------------------------------
prof <- profile_loglik(mle_fish, df_fish, model_fish, param = 1, n_grid = 100)

plot(prof$lambda, prof$support, type = "l", lwd = 2,
     xlab = expression(lambda), ylab = "Support",
     main = "Profile Support Function")
abline(h = -log(8), lty = 2, col = "red")
abline(v = coef(mle_fish), lty = 3, col = "blue")
legend("topright",
       legend = c("Support S(lambda)", "1/8 cutoff", "MLE"),
       lty = c(1, 2, 3), col = c("black", "red", "blue"), lwd = c(2, 1, 1))

## ----monte-carlo, cache=TRUE--------------------------------------------------
set.seed(42)
true_lambda <- 3.0
n_obs <- 100
n_sim <- 1000

model_mc <- exponential_lifetime("t")
gen <- rdata(model_mc)

# Simulate n_sim datasets and fit each
mle_vals <- replicate(n_sim, {
  sim_df <- gen(true_lambda, n_obs)
  coef(fit(model_mc)(sim_df))
})

cat("Monte Carlo results (", n_sim, "simulations, n=", n_obs, "):\n")
cat("  True lambda:       ", true_lambda, "\n")
cat("  Mean of MLEs:      ", mean(mle_vals), "\n")
cat("  Bias:              ", mean(mle_vals) - true_lambda, "\n")
cat("  Empirical SE:      ", sd(mle_vals), "\n")
cat("  Theoretical SE:    ", true_lambda / sqrt(n_obs), "\n")
cat("  (SE = lambda/sqrt(n) from Fisher information)\n")

## ----crossval-----------------------------------------------------------------
set.seed(42)
df_cv <- data.frame(t = rexp(100, rate = 2.0))

model_analytical <- exponential_lifetime("t")
model_generic <- likelihood_name("exp", ob_col = "t",
                                  censor_col = "censor")

# Need a censor column for likelihood_name
df_cv_generic <- data.frame(t = df_cv$t, censor = rep("exact", 100))

lambda_test <- 2.3
ll_analytical <- loglik(model_analytical)(df_cv, lambda_test)
ll_generic <- loglik(model_generic)(df_cv_generic, lambda_test)

cat("Log-likelihood at lambda =", lambda_test, ":\n")
cat("  Analytical model:", ll_analytical, "\n")
cat("  Generic model:   ", ll_generic, "\n")
cat("  Match:", all.equal(ll_analytical, ll_generic), "\n")

## ----session------------------------------------------------------------------
sessionInfo()

