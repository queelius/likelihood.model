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

## ----basic-normal-------------------------------------------------------------
# Generate data from N(5, 2)
set.seed(42)
n <- 200
df <- data.frame(x = rnorm(n, mean = 5, sd = 2))

# Create model -- no censor_col means all exact
model <- likelihood_name("norm", ob_col = "x")
print(model)

# Fit the MLE
mle <- fit(model)(df, par = c(mean = 0, sd = 1))
summary(mle)

## ----basic-inference----------------------------------------------------------
cat("Estimates:", coef(mle), "\n")
cat("Standard errors:", se(mle), "\n")
cat("AIC:", aic(mle), "\n")
confint(mle)

## ----right-censor-------------------------------------------------------------
# Simulate Weibull lifetimes with type-I right-censoring at time 1.2
set.seed(123)
n <- 150
true_shape <- 2
true_scale <- 1.5
censor_time <- 1.2

raw_t <- rweibull(n, true_shape, true_scale)
df_right <- data.frame(
  x = pmin(raw_t, censor_time),
  censor = ifelse(raw_t > censor_time, "right", "exact")
)

cat("Censoring rate:", mean(df_right$censor == "right") * 100, "%\n")

# Fit with proper censoring handling
model_right <- likelihood_name("weibull", ob_col = "x", censor_col = "censor")
mle_right <- fit(model_right)(df_right, par = c(shape = 1.5, scale = 1))

cat("\nMLE (with censoring):\n")
cat("  shape:", coef(mle_right)[1], "(true:", true_shape, ")\n")
cat("  scale:", coef(mle_right)[2], "(true:", true_scale, ")\n")
confint(mle_right)

## ----left-censor--------------------------------------------------------------
# Simulate concentrations from a log-normal distribution
# with detection limit at 0.5
set.seed(456)
n <- 100
true_meanlog <- 0
true_sdlog <- 0.8
detect_limit <- 0.5

raw_conc <- rlnorm(n, true_meanlog, true_sdlog)
df_left <- data.frame(
  x = ifelse(raw_conc < detect_limit, detect_limit, raw_conc),
  censor = ifelse(raw_conc < detect_limit, "left", "exact")
)

cat("Below detection limit:", mean(df_left$censor == "left") * 100, "%\n")

model_left <- likelihood_name("lnorm", ob_col = "x", censor_col = "censor")
mle_left <- fit(model_left)(df_left, par = c(meanlog = 0, sdlog = 1))

cat("\nMLE (accounting for detection limit):\n")
cat("  meanlog:", coef(mle_left)[1], "(true:", true_meanlog, ")\n")
cat("  sdlog:  ", coef(mle_left)[2], "(true:", true_sdlog, ")\n")

## ----interval-censor----------------------------------------------------------
# Simulate inspection data: items are checked at fixed intervals
# and we only know failure occurred between two inspection times
set.seed(789)
n <- 120
true_shape <- 3
true_scale <- 5

raw_t <- rweibull(n, true_shape, true_scale)
# Inspection every 1 unit of time
inspection_lower <- floor(raw_t)
inspection_upper <- ceiling(raw_t)
# Ensure lower < upper
inspection_upper <- ifelse(inspection_lower == inspection_upper,
                           inspection_upper + 1, inspection_upper)

df_interval <- data.frame(
  x = inspection_lower,
  x_upper = inspection_upper,
  censor = rep("interval", n)
)

head(df_interval)

model_interval <- likelihood_name("weibull", ob_col = "x",
                                  censor_col = "censor",
                                  ob_col_upper = "x_upper")
mle_interval <- fit(model_interval)(df_interval, par = c(shape = 2, scale = 4))

cat("MLE from interval-censored data:\n")
cat("  shape:", coef(mle_interval)[1], "(true:", true_shape, ")\n")
cat("  scale:", coef(mle_interval)[2], "(true:", true_scale, ")\n")

## ----mixed-censor-------------------------------------------------------------
set.seed(101)
n <- 200
true_mean <- 10
true_sd <- 3

raw_x <- rnorm(n, true_mean, true_sd)

# Simulate mixed censoring:
# - 60% exact observations
# - 15% right-censored (above upper limit)
# - 15% left-censored (below detection limit)
# - 10% interval-censored (binned to nearest integer)
detect_lower <- 5
detect_upper <- 14

censor_type <- character(n)
x_obs <- numeric(n)
x_upper <- rep(NA_real_, n)

for (i in seq_len(n)) {
  if (raw_x[i] < detect_lower) {
    censor_type[i] <- "left"
    x_obs[i] <- detect_lower
  } else if (raw_x[i] > detect_upper) {
    censor_type[i] <- "right"
    x_obs[i] <- detect_upper
  } else if (runif(1) < 0.15) {
    # Some exact observations get binned
    censor_type[i] <- "interval"
    x_obs[i] <- floor(raw_x[i])
    x_upper[i] <- ceiling(raw_x[i])
    if (x_obs[i] == x_upper[i]) x_upper[i] <- x_upper[i] + 1
  } else {
    censor_type[i] <- "exact"
    x_obs[i] <- raw_x[i]
  }
}

df_mixed <- data.frame(x = x_obs, x_upper = x_upper, censor = censor_type)
cat("Censoring breakdown:\n")
print(table(df_mixed$censor))

model_mixed <- likelihood_name("norm", ob_col = "x", censor_col = "censor",
                               ob_col_upper = "x_upper")
mle_mixed <- fit(model_mixed)(df_mixed, par = c(mean = 8, sd = 2))

cat("\nMLE with mixed censoring:\n")
cat("  mean:", coef(mle_mixed)[1], "(true:", true_mean, ")\n")
cat("  sd:  ", coef(mle_mixed)[2], "(true:", true_sd, ")\n")
confint(mle_mixed)

## ----multi-dist---------------------------------------------------------------
set.seed(202)

# Exponential
df_exp <- data.frame(x = rexp(200, rate = 2.5))
mle_exp <- fit(likelihood_name("exp", "x"))(df_exp, par = c(rate = 1))
cat("Exponential: rate =", coef(mle_exp), "(true: 2.5)\n")

# Gamma
df_gam <- data.frame(x = rgamma(200, shape = 3, rate = 2))
mle_gam <- fit(likelihood_name("gamma", "x"))(df_gam, par = c(shape = 1, rate = 1))
cat("Gamma: shape =", coef(mle_gam)[1], "rate =", coef(mle_gam)[2],
    "(true: 3, 2)\n")

# Log-normal
df_lnorm <- data.frame(x = rlnorm(200, meanlog = 1, sdlog = 0.5))
mle_lnorm <- fit(likelihood_name("lnorm", "x"))(df_lnorm,
                                                  par = c(meanlog = 0, sdlog = 1))
cat("Log-normal: meanlog =", coef(mle_lnorm)[1], "sdlog =", coef(mle_lnorm)[2],
    "(true: 1, 0.5)\n")

## ----fisherian----------------------------------------------------------------
# Fit a normal model
set.seed(303)
df <- data.frame(x = rnorm(150, mean = 5, sd = 2))
model <- likelihood_name("norm", ob_col = "x")
mle <- fit(model)(df, par = c(0, 1))

# Support function: S(theta) = logL(theta) - logL(theta_hat)
# At the MLE, support is always 0
s_mle <- support(mle, coef(mle), df, model)
cat("Support at MLE:", s_mle, "\n")

# At other values, support is negative
s_alt <- support(mle, c(4, 3), df, model)
cat("Support at (4, 3):", s_alt, "\n")

# Relative likelihood: R(theta) = L(theta)/L(theta_hat) = exp(S(theta))
rl <- relative_likelihood(mle, c(4, 3), df, model)
cat("Relative likelihood at (4, 3):", rl, "\n")

## ----likelihood-interval------------------------------------------------------
# 1/8 likelihood interval for the mean (parameter 1)
li <- likelihood_interval(mle, df, model, k = 8, param = 1)
cat("1/8 likelihood interval for mean:\n")
print(li)

# Compare with Wald CI
cat("\nWald 95% CI:\n")
print(confint(mle))

## ----profile------------------------------------------------------------------
prof <- profile_loglik(mle, df, model, param = 1, n_grid = 30)

plot(prof[[1]], prof$relative_likelihood, type = "l",
     xlab = "Mean", ylab = "Relative likelihood",
     main = "Profile likelihood for mean parameter")
abline(h = 1/8, lty = 2, col = "gray")
abline(v = coef(mle)[1], lty = 2, col = "blue")
legend("topright", legend = c("Profile R(mean)", "1/8 cutoff", "MLE"),
       lty = c(1, 2, 2), col = c("black", "gray", "blue"), cex = 0.8)

## ----comparison---------------------------------------------------------------
set.seed(404)
x <- rexp(150, rate = 2.0)

# Generic approach
df_gen <- data.frame(t = x, censor = rep("exact", length(x)))
model_gen <- likelihood_name("exp", ob_col = "t", censor_col = "censor")
ll_gen <- loglik(model_gen)(df_gen, 2.0)

# Specialized approach
df_spec <- data.frame(t = x)
model_spec <- exponential_lifetime("t")
ll_spec <- loglik(model_spec)(df_spec, 2.0)

cat("Generic loglik:     ", ll_gen, "\n")
cat("Specialized loglik: ", ll_spec, "\n")
cat("Match:", all.equal(ll_gen, ll_spec), "\n")

## ----comparison-fit-----------------------------------------------------------
# Generic: uses optim
mle_gen <- fit(model_gen)(df_gen, par = c(rate = 1))

# Specialized: closed-form MLE, no optimization needed
mle_spec <- fit(model_spec)(df_spec)

cat("Generic MLE:     ", coef(mle_gen), "\n")
cat("Specialized MLE: ", coef(mle_spec), "\n")

## ----param-naming-------------------------------------------------------------
model <- likelihood_name("norm", ob_col = "x")
ll <- loglik(model)
df <- data.frame(x = rnorm(50, 5, 2))

# Both produce identical results:
ll_named <- ll(df, c(mean = 5, sd = 2))
ll_unnamed <- ll(df, c(5, 2))
cat("Named:", ll_named, "\n")
cat("Unnamed:", ll_unnamed, "\n")
cat("Match:", ll_named == ll_unnamed, "\n")

## ----param-formals------------------------------------------------------------
# Gamma has formals: shape, rate (or scale)
# dnorm has formals: mean, sd
cat("dnorm formals:", paste(names(formals(dnorm)), collapse = ", "), "\n")
cat("dgamma formals:", paste(names(formals(dgamma)), collapse = ", "), "\n")
cat("dweibull formals:", paste(names(formals(dweibull)), collapse = ", "), "\n")

## ----param-error, error=TRUE--------------------------------------------------
try({
ll(df, c(0, 1, 0.5))  # 3 params for 2-param distribution
})

## ----lrt----------------------------------------------------------------------
set.seed(505)
df_test <- data.frame(x = rweibull(200, shape = 1.8, scale = 2))

# Fit both models
model_exp <- likelihood_name("exp", ob_col = "x")
model_weib <- likelihood_name("weibull", ob_col = "x")

mle_exp_test <- fit(model_exp)(df_test, par = c(rate = 0.5))
mle_weib_test <- fit(model_weib)(df_test, par = c(shape = 1, scale = 1))

cat("Exponential loglik:", loglik_val(mle_exp_test), "\n")
cat("Weibull loglik:    ", loglik_val(mle_weib_test), "\n")

# LRT: exponential (1 param) vs Weibull (2 params)
lrt_result <- lrt(model_exp, model_weib, df_test,
                  null_par = coef(mle_exp_test),
                  alt_par = coef(mle_weib_test),
                  dof = 1)
print(lrt_result)

