library(testthat)
library(likelihood.model)

# Helper function to get parameter estimates from MLE result
# fisher_mle stores estimates in $par, accessible via coef()
get_mle_params <- function(mle_result) {
  # Use coef() for fisher_mle objects
  result <- coef(mle_result)
  if (!is.null(result)) {
    return(result)
  }
  # Fallback: try $par directly
  if (!is.null(mle_result$par)) {
    return(mle_result$par)
  }
  return(NULL)
}

# =============================================================================
# IS_LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("is_likelihood_model correctly identifies likelihood models", {
  # exponential_lifetime is a likelihood model
  model_exp <- exponential_lifetime("t")
  expect_true(is_likelihood_model(model_exp))

  # Non-likelihood model objects return FALSE
  expect_false(is_likelihood_model(list()))
  expect_false(is_likelihood_model(data.frame()))
  expect_false(is_likelihood_model(NULL))
  expect_false(is_likelihood_model("not a model"))
  expect_false(is_likelihood_model(42))
})

# =============================================================================
# FIT.LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("fit.likelihood_model finds correct MLE for exponential data", {
  set.seed(999)
  true_rate <- 2.0
  n <- 200
  df <- data.frame(t = rexp(n, rate = true_rate))

  model <- exponential_lifetime("t")
  result <- fit(model)(df)

  estimated <- get_mle_params(result)

  # Closed-form MLE: lambda_hat = n / sum(t)
  expect_equal(unname(estimated), n / sum(df$t), tolerance = 1e-12)

  # Should be close to true value
  expect_lt(abs(estimated[1] - true_rate), 0.5)
})

# =============================================================================
# SAMPLER.LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("sampler.likelihood_model generates bootstrap samples", {
  skip_if_not_installed("boot")

  set.seed(1002)
  df <- data.frame(t = rexp(200, rate = 2))

  model <- exponential_lifetime("t")

  # Create sampler
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))

  # Generate bootstrap samples
  boot_result <- model_sampler(n = 10)

  # Should return a fisher_boot object
  expect_true(inherits(boot_result, "fisher_boot"))
  expect_true(inherits(boot_result, "fisher_mle"))

  # Bootstrap estimates should be stored in replicates
  expect_equal(nrow(boot_result$replicates), 10)
  expect_equal(ncol(boot_result$replicates), 1)  # rate only

  # coef() should return original MLE
  expect_equal(length(coef(boot_result)), 1)
})

test_that("sampler.likelihood_model bootstrap distribution is centered near MLE", {
  skip_if_not_installed("boot")

  set.seed(1003)
  true_rate <- 2.0
  df <- data.frame(t = rexp(200, rate = true_rate))

  model <- exponential_lifetime("t")

  # First fit to get MLE
  mle_result <- fit(model)(df)
  mle_params <- get_mle_params(mle_result)

  # Create sampler and generate bootstrap samples
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_result <- model_sampler(n = 50)

  # Bootstrap mean should be close to MLE
  boot_mean <- colMeans(boot_result$replicates)
  expect_equal(unname(boot_mean[1]), unname(mle_params[1]), tolerance = 0.5)
})

# =============================================================================
# LRT (LIKELIHOOD RATIO TEST) TESTS
# =============================================================================

test_that("lrt requires likelihood models as inputs", {
  df <- data.frame(t = rexp(10, rate = 1))

  expect_error(lrt(list(), exponential_lifetime("t"),
                   df, null_par = c(1), alt_par = c(2)))
})

test_that("lrt computes correct test statistic", {
  set.seed(1004)
  df <- data.frame(t = rexp(100, rate = 1))

  model <- exponential_lifetime("t")
  ll_func <- loglik(model)

  null_par <- c(1.0)
  alt_par <- c(1.05)

  null_ll <- ll_func(df, null_par)
  alt_ll <- ll_func(df, alt_par)

  result <- lrt(model, model, df,
                null_par = null_par, alt_par = alt_par, dof = 1)

  expected_stat <- -2 * (null_ll - alt_ll)
  expect_equal(result$stat, expected_stat, tolerance = 1e-10)
  expect_equal(result$dof, 1)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("lrt p-value is computed correctly", {
  set.seed(1005)
  df <- data.frame(t = rexp(200, rate = 2))

  model <- exponential_lifetime("t")

  null_par <- c(2.0)
  mle <- fit(model)(df)
  alt_par <- coef(mle)

  result <- lrt(model, model, df,
                null_par = null_par, alt_par = alt_par, dof = 1)

  # Test statistic should be non-negative
  expect_true(result$stat >= 0)
  # p-value should be in [0, 1]
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("lrt auto-calculates degrees of freedom", {
  set.seed(1006)
  df <- data.frame(t = rexp(50, rate = 2))

  model <- exponential_lifetime("t")

  result <- lrt(model, model, df,
                null_par = c(1),
                alt_par = c(1, 0.5),
                dof = NULL)

  expect_equal(result$dof, 1)
  expect_true(!is.null(result$stat))
  expect_true(!is.null(result$p.value))
})

# =============================================================================
# EDGE CASES AND ERROR HANDLING
# =============================================================================

test_that("loglik handles data with single observation", {
  df <- data.frame(t = 5.0)

  model <- exponential_lifetime("t")
  ll_func <- loglik(model)

  expected_ll <- dexp(5.0, rate = 2, log = TRUE)
  computed_ll <- ll_func(df, c(2))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

# =============================================================================
# PRINT.LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("print.likelihood_model outputs model information", {
  model <- exponential_lifetime("t")

  # Capture output
  output <- capture.output(print(model))

  # Check that key information is printed
  expect_true(any(grepl("Likelihood model", output)))
  expect_true(any(grepl("Observation column", output)))
})

test_that("print.likelihood_model with show.loglik=TRUE shows loglik function", {
  model <- exponential_lifetime("t")

  output <- capture.output(print(model, show.loglik = TRUE))
  expect_true(any(grepl("Log-likelihood function", output)))
})

test_that("print.likelihood_model returns model invisibly", {
  model <- exponential_lifetime("t")

  result <- print(model)
  expect_identical(result, model)
})

# =============================================================================
# LRT WITH NULL DOF TESTS
# =============================================================================

test_that("lrt computes dof automatically when NULL", {
  set.seed(1014)
  df <- data.frame(t = rexp(100, rate = 2))

  model <- exponential_lifetime("t")

  null_par <- c(1.5)
  alt_par <- c(2.0)

  # When dof is NULL, it should be computed as length(alt_par) - length(null_par)
  result <- lrt(model, model, df,
                null_par = null_par, alt_par = alt_par, dof = NULL)

  # dof should be 0 (1 - 1)
  expect_equal(result$dof, 0)
  expect_true(!is.null(result$stat))
  expect_true(!is.null(result$p.value))
})

# =============================================================================
# FIT WITH SANN METHOD TESTS
# =============================================================================

test_that("fit.likelihood_model works with SANN method (no gradient)", {
  set.seed(1015)
  df <- data.frame(t = rexp(100, rate = 2))

  # Use the generic fit.likelihood_model (not exponential_lifetime's override)
  # by creating a mock model that delegates to exponential loglik
  mock_model <- structure(
    list(ob_col = "t", censor_col = NULL),
    class = c("mock_exp_sann", "likelihood_model")
  )
  loglik.mock_exp_sann <<- function(model, ...) {
    function(df, par, ...) {
      lambda <- par[1]
      if (lambda <= 0) return(-.Machine$double.xmax / 2)
      n <- nrow(df)
      n * log(lambda) - lambda * sum(df[[model$ob_col]])
    }
  }

  solver <- fit(mock_model)
  result <- solver(df, par = c(1), method = "SANN",
                   control = list(maxit = 1000))

  params_val <- get_mle_params(result)
  expect_true(!is.null(params_val))
  expect_true(!any(is.na(params_val)))

  # Parameters should be reasonable
  expect_lt(abs(params_val[1] - 2), 1.5)

  rm(loglik.mock_exp_sann, envir = .GlobalEnv)
})

# =============================================================================
# FIM GENERIC TESTS
# =============================================================================

test_that("fim generic dispatches correctly", {
  # fim is a generic that requires a method
  # Test that it exists and errors appropriately for objects without methods
  expect_error(fim(list()), "no applicable method")
})

# =============================================================================
# FISHER_MLE CLASS TESTS
# =============================================================================

test_that("fisher_mle creates object with correct structure", {
  result <- fisher_mle(
    par = c(mean = 0, sd = 1),
    loglik_val = -50.5,
    hessian = matrix(c(-100, 0, 0, -50), 2, 2),
    nobs = 100
  )

  expect_true(inherits(result, "fisher_mle"))
  expect_true(inherits(result, "mle_fit"))

  # Test base R generic methods

  expect_equal(coef(result), c(mean = 0, sd = 1))
  expect_true(!is.null(vcov(result)))
  expect_equal(nobs(result), 100)

  # logLik should return a logLik object
  ll <- logLik(result)
  expect_true(inherits(ll, "logLik"))
  expect_equal(as.numeric(ll), -50.5)
})

test_that("fisher_mle confint computes Wald confidence intervals", {
  result <- fisher_mle(
    par = c(mean = 5, sd = 2),
    vcov = matrix(c(0.04, 0, 0, 0.02), 2, 2),  # SE = 0.2 and ~0.14
    loglik_val = -100,
    nobs = 100
  )

  ci <- confint(result, level = 0.95)

  # Check structure
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)

  # Intervals should be centered around estimates
  expect_lt(ci[1, 1], 5)  # lower < estimate
  expect_gt(ci[1, 2], 5)  # upper > estimate
})

test_that("fisher_mle se returns standard errors", {
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    vcov = matrix(c(0.25, 0, 0, 0.16), 2, 2),  # SE = 0.5 and 0.4
    loglik_val = -50
  )

  se_vals <- se(result)
  expect_equal(se_vals, c(0.5, 0.4))
})

test_that("fisher_mle AIC and BIC compute correctly via stats generics", {
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    loglik_val = -100,
    nobs = 50
  )

  # AIC = -2*logL + 2*k = -2*(-100) + 2*2 = 204
  expect_equal(stats::AIC(result), 204)

  # BIC = -2*logL + k*log(n) = 200 + 2*log(50)
  expected_bic <- 200 + 2 * log(50)
  expect_equal(stats::BIC(result), expected_bic)
})

test_that("fisher_mle summary produces correct output", {
  result <- fisher_mle(
    par = c(mean = 5, sd = 2),
    vcov = matrix(c(0.04, 0, 0, 0.02), 2, 2),
    loglik_val = -100,
    nobs = 100
  )

  summ <- summary(result)
  expect_true(inherits(summ, "summary_fisher_mle"))
  expect_true(!is.null(summ$coefficients))
  expect_equal(summ$loglik, -100)
  expect_equal(summ$nobs, 100)
})

test_that("fit.likelihood_model returns fisher_mle object", {
  set.seed(2001)
  df <- data.frame(t = rexp(100, rate = 2))

  model <- exponential_lifetime("t")
  result <- fit(model)(df)

  expect_true(inherits(result, "fisher_mle"))
  expect_true(inherits(result, "mle_fit"))
  expect_true(result$converged)
  expect_equal(result$nobs, 100)
})

# =============================================================================
# FISHERIAN INFERENCE TESTS
# =============================================================================

test_that("support function computes log relative likelihood", {
  set.seed(2002)
  df <- data.frame(t = rexp(100, rate = 2))

  model <- exponential_lifetime("t")
  mle_result <- fit(model)(df)

  # Support at MLE should be 0
  s_at_mle <- support(mle_result, coef(mle_result), df, model)
  expect_equal(unname(s_at_mle), 0, tolerance = 1e-10)

  # Support at other values should be negative
  s_away <- support(mle_result, c(lambda = 0.5), df, model)
  expect_true(s_away < 0)
})

test_that("relative_likelihood computes exp(support)", {
  set.seed(2003)
  df <- data.frame(t = rexp(50, rate = 2))

  model <- exponential_lifetime("t")
  mle_result <- fit(model)(df)

  # Relative likelihood at MLE should be 1
  rl_at_mle <- relative_likelihood(mle_result, coef(mle_result), df, model)
  expect_equal(unname(rl_at_mle), 1, tolerance = 1e-10)

  # Relative likelihood elsewhere should be between 0 and 1
  rl_away <- relative_likelihood(mle_result, c(lambda = 0.5), df, model)
  expect_true(rl_away > 0 && rl_away < 1)
})

test_that("likelihood_interval works for single parameter model", {
  set.seed(2004)
  df <- data.frame(t = rexp(200, rate = 2))

  model <- exponential_lifetime("t")
  mle_result <- fit(model)(df)

  # Compute 1/8 likelihood interval
  li <- likelihood_interval(mle_result, df, model, k = 8)

  expect_true(inherits(li, "likelihood_interval"))
  expect_equal(nrow(li), 1)
  expect_equal(ncol(li), 2)

  # Interval should contain MLE
  mle_val <- unname(coef(mle_result))
  expect_true(li[1, 1] < mle_val && mle_val < li[1, 2])
})

test_that("profile_loglik computes profile for single parameter", {
  set.seed(2005)
  df <- data.frame(t = rexp(100, rate = 2))

  model <- exponential_lifetime("t")
  mle_result <- fit(model)(df)

  # Profile for first parameter (lambda)
  prof <- profile_loglik(mle_result, df, model, param = 1, n_grid = 20)

  expect_true(inherits(prof, "profile_loglik"))
  expect_equal(nrow(prof), 20)
  expect_true("loglik" %in% names(prof))
  expect_true("support" %in% names(prof))
  expect_true("relative_likelihood" %in% names(prof))

  # Maximum profile loglik should be near MLE
  max_idx <- which.max(prof$loglik)
  expect_equal(unname(prof[[1]][max_idx]), unname(coef(mle_result)[1]), tolerance = 0.5)
})

test_that("deviance computes correctly for fisher_mle", {
  result <- fisher_mle(
    par = c(a = 1),
    loglik_val = -50,
    nobs = 100
  )

  # Single model deviance = -2 * logL = 100
  expect_equal(deviance(result), 100)

  # Comparison deviance
  null_result <- fisher_mle(
    par = c(a = 0.5),
    loglik_val = -55,
    nobs = 100
  )

  # Deviance = 2*(logL_full - logL_null) = 2*(-50 - (-55)) = 10
  expect_equal(deviance(result, null_result), 10)
})

test_that("evidence computes log likelihood ratio", {
  set.seed(2006)
  df <- data.frame(t = rexp(50, rate = 2))

  model <- exponential_lifetime("t")

  # Evidence for theta1 vs theta2
  ev <- evidence(model, df, c(lambda = 2), c(lambda = 0.5))

  # Should strongly favor true parameters (2) over (0.5)
  expect_true(ev > 0)
  expect_true(ev > log(8))  # Strong evidence
})

# =============================================================================
# FISHER_BOOT CLASS TESTS
# =============================================================================

test_that("fisher_boot inherits from fisher_mle", {
  skip_if_not_installed("boot")

  set.seed(2007)
  df <- data.frame(t = rexp(100, rate = 2))

  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_result <- model_sampler(n = 20)

  expect_true(inherits(boot_result, "fisher_boot"))
  expect_true(inherits(boot_result, "fisher_mle"))
  expect_true(inherits(boot_result, "mle_fit"))

  # Base methods should work
  expect_equal(length(coef(boot_result)), 1)
  expect_true(is.matrix(vcov(boot_result)))
})

test_that("fisher_boot bias computes from bootstrap replicates", {
  skip_if_not_installed("boot")

  set.seed(2008)
  df <- data.frame(t = rexp(200, rate = 2))

  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_result <- model_sampler(n = 50)

  # Bias should be computed from replicates
  b <- bias(boot_result)
  expect_equal(length(b), 1)

  # Expected: mean(replicates) - original MLE
  expected_bias <- colMeans(boot_result$replicates) - coef(boot_result)
  expect_equal(b, expected_bias)
})

test_that("sampler generic works on fisher_mle for asymptotic sampling", {
  result <- fisher_mle(
    par = c(mean = 5, sd = 2),
    vcov = matrix(c(0.04, 0, 0, 0.02), 2, 2),
    loglik_val = -100,
    nobs = 100
  )

  samp_fn <- sampler(result)
  samples <- samp_fn(1000)

  expect_equal(dim(samples), c(1000, 2))

  # Samples should be centered near MLE
  expect_equal(mean(samples[, 1]), 5, tolerance = 0.1)
  expect_equal(mean(samples[, 2]), 2, tolerance = 0.1)
})

# =============================================================================
# EDGE CASE AND BUG FIX VALIDATION TESTS
# =============================================================================

test_that("lrt with dof=NULL computes dof from parameter lengths", {
  set.seed(3001)
  df <- data.frame(t = rexp(100, rate = 2))

  model <- exponential_lifetime("t")

  # alt_par has 2 elements, null_par has 1 => dof should be 1
  null_par <- c(2)
  alt_par <- c(2, 0.5)

  result <- lrt(model, model, df,
                null_par = null_par, alt_par = alt_par, dof = NULL)

  expect_equal(result$dof, 1)
  expect_true(is.finite(result$stat))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("evidence returns 0 when comparing same parameters", {
  set.seed(3003)
  df <- data.frame(t = rexp(50, rate = 2))

  model <- exponential_lifetime("t")

  # Evidence for identical parameters should be exactly 0
  ev <- evidence(model, df, c(lambda = 2), c(lambda = 2))
  expect_equal(unname(ev), 0)
})

test_that("fisher_mle with NULL vcov: se returns NA, sampler errors", {
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    loglik_val = -50,
    nobs = 100
    # vcov and hessian both NULL
  )

  # se should return NA values
  se_vals <- se(result)
  expect_equal(length(se_vals), 2)
  expect_true(all(is.na(se_vals)))

  # sampler should error cleanly
  expect_error(sampler(result), "Cannot create sampler: vcov is NULL")
})

# =============================================================================
# EXPONENTIAL_LIFETIME TESTS
# =============================================================================

test_that("exponential_lifetime creates model with correct structure", {
  model <- exponential_lifetime("t")
  expect_true(is_likelihood_model(model))
  expect_equal(model$ob_col, "t")
  expect_null(model$censor_col)
  expect_true("exponential_lifetime" %in% class(model))

  model_c <- exponential_lifetime("t", censor_col = "status")
  expect_equal(model_c$censor_col, "status")
})

test_that("loglik matches manual d*log(lambda) - lambda*T", {
  set.seed(4001)
  lambda <- 2.0
  x <- rexp(100, rate = lambda)
  df <- data.frame(t = x)

  model <- exponential_lifetime("t")
  ll_func <- loglik(model)

  # Manual calculation
  n <- length(x)
  total_time <- sum(x)
  manual_ll <- n * log(lambda) - lambda * total_time

  expect_equal(ll_func(df, lambda), manual_ll, tolerance = 1e-12)
})

test_that("loglik matches sum of dexp(x, rate, log=TRUE)", {
  set.seed(4002)
  lambda <- 1.5
  x <- rexp(80, rate = lambda)
  df <- data.frame(t = x)

  model <- exponential_lifetime("t")
  ll_func <- loglik(model)

  expected_ll <- sum(dexp(x, rate = lambda, log = TRUE))
  expect_equal(ll_func(df, lambda), expected_ll, tolerance = 1e-10)
})

test_that("closed-form MLE equals d/T exactly", {
  set.seed(4003)
  lambda_true <- 3.0
  n <- 200
  x <- rexp(n, rate = lambda_true)
  df <- data.frame(t = x)

  model <- exponential_lifetime("t")
  result <- fit(model)(df)

  # Closed-form MLE: lambda_hat = n / sum(x)
  expected_mle <- n / sum(x)
  expect_equal(unname(coef(result)), expected_mle, tolerance = 1e-12)

  # Should be close to true value
  expect_lt(abs(coef(result) - lambda_true), 0.5)

  # MLE should have converged
  expect_true(result$converged)
  expect_equal(result$nobs, n)
})

test_that("analytical score matches numDeriv::grad", {
  set.seed(4004)
  lambda <- 2.5
  x <- rexp(50, rate = lambda)
  df <- data.frame(t = x)

  model <- exponential_lifetime("t")
  ll_func <- loglik(model)
  score_func <- score(model)

  # Test at several parameter values
  for (lam in c(1.0, 2.0, 3.0, 5.0)) {
    analytical <- score_func(df, lam)
    numerical <- numDeriv::grad(function(p) ll_func(df, p), lam)
    expect_equal(unname(analytical), numerical, tolerance = 1e-6,
                 info = paste("Score mismatch at lambda =", lam))
  }
})

test_that("analytical hessian matches numDeriv::hessian", {
  set.seed(4005)
  lambda <- 2.0
  x <- rexp(100, rate = lambda)
  df <- data.frame(t = x)

  model <- exponential_lifetime("t")
  ll_func <- loglik(model)
  hess_func <- hess_loglik(model)

  for (lam in c(1.0, 2.0, 4.0)) {
    analytical <- hess_func(df, lam)
    numerical <- numDeriv::hessian(function(p) ll_func(df, p), lam)
    expect_equal(as.numeric(analytical), as.numeric(numerical), tolerance = 1e-4,
                 info = paste("Hessian mismatch at lambda =", lam))
  }
})

test_that("right-censored MLE is correct", {
  set.seed(4006)
  lambda_true <- 2.0
  n <- 200
  censor_time <- 0.4

  x_raw <- rexp(n, rate = lambda_true)
  x_obs <- pmin(x_raw, censor_time)
  status <- ifelse(x_raw > censor_time, "right", "exact")

  df <- data.frame(t = x_obs, status = status)
  model <- exponential_lifetime("t", censor_col = "status")

  result <- fit(model)(df)
  lambda_hat <- unname(coef(result))

  # Verify closed-form: d / T
  d <- sum(status == "exact")
  total_time <- sum(x_obs)
  expect_equal(lambda_hat, d / total_time, tolerance = 1e-12)

  # Should be reasonable (right-censoring still gives consistent estimator)
  expect_lt(abs(lambda_hat - lambda_true), 0.8)
})

test_that("score is exactly zero at MLE", {
  set.seed(4007)
  x <- rexp(100, rate = 1.5)
  df <- data.frame(t = x)

  model <- exponential_lifetime("t")
  result <- fit(model)(df)

  # Score at MLE should be exactly 0 (closed-form solution)
  score_func <- score(model)
  s <- score_func(df, coef(result))
  expect_equal(unname(s), 0, tolerance = 1e-12)

  # Also stored in the MLE object
  expect_equal(unname(score_val(result)), 0, tolerance = 1e-12)
})

test_that("FIM equals n/lambda^2", {
  model <- exponential_lifetime("t")
  fim_func <- fim(model)

  # FIM for n observations at lambda
  for (lam in c(0.5, 1.0, 2.0, 5.0)) {
    for (n_obs in c(10, 100, 1000)) {
      expected_fim <- matrix(n_obs / lam^2, 1, 1)
      computed_fim <- fim_func(lam, n_obs)
      expect_equal(as.numeric(computed_fim), as.numeric(expected_fim),
                   tolerance = 1e-12,
                   info = paste("FIM mismatch at lambda =", lam, "n =", n_obs))
    }
  }
})

test_that("full Fisherian inference works with exponential_lifetime", {
  set.seed(4008)
  lambda_true <- 2.0
  n <- 200
  x <- rexp(n, rate = lambda_true)
  df <- data.frame(t = x)

  model <- exponential_lifetime("t")
  mle_result <- fit(model)(df)

  # Support at MLE should be 0
  s_at_mle <- support(mle_result, coef(mle_result), df, model)
  expect_equal(unname(s_at_mle), 0, tolerance = 1e-10)

  # Support away from MLE should be negative
  s_away <- support(mle_result, c(lambda = 0.5), df, model)
  expect_true(unname(s_away) < 0)

  # Relative likelihood at MLE is 1
  rl_at_mle <- relative_likelihood(mle_result, coef(mle_result), df, model)
  expect_equal(unname(rl_at_mle), 1, tolerance = 1e-10)

  # Likelihood interval should contain MLE
  li <- likelihood_interval(mle_result, df, model, k = 8)
  mle_val <- unname(coef(mle_result))
  expect_true(!is.na(li[1, 1]) && !is.na(li[1, 2]))
  expect_true(li[1, 1] < mle_val && mle_val < li[1, 2])
})

test_that("exponential_lifetime errors with no exact observations", {
  df <- data.frame(t = c(1, 2, 3), status = rep("right", 3))
  model <- exponential_lifetime("t", censor_col = "status")

  expect_error(fit(model)(df), "no exact")
})

test_that("exponential_lifetime rdata generates correct data", {
  model <- exponential_lifetime("t")
  rdata_fn <- rdata(model)

  set.seed(4010)
  df <- rdata_fn(theta = 2.0, n = 100)
  expect_equal(nrow(df), 100)
  expect_true("t" %in% names(df))
  expect_true(all(df$t > 0))

  # With censoring
  model_c <- exponential_lifetime("t", censor_col = "status")
  rdata_fn_c <- rdata(model_c)
  df_c <- rdata_fn_c(theta = 2.0, n = 100, censor_time = 0.5)
  expect_equal(nrow(df_c), 100)
  expect_true("t" %in% names(df_c))
  expect_true("status" %in% names(df_c))
  expect_true(all(df_c$t <= 0.5))
  expect_true(all(df_c$status %in% c("exact", "right")))
})

test_that("assumptions include censoring note when censor_col is set", {
  model <- exponential_lifetime("t")
  a <- assumptions(model)
  expect_false("non-informative right censoring" %in% a)

  model_c <- exponential_lifetime("t", censor_col = "status")
  a_c <- assumptions(model_c)
  expect_true("non-informative right censoring" %in% a_c)
})

# =============================================================================
# ALGEBRAIC.MLE INTERFACE COMPATIBILITY TESTS
# =============================================================================

test_that("params() returns same as coef() on fisher_mle", {
  result <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    loglik_val = -50,
    hessian = matrix(c(-100, 0, 0, -50), 2, 2),
    nobs = 100
  )

  expect_equal(params(result), coef(result))
  expect_equal(params(result), c(shape = 2.0, scale = 1.5))
})

test_that("nparams() returns correct count", {
  # Multivariate
  result2 <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    loglik_val = -50,
    nobs = 100
  )
  expect_equal(nparams(result2), 2L)

  # Univariate
  result1 <- fisher_mle(
    par = c(lambda = 3.0),
    loglik_val = -30,
    nobs = 50
  )
  expect_equal(nparams(result1), 1L)
})

test_that("observed_fim() returns -hessian (positive definite)", {
  H <- matrix(c(-100, -5, -5, -50), 2, 2)
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    loglik_val = -50,
    hessian = H,
    nobs = 100
  )

  fim_mat <- observed_fim(result)
  expect_equal(fim_mat, -H)

  # Should be positive definite (all eigenvalues > 0)
  evals <- eigen(fim_mat)$values
  expect_true(all(evals > 0))
})

test_that("observed_fim() returns NULL when hessian is NULL", {
  result <- fisher_mle(
    par = c(a = 1),
    loglik_val = -50,
    nobs = 100
    # hessian and vcov both NULL
  )

  expect_null(observed_fim(result))
})

test_that("obs() returns NULL for fisher_mle (by design)", {
  result <- fisher_mle(
    par = c(a = 1),
    loglik_val = -50,
    nobs = 100
  )
  expect_null(obs(result))
})

test_that("mse() equals vcov under zero asymptotic bias", {
  V <- matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2)
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    vcov = V,
    loglik_val = -50,
    nobs = 100
  )

  # Bias is zero for fisher_mle, so MSE = vcov
  expect_equal(mse(result), V)
})

test_that("mse() works for univariate fisher_mle", {
  result <- fisher_mle(
    par = c(lambda = 2.0),
    vcov = matrix(0.04, 1, 1),
    loglik_val = -30,
    nobs = 50
  )

  # Scalar case: MSE = var + 0^2 = var
  expect_equal(mse(result), matrix(0.04, 1, 1))
})

test_that("params()/nparams() work on fit() output from real models", {
  set.seed(5001)
  df <- data.frame(t = rexp(100, rate = 2))
  model <- exponential_lifetime("t")
  result <- fit(model)(df)

  p <- params(result)
  expect_equal(length(p), 1)
  expect_equal(p, coef(result))
  expect_equal(nparams(result), 1L)
})

test_that("observed_fim() is positive definite on real fit output", {
  set.seed(5002)
  df <- data.frame(t = rexp(200, rate = 2))
  model <- exponential_lifetime("t")
  result <- fit(model)(df)

  fim_mat <- observed_fim(result)
  expect_true(!is.null(fim_mat))
  expect_true(is.matrix(fim_mat))

  evals <- eigen(fim_mat)$values
  expect_true(all(evals > 0),
              info = "Observed FIM should be positive definite at MLE")
})

test_that("rmap() fallthrough works on fisher_mle", {
  skip_if_not_installed("algebraic.mle")

  result <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    vcov = matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2),
    loglik_val = -50,
    nobs = 100
  )

  # rmap transforms parameters via delta method
  # Transform: mean lifetime = scale * gamma(1 + 1/shape)
  mapped <- algebraic.mle::rmap(result, function(p) {
    c(mean_life = p[2] * gamma(1 + 1 / p[1]))
  })

  expect_true(inherits(mapped, "mle_fit"))
  expect_equal(length(params(mapped)), 1)
  expect_true(params(mapped) > 0)
})

test_that("marginal() fallthrough works on fisher_mle", {
  skip_if_not_installed("algebraic.mle")

  result <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    vcov = matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2),
    loglik_val = -50,
    nobs = 100
  )

  # Extract marginal for first parameter
  m <- algebraic.mle::marginal(result, 1)
  expect_true(inherits(m, "mle_fit"))
  expect_equal(nparams(m), 1L)
  expect_equal(unname(params(m)), 2.0)
})

test_that("expectation() fallthrough works on fisher_mle", {
  skip_if_not_installed("algebraic.mle")

  result <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    vcov = matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2),
    loglik_val = -50,
    nobs = 100
  )

  set.seed(5003)
  # Monte Carlo expectation of shape parameter
  # Use qualified name to avoid conflict with testthat::expectation
  # n goes in control list per algebraic.mle API
  e <- algebraic.mle::expectation(result, function(p) p[1],
                                  control = list(n = 5000L))
  expect_equal(e, 2.0, tolerance = 0.1)
})


# =============================================================================
# FISHER_MLE PRINT AND SUMMARY COVERAGE TESTS
# =============================================================================

test_that("print.fisher_mle outputs expected text", {
  result <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    vcov = matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2),
    loglik_val = -50,
    nobs = 100
  )

  out <- capture.output(print(result))
  expect_true(any(grepl("Maximum Likelihood Estimate", out)))
  expect_true(any(grepl("Log-likelihood", out)))
  expect_true(any(grepl("Observations:", out)))
})

test_that("print.fisher_mle shows warning when not converged", {
  result <- fisher_mle(
    par = c(a = 1),
    vcov = matrix(0.1, 1, 1),
    loglik_val = -50,
    nobs = 10,
    converged = FALSE
  )

  out <- capture.output(print(result))
  expect_true(any(grepl("WARNING.*converge", out)))
})

test_that("print.fisher_mle works without nobs", {
  result <- fisher_mle(
    par = c(a = 1),
    vcov = matrix(0.1, 1, 1),
    loglik_val = -50
  )

  out <- capture.output(print(result))
  expect_false(any(grepl("Observations:", out)))
})

test_that("print.summary_fisher_mle outputs expected text", {
  result <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    vcov = matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2),
    loglik_val = -50,
    nobs = 100
  )

  out <- capture.output(print(summary(result)))
  expect_true(any(grepl("Maximum Likelihood Estimate", out)))
  expect_true(any(grepl("AIC", out)))
  expect_true(any(grepl("Log-likelihood", out)))
  expect_true(any(grepl("Number of observations", out)))
})

test_that("print.summary_fisher_mle shows warning when not converged", {
  result <- fisher_mle(
    par = c(a = 1),
    vcov = matrix(0.1, 1, 1),
    loglik_val = -50,
    nobs = 10,
    converged = FALSE
  )

  out <- capture.output(print(summary(result)))
  expect_true(any(grepl("WARNING.*converge", out)))
})

test_that("print.summary_fisher_mle works without nobs", {
  result <- fisher_mle(
    par = c(a = 1),
    vcov = matrix(0.1, 1, 1),
    loglik_val = -50
  )

  out <- capture.output(print(summary(result)))
  expect_false(any(grepl("Number of observations", out)))
})


# =============================================================================
# FISHER_MLE ACCESSOR COVERAGE TESTS
# =============================================================================

test_that("logLik.fisher_mle returns the log-likelihood value", {
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    vcov = matrix(c(0.04, 0, 0, 0.02), 2, 2),
    loglik_val = -123.45,
    nobs = 50
  )
  expect_equal(as.numeric(logLik(result)), -123.45)
})

test_that("confint.fisher_mle works with character parm", {
  result <- fisher_mle(
    par = c(shape = 2.0, scale = 1.5),
    vcov = matrix(c(0.04, 0.01, 0.01, 0.02), 2, 2),
    loglik_val = -50,
    nobs = 100
  )

  ci_all <- confint(result)
  ci_shape <- confint(result, parm = "shape")
  expect_equal(nrow(ci_shape), 1)
  expect_equal(ci_shape[1, ], ci_all["shape", ])
})

test_that("fisher_mle warns when hessian is singular", {
  # Singular hessian: second row is a multiple of first
  H <- matrix(c(-1, -2, -2, -4), 2, 2)
  expect_warning(
    result <- fisher_mle(
      par = c(a = 1, b = 2),
      loglik_val = -50,
      hessian = H,
      nobs = 100
    ),
    "not invertible"
  )
  expect_null(result$vcov)
})

test_that("BIC works via stats::BIC when nobs is available", {
  result <- fisher_mle(
    par = c(a = 1),
    vcov = matrix(0.1, 1, 1),
    loglik_val = -50,
    nobs = 100
  )
  # BIC = -2*logL + k*log(n) = 100 + 1*log(100)
  expect_equal(stats::BIC(result), 100 + log(100))
})


# =============================================================================
# FISHER_BOOT TESTS
# =============================================================================

test_that("confint.fisher_boot computes percentile CI", {
  skip_if_not_installed("boot")
  set.seed(6001)
  df <- data.frame(t = rexp(100, rate = 2))
  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, level = 0.95, type = "perc")
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
  expect_true(ci[1, 1] < ci[1, 2])  # lower < upper
})

test_that("confint.fisher_boot computes basic CI", {
  skip_if_not_installed("boot")
  set.seed(6002)
  df <- data.frame(t = rexp(100, rate = 2))
  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, level = 0.95, type = "basic")
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
})

test_that("confint.fisher_boot computes norm CI", {
  skip_if_not_installed("boot")
  set.seed(6003)
  df <- data.frame(t = rexp(100, rate = 2))
  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, level = 0.95, type = "norm")
  expect_equal(nrow(ci), 1)
  expect_equal(ncol(ci), 2)
})

test_that("confint.fisher_boot with integer parm", {
  skip_if_not_installed("boot")
  set.seed(6004)
  df <- data.frame(t = rexp(100, rate = 2))
  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, parm = 1, level = 0.95, type = "perc")
  expect_equal(nrow(ci), 1)
})

test_that("print.fisher_boot outputs expected text", {
  skip_if_not_installed("boot")
  set.seed(6005)
  df <- data.frame(t = rexp(100, rate = 2))
  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_obj <- model_sampler(n = 50)

  out <- capture.output(print(boot_obj))
  expect_true(any(grepl("Bootstrap MLE", out)))
  expect_true(any(grepl("Bootstrap replicates:", out)))
  expect_true(any(grepl("Bootstrap SE:", out)))
})

test_that("sampler.fisher_boot resamples from bootstrap distribution", {
  skip_if_not_installed("boot")
  set.seed(6006)
  df <- data.frame(t = rexp(100, rate = 2))
  model <- exponential_lifetime("t")
  model_sampler <- sampler(model, df = df, par = c(lambda = 2))
  boot_obj <- model_sampler(n = 200)

  samp_fn <- sampler(boot_obj)
  samples <- samp_fn(100)
  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 1)
  expect_true(all(is.finite(samples)))
})


# =============================================================================
# SAMPLER.FISHER_MLE UNIVARIATE COVERAGE
# =============================================================================

test_that("sampler.fisher_mle works for univariate case", {
  result <- fisher_mle(
    par = c(rate = 2.0),
    vcov = matrix(0.04, 1, 1),
    loglik_val = -30,
    nobs = 50
  )

  samp_fn <- sampler(result)
  set.seed(6007)
  samples <- samp_fn(500)
  expect_equal(ncol(samples), 1)
  expect_equal(nrow(samples), 500)
  expect_equal(mean(samples), 2.0, tolerance = 0.1)
})

test_that("sampler.fisher_mle errors when vcov is NULL", {
  result <- fisher_mle(
    par = c(rate = 2.0),
    loglik_val = -30,
    nobs = 50
  )
  expect_error(sampler(result), "vcov is NULL")
})


# =============================================================================
# FISHERIAN INFERENCE COVERAGE TESTS
# =============================================================================

test_that("likelihood_interval works for single parameter model", {
  set.seed(6009)
  df <- data.frame(x = rexp(100, rate = 2))
  model <- exponential_lifetime("x")
  result <- fit(model)(df, par = c(rate = 1))

  li <- likelihood_interval(result, data = df, model = model, k = 8)
  expect_true(!is.null(li))
})

test_that("print.likelihood_interval outputs expected text", {
  set.seed(6010)
  df <- data.frame(x = rexp(100, rate = 2))
  model <- exponential_lifetime("x")
  result <- fit(model)(df, par = c(rate = 1))

  li <- likelihood_interval(result, data = df, model = model, k = 8)
  out <- capture.output(print(li))
  expect_true(any(grepl("Likelihood Interval", out)))
})

test_that("profile_loglik errors for > 2 parameters", {
  result <- fisher_mle(
    par = c(a = 1, b = 2, c = 3),
    vcov = diag(3) * 0.01,
    loglik_val = -50,
    nobs = 100
  )
  # Use a mock model for this error test
  mock_model <- structure(list(), class = c("mock_3param", "likelihood_model"))
  loglik.mock_3param <<- function(model, ...) function(df, par, ...) 0
  df <- data.frame(x = 1:10)
  expect_error(
    profile_loglik(result, data = df, model = mock_model, param = 1:3),
    "1 or 2 parameters"
  )
  rm(loglik.mock_3param, envir = .GlobalEnv)
})

test_that("print.profile_loglik outputs expected text", {
  set.seed(6013)
  df <- data.frame(t = rexp(50, rate = 2))
  model <- exponential_lifetime("t")
  result <- fit(model)(df)

  pl <- profile_loglik(result, data = df, model = model,
                       param = 1, n_grid = 10)
  out <- capture.output(print(pl))
  expect_true(any(grepl("Profile Log-Likelihood", out)))
  expect_true(any(grepl("Grid points:", out)))
})

test_that("deviance.fisher_mle errors on non-fisher_mle null model", {
  result <- fisher_mle(
    par = c(a = 1),
    vcov = matrix(0.1, 1, 1),
    loglik_val = -50,
    nobs = 100
  )
  expect_error(deviance(result, null_model = list(loglik = -60)),
               "fisher_mle")
})


# =============================================================================
# FIM AND MOCK MODEL COVERAGE TESTS
# =============================================================================

# Helper: register mock exponential loglik/rdata methods in the global env.
# Returns a model object that falls through to fim.likelihood_model (unlike
# exponential_lifetime, which has its own fim method).
make_mock_exp <- function(class_name) {
  assign(paste0("loglik.", class_name), function(model, ...) {
    function(df, par, ...) {
      lambda <- par[1]
      if (lambda <= 0) return(-.Machine$double.xmax / 2)
      nrow(df) * log(lambda) - lambda * sum(df$x)
    }
  }, envir = .GlobalEnv)

  assign(paste0("rdata.", class_name), function(model, ...) {
    function(theta, n, ...) {
      data.frame(x = rexp(n, rate = theta[1]))
    }
  }, envir = .GlobalEnv)

  structure(list(), class = c(class_name, "likelihood_model"))
}

cleanup_mock_exp <- function(class_name) {
  rm(list = paste0(c("loglik.", "rdata."), class_name), envir = .GlobalEnv)
}

test_that("fim.likelihood_model computes FIM via Monte Carlo (neg Hessian)", {
  mock_model <- make_mock_exp("mock_exp")

  set.seed(6014)
  fim_fn <- fim(mock_model)
  # True FIM for exp(rate=2) with n=50 is n/rate^2 = 50/4 = 12.5
  fim_mat <- fim_fn(theta = c(2), n_obs = 50, n_samples = 500)
  expect_true(is.matrix(fim_mat))
  expect_equal(nrow(fim_mat), 1)
  expect_equal(ncol(fim_mat), 1)
  expect_equal(fim_mat[1, 1], 12.5, tolerance = 2)

  cleanup_mock_exp("mock_exp")
})

test_that("fim.likelihood_model scales linearly with n_obs", {
  mock_model <- make_mock_exp("mock_exp2")

  set.seed(7001)
  fim_fn <- fim(mock_model)
  fim_10 <- fim_fn(theta = c(2), n_obs = 10, n_samples = 2000)
  fim_50 <- fim_fn(theta = c(2), n_obs = 50, n_samples = 2000)

  # FIM should scale as n_obs / lambda^2, so ratio ~ 50/10 = 5
  expect_equal(fim_50[1, 1] / fim_10[1, 1], 5, tolerance = 0.5)

  cleanup_mock_exp("mock_exp2")
})

test_that("fim.likelihood_model adds parameter names from theta", {
  mock_model <- make_mock_exp("mock_exp3")

  set.seed(7002)
  fim_fn <- fim(mock_model)
  fim_mat <- fim_fn(theta = c(rate = 2), n_obs = 10, n_samples = 100)
  expect_equal(rownames(fim_mat), "rate")
  expect_equal(colnames(fim_mat), "rate")

  cleanup_mock_exp("mock_exp3")
})


# =============================================================================
# BIAS.FISHER_MLE TESTS
# =============================================================================

test_that("bias.fisher_mle returns zeros without model", {
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    vcov = diag(2) * 0.01,
    loglik_val = -50,
    nobs = 100
  )
  b <- bias(result)
  expect_equal(b, c(0, 0))

  # With theta but no model, still zeros
  b2 <- bias(result, theta = c(1, 2))
  expect_equal(b2, c(0, 0))
})

test_that("bias.fisher_mle estimates bias via MC when model provided", {
  model <- exponential_lifetime("t")
  set.seed(7020)
  df <- data.frame(t = rexp(200, rate = 2))
  result <- fit(model)(df)

  # MC bias for exponential MLE should be near zero for large n
  b <- bias(result, theta = c(lambda = 2), model = model, n_sim = 200)
  expect_equal(length(b), 1)
  # Bias should be small relative to the parameter value
  expect_true(abs(b[1]) < 0.5)
})


# =============================================================================
# FIM ... FORWARDING + NaN GUARD TESTS
# =============================================================================

test_that("fim.likelihood_model forwards ... to hess_fn", {
  mock_model <- structure(list(), class = c("mock_dots_fim", "likelihood_model"))

  loglik.mock_dots_fim <<- function(model, ...) {
    function(df, par, multiplier = 1, ...) {
      multiplier * (nrow(df) * log(par[1]) - par[1] * sum(df$x))
    }
  }

  rdata.mock_dots_fim <<- function(model, ...) {
    function(theta, n, ...) {
      data.frame(x = rexp(n, rate = theta[1]))
    }
  }

  set.seed(8002)
  fim_fn <- fim(mock_model)
  # With multiplier=1 (default), FIM for exp(rate=2), n=10 is 10/4 = 2.5
  fim_default <- fim_fn(theta = c(2), n_obs = 10, n_samples = 500)
  # With multiplier=2, FIM should be ~2x
  fim_doubled <- fim_fn(theta = c(2), n_obs = 10, n_samples = 500, multiplier = 2)

  expect_equal(fim_doubled[1, 1] / fim_default[1, 1], 2, tolerance = 0.5)

  rm(loglik.mock_dots_fim, rdata.mock_dots_fim, envir = .GlobalEnv)
})

test_that("fim.likelihood_model handles NaN from hess_loglik gracefully", {
  mock_model <- structure(list(), class = c("mock_nan_fim", "likelihood_model"))

  call_count <- 0L
  hess_loglik.mock_nan_fim <<- function(model, ...) {
    function(df, par, ...) {
      # Return NaN matrix every 3rd call
      call_count <<- call_count + 1L
      if (call_count %% 3 == 0) return(matrix(NaN, 1, 1))
      matrix(-nrow(df) / par[1]^2, 1, 1)
    }
  }

  rdata.mock_nan_fim <<- function(model, ...) {
    function(theta, n, ...) {
      data.frame(x = rexp(n, rate = theta[1]))
    }
  }

  set.seed(8003)
  fim_fn <- fim(mock_model)
  # Should still produce a valid FIM despite some NaN samples
  expect_warning(
    fim_mat <- fim_fn(theta = c(2), n_obs = 10, n_samples = 100),
    "samples failed"
  )
  expect_true(is.matrix(fim_mat))
  expect_true(all(is.finite(fim_mat)))

  rm(hess_loglik.mock_nan_fim, rdata.mock_nan_fim, envir = .GlobalEnv)
  call_count <- NULL
})

# =============================================================================
# MSE TESTS
# =============================================================================

test_that("mse.fisher_mle returns Vcov when bias is zero (asymptotic)", {
  model <- exponential_lifetime("t")
  set.seed(8005)
  df <- data.frame(t = rexp(100, rate = 2))
  result <- fit(model)(df)

  # Without MC, bias = 0, so MSE = Vcov
  mse_val <- mse(result, theta = c(lambda = 2))
  expect_equal(mse_val, vcov(result))
})

test_that("mse.fisher_mle forwards model and n_sim to bias for MC-based MSE", {
  model <- exponential_lifetime("t")
  set.seed(8005)
  df <- data.frame(t = rexp(100, rate = 2))
  result <- fit(model)(df)

  # With model, MSE accounts for finite-sample bias
  mse_mc <- mse(result, theta = c(lambda = 2), model = model, n_sim = 100)
  expect_true(is.finite(mse_mc))
})

# =============================================================================
# SAMPLER ... SCOPING TEST
# =============================================================================

test_that("sampler.likelihood_model bootstrap works correctly", {
  model <- exponential_lifetime("t")
  set.seed(8006)
  df <- data.frame(t = rexp(50, rate = 2))

  boot_sampler <- sampler(model, df, par = c(lambda = 2))
  boot_result <- boot_sampler(50)

  expect_s3_class(boot_result, "fisher_boot")
  expect_equal(ncol(boot_result$replicates), 1)
  expect_equal(nrow(boot_result$replicates), 50)
})

# =============================================================================
# DEFAULT METHOD ERROR MESSAGES
# =============================================================================

test_that("loglik.likelihood_model gives clear error for unimplemented class", {
  fake_model <- structure(list(), class = c("fake_model", "likelihood_model"))
  expect_error(loglik(fake_model), "must implement loglik")
})

test_that("rdata.likelihood_model gives clear error", {
  fake_model <- structure(list(), class = c("fake_model2", "likelihood_model"))
  loglik.fake_model2 <<- function(model, ...) function(df, par, ...) 0
  expect_error(rdata(fake_model), "does not implement rdata")
  rm(loglik.fake_model2, envir = .GlobalEnv)
})

# =============================================================================
# LRT RESULT CLASS + PRINT
# =============================================================================

test_that("lrt returns lrt_result class with print method", {
  set.seed(8007)
  df <- data.frame(t = rexp(50, rate = 2))
  model <- exponential_lifetime("t")

  result <- lrt(model, model, df,
                null_par = c(2), alt_par = c(2), dof = 0)

  expect_s3_class(result, "lrt_result")
  expect_true("stat" %in% names(result))
  expect_true("p.value" %in% names(result))

  out <- capture.output(print(result))
  expect_true(any(grepl("Likelihood Ratio Test", out)))
})

# =============================================================================
# EVIDENCE IS A GENERIC
# =============================================================================

test_that("evidence dispatches via S3", {
  set.seed(8008)
  df <- data.frame(t = rexp(50, rate = 2))
  model <- exponential_lifetime("t")

  e <- evidence(model, data = df, theta1 = c(lambda = 2), theta2 = c(lambda = 1))
  expect_true(is.numeric(e))
  expect_equal(length(e), 1)

  # Verify it matches manual computation
  ll <- loglik(model)
  expected <- ll(df, c(lambda = 2)) - ll(df, c(lambda = 1))
  expect_equal(e, expected)
})


# =============================================================================
# COVERAGE RECOVERY TESTS
# =============================================================================
# These tests target multivariate code paths that require a 2+ parameter model.
#
# Mock normal model: loglik.mock_norm, score.mock_norm, hess_loglik.mock_norm
# are registered in the global environment and cleaned up after each test
# block that needs them.

# Helper: set up and tear down mock_norm model
make_mock_norm <- function() {
  model <- structure(list(ob_col = "x"), class = c("mock_norm", "likelihood_model"))

  loglik.mock_norm <<- function(model, ...) {
    function(df, par, ...) {
      mu <- par[1]; sigma <- par[2]
      if (sigma <= 0) return(-.Machine$double.xmax / 2)
      sum(dnorm(df[[model$ob_col]], mu, sigma, log = TRUE))
    }
  }

  model
}

cleanup_mock_norm <- function() {
  to_rm <- intersect(c("loglik.mock_norm"), ls(envir = .GlobalEnv))
  if (length(to_rm)) rm(list = to_rm, envir = .GlobalEnv)
}

# ---------------------------------------------------------------------------
# sampler.fisher_mle: multivariate branch (mvtnorm::rmvnorm) — line 543
# ---------------------------------------------------------------------------

test_that("sampler.fisher_mle uses mvtnorm for 2-param model", {
  skip_if_not_installed("mvtnorm")

  result <- fisher_mle(
    par = c(mu = 5.0, sigma = 2.0),
    vcov = matrix(c(0.04, 0.0, 0.0, 0.02), 2, 2),
    loglik_val = -150,
    nobs = 100
  )

  samp_fn <- sampler(result)
  set.seed(9001)
  samples <- samp_fn(500)

  expect_equal(dim(samples), c(500, 2))
  expect_equal(mean(samples[, 1]), 5.0, tolerance = 0.1)
  expect_equal(mean(samples[, 2]), 2.0, tolerance = 0.1)
})

# ---------------------------------------------------------------------------
# likelihood_interval: multivariate (profile) branch — lines 149-213
# ---------------------------------------------------------------------------

test_that("likelihood_interval computes profile interval for 2-param model", {
  skip_if_not_installed("mvtnorm")
  model <- make_mock_norm()
  on.exit(cleanup_mock_norm())

  set.seed(9002)
  df <- data.frame(x = rnorm(100, mean = 3, sd = 1.5))

  # Fit via fit.likelihood_model (optim-based)
  solver <- fit(model)
  result <- solver(df, par = c(mu = 0, sigma = 1))

  # Compute profile likelihood interval for mu (param = 1)
  li <- likelihood_interval(result, data = df, model = model, k = 8, param = 1)

  expect_true(inherits(li, "likelihood_interval"))
  expect_equal(nrow(li), 1)
  expect_equal(ncol(li), 2)

  # Interval should contain the MLE for mu
  mu_hat <- coef(result)[1]
  expect_true(li[1, 1] < mu_hat && mu_hat < li[1, 2])
})

test_that("likelihood_interval computes profile intervals for both params", {
  skip_if_not_installed("mvtnorm")
  model <- make_mock_norm()
  on.exit(cleanup_mock_norm())

  set.seed(9003)
  df <- data.frame(x = rnorm(80, mean = 2, sd = 1))

  solver <- fit(model)
  result <- solver(df, par = c(mu = 0, sigma = 1))

  # Profile interval for all params (NULL => both)
  li <- likelihood_interval(result, data = df, model = model, k = 8)

  expect_equal(nrow(li), 2)
  mu_hat <- coef(result)[1]
  sigma_hat <- coef(result)[2]
  expect_true(li[1, 1] < mu_hat && mu_hat < li[1, 2])
  expect_true(li[2, 1] < sigma_hat && sigma_hat < li[2, 2])
})

test_that("likelihood_interval handles character param name for 2-param model", {
  skip_if_not_installed("mvtnorm")
  model <- make_mock_norm()
  on.exit(cleanup_mock_norm())

  set.seed(9004)
  df <- data.frame(x = rnorm(60, mean = 1, sd = 2))

  solver <- fit(model)
  result <- solver(df, par = c(mu = 0, sigma = 1))

  # Use character param name
  li <- likelihood_interval(result, data = df, model = model, k = 8,
                            param = "mu")
  expect_equal(nrow(li), 1)
  mu_hat <- coef(result)["mu"]
  expect_true(li[1, 1] < mu_hat && mu_hat < li[1, 2])
})

# ---------------------------------------------------------------------------
# profile_loglik: 1D profile over a 2-param model — lines 310-329
# ---------------------------------------------------------------------------

test_that("profile_loglik 1D profile over nuisance parameter works", {
  skip_if_not_installed("mvtnorm")
  model <- make_mock_norm()
  on.exit(cleanup_mock_norm())

  set.seed(9005)
  df <- data.frame(x = rnorm(80, mean = 3, sd = 1.5))

  solver <- fit(model)
  result <- solver(df, par = c(mu = 0, sigma = 1))

  # Profile over mu (param=1), optimising out sigma
  prof <- profile_loglik(result, data = df, model = model, param = 1, n_grid = 15)

  expect_true(inherits(prof, "profile_loglik"))
  expect_equal(nrow(prof), 15)
  expect_true("loglik" %in% names(prof))
  expect_true("support" %in% names(prof))
  expect_true("relative_likelihood" %in% names(prof))

  # Profile should peak near MLE
  mu_hat <- unname(coef(result)[1])
  max_idx <- which.max(prof$loglik)
  expect_equal(unname(prof[[1]][max_idx]), mu_hat, tolerance = 0.5)
})

# ---------------------------------------------------------------------------
# profile_loglik: 2D profile — lines 332-358
# ---------------------------------------------------------------------------

test_that("profile_loglik 2D profile produces correct shape data frame", {
  skip_if_not_installed("mvtnorm")
  model <- make_mock_norm()
  on.exit(cleanup_mock_norm())

  set.seed(9006)
  df <- data.frame(x = rnorm(60, mean = 2, sd = 1))

  solver <- fit(model)
  result <- solver(df, par = c(mu = 0, sigma = 1))

  prof <- profile_loglik(result, data = df, model = model,
                         param = 1:2, n_grid = 8)

  expect_true(inherits(prof, "profile_loglik"))
  expect_equal(nrow(prof), 8 * 8)   # expand.grid of 8x8 grid
  expect_true("loglik" %in% names(prof))
  expect_true("support" %in% names(prof))
  expect_true("relative_likelihood" %in% names(prof))
})

# ---------------------------------------------------------------------------
# print.profile_loglik: 2D case — covered when param has length 2
# ---------------------------------------------------------------------------

test_that("print.profile_loglik prints correctly for 2D profile", {
  skip_if_not_installed("mvtnorm")
  model <- make_mock_norm()
  on.exit(cleanup_mock_norm())

  set.seed(9007)
  df <- data.frame(x = rnorm(40, mean = 1, sd = 1))

  solver <- fit(model)
  result <- solver(df, par = c(mu = 0, sigma = 1))

  prof <- profile_loglik(result, data = df, model = model,
                         param = 1:2, n_grid = 5)

  out <- capture.output(print(prof))
  expect_true(any(grepl("Profile Log-Likelihood", out)))
  expect_true(any(grepl("Grid points:", out)))
  # Should mention both parameter names
  expect_true(any(grepl("mu", out)))
})

# ---------------------------------------------------------------------------
# bias.fisher_mle: nobs NULL error — line 294
# ---------------------------------------------------------------------------

test_that("bias.fisher_mle errors when nobs is NULL and model is provided", {
  model <- exponential_lifetime("t")

  result <- fisher_mle(
    par = c(lambda = 2.0),
    vcov = matrix(0.04, 1, 1),
    loglik_val = -50
    # nobs intentionally omitted (NULL)
  )

  expect_error(
    bias(result, theta = c(lambda = 2), model = model),
    "nobs not available"
  )
})

# ---------------------------------------------------------------------------
# bias.fisher_mle: warning when most MC replicates fail — lines 314-317
# ---------------------------------------------------------------------------

test_that("bias.fisher_mle warns when most MC replicates fail", {
  # Create a model whose fit() always errors
  bad_model <- structure(
    list(ob_col = "t"),
    class = c("mock_bad_fit", "likelihood_model")
  )

  loglik.mock_bad_fit <<- function(model, ...) {
    function(df, par, ...) {
      lambda <- par[1]
      if (lambda <= 0) return(-.Machine$double.xmax / 2)
      nrow(df) * log(lambda) - lambda * sum(df[[model$ob_col]])
    }
  }

  rdata.mock_bad_fit <<- function(model, ...) {
    function(theta, n, ...) {
      data.frame(t = rexp(n, rate = theta[1]))
    }
  }

  # Override fit to always error
  fit.mock_bad_fit <<- function(object, ...) {
    function(df, par = NULL, ...) {
      stop("Intentional fit failure for testing")
    }
  }

  on.exit({
    rm(list = intersect(
      c("loglik.mock_bad_fit", "rdata.mock_bad_fit", "fit.mock_bad_fit"),
      ls(envir = .GlobalEnv)
    ), envir = .GlobalEnv)
  })

  set.seed(9008)
  result <- fisher_mle(
    par = c(lambda = 2.0),
    vcov = matrix(0.04, 1, 1),
    loglik_val = -50,
    nobs = 50
  )

  # All replicates will fail, so we get a warning + NA return
  expect_warning(
    b <- bias(result, theta = c(lambda = 2), model = bad_model, n_sim = 20),
    "MC replicates succeeded"
  )
  expect_true(all(is.na(b)))
})

# ---------------------------------------------------------------------------
# print.fisher_boot: bias line is covered when boot object has 2+ params
# ---------------------------------------------------------------------------

test_that("print.fisher_boot covers format(bias(x)) for multivariate bootstrap", {
  skip_if_not_installed("boot")

  # Build a 2-param mock model that fit.likelihood_model can use
  model <- make_mock_norm()
  on.exit(cleanup_mock_norm())

  set.seed(9009)
  df <- data.frame(x = rnorm(50, mean = 2, sd = 1))

  boot_sampler <- sampler(model, df, par = c(mu = 0, sigma = 1))
  boot_obj <- boot_sampler(30)

  out <- capture.output(print(boot_obj))
  expect_true(any(grepl("Bootstrap MLE", out)))
  expect_true(any(grepl("Bootstrap bias:", out)))
  # Should show 2 bias values
  expect_equal(length(bias(boot_obj)), 2)
})

# ---------------------------------------------------------------------------
# core-generics.R line 247: "All Monte Carlo samples failed" error in fim
# ---------------------------------------------------------------------------

test_that("fim.likelihood_model errors when all MC samples fail", {
  always_fail_model <- structure(
    list(),
    class = c("mock_always_fail_fim", "likelihood_model")
  )

  loglik.mock_always_fail_fim <<- function(model, ...) {
    function(df, par, ...) 0
  }

  hess_loglik.mock_always_fail_fim <<- function(model, ...) {
    function(df, par, ...) stop("Always fails")
  }

  rdata.mock_always_fail_fim <<- function(model, ...) {
    function(theta, n, ...) data.frame(x = 1)
  }

  on.exit({
    rm(list = intersect(
      c("loglik.mock_always_fail_fim",
        "hess_loglik.mock_always_fail_fim",
        "rdata.mock_always_fail_fim"),
      ls(envir = .GlobalEnv)
    ), envir = .GlobalEnv)
  })

  fim_fn <- fim(always_fail_model)
  expect_error(
    fim_fn(theta = c(1), n_obs = 10, n_samples = 5),
    "All Monte Carlo samples failed"
  )
})

# ---------------------------------------------------------------------------
# exponential_lifetime loglik: negative lambda guard — line 76
# ---------------------------------------------------------------------------

test_that("loglik.exponential_lifetime returns sentinel for non-positive lambda", {
  model <- exponential_lifetime("t")
  ll_fn <- loglik(model)
  df <- data.frame(t = c(1, 2, 3))

  # Negative lambda
  val_neg <- ll_fn(df, c(-1))
  expect_equal(val_neg, -.Machine$double.xmax / 2)

  # Zero lambda
  val_zero <- ll_fn(df, c(0))
  expect_equal(val_zero, -.Machine$double.xmax / 2)
})

# ---------------------------------------------------------------------------
# fit.likelihood_model: auto-switch 1-param to BFGS — line 61
# (the SANN test uses a custom model; this test exercises the 1-param
#  Nelder-Mead -> BFGS switch via fit.likelihood_model directly)
# ---------------------------------------------------------------------------

test_that("fit.likelihood_model auto-switches 1-param Nelder-Mead to BFGS", {
  # Use make_mock_norm but override to 1-param normal (fixed sigma)
  mock_1p <- structure(list(), class = c("mock_1param_bfgs", "likelihood_model"))

  loglik.mock_1param_bfgs <<- function(model, ...) {
    function(df, par, ...) {
      mu <- par[1]
      sum(dnorm(df$x, mu, sd = 1, log = TRUE))
    }
  }

  on.exit(rm(list = intersect("loglik.mock_1param_bfgs", ls(.GlobalEnv)),
             envir = .GlobalEnv))

  set.seed(9010)
  df <- data.frame(x = rnorm(50, mean = 3, sd = 1))

  # Deliberately pass method = "Nelder-Mead" with 1 param; should switch to BFGS
  solver <- fit(mock_1p)
  result <- solver(df, par = c(mu = 0), method = "Nelder-Mead")

  expect_true(inherits(result, "fisher_mle"))
  expect_equal(unname(coef(result)[1]), 3, tolerance = 0.3)
})
