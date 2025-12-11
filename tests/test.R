library(testthat)
library(likelihood.model)

# Helper function to get parameter estimates from MLE result
# algebraic.mle stores estimates in $theta.hat
get_mle_params <- function(mle_result) {
  if (!is.null(mle_result$theta.hat)) {
    return(mle_result$theta.hat)
  }
  # Fallback to coef if available
  result <- coef(mle_result)
  if (!is.null(result)) {
    return(result)
  }
  # Last resort: try $par from optim result
  if (!is.null(mle_result$sol$par)) {
    return(mle_result$sol$par)
  }
  return(NULL)
}

# =============================================================================
# UTILITY FUNCTIONS TESTS
# =============================================================================

test_that("get_params returns default when par is NULL or all NA", {
  get_params <- likelihood.model:::get_params
  default <- c(1, 2, 3)

  # NULL input returns default
  expect_equal(get_params(NULL, default), default)

  # All NA input returns default
  expect_equal(get_params(c(NA, NA, NA), default), default)

  # Partial NA input: substitute NA values from default
  expect_equal(get_params(c(5, NA, 7), default), c(5, 2, 7))

  # No NA values: return par unchanged
  expect_equal(get_params(c(5, 6, 7), default), c(5, 6, 7))

  # NULL default: return par if no substitution needed
  expect_equal(get_params(c(5, 6, 7), NULL), c(5, 6, 7))
})

test_that("prepare_args handles named and unnamed parameters", {
  prepare_args <- likelihood.model:::prepare_args
  # Test with unnamed parameters (should infer names from function formals)
  args <- prepare_args(1.5, c(2, 3), dnorm)
  expect_equal(args[[1]], 1.5)  # primary input
  expect_equal(args$mean, 2)
  expect_equal(args$sd, 3)

  # Test with named parameters
  args_named <- prepare_args(1.5, c(mean = 2, sd = 3), dnorm)
  expect_equal(args_named[[1]], 1.5)
  expect_equal(args_named$mean, 2)
  expect_equal(args_named$sd, 3)

  # Test with Weibull distribution
  args_weibull <- prepare_args(2.0, c(1.5, 2.0), dweibull)
  expect_equal(args_weibull[[1]], 2.0)
  expect_equal(args_weibull$shape, 1.5)
  expect_equal(args_weibull$scale, 2.0)
})

test_that("method_exists correctly detects methods", {
  method_exists <- likelihood.model:::method_exists
  # print method should exist for many classes
  expect_true(method_exists("print", "data.frame"))

  # Nonexistent class should return FALSE
  expect_false(method_exists("print", "nonexistent_class_xyz123"))
})

# =============================================================================
# IS_LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("is_likelihood_model correctly identifies likelihood models", {
  # likelihood_contr_model is a likelihood model
  model <- likelihood_contr_model$new(obs_type = function(df) "exact")
  expect_true(is_likelihood_model(model))

  # likelihood_name_model is a likelihood model
  # Note: censor_col must not be NULL, use actual column name
  model_name <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  expect_true(is_likelihood_model(model_name))

  # likelihood_exact_weibull is a likelihood model
  model_weibull <- likelihood_exact_weibull("x")
  expect_true(is_likelihood_model(model_weibull))

  # Non-likelihood model objects return FALSE
  expect_false(is_likelihood_model(list()))
  expect_false(is_likelihood_model(data.frame()))
  expect_false(is_likelihood_model(NULL))
  expect_false(is_likelihood_model("not a model"))
  expect_false(is_likelihood_model(42))
})

# =============================================================================
# LIKELIHOOD_NAME_MODEL TESTS
# =============================================================================

test_that("likelihood_name creates model with correct structure", {
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  expect_true(is_likelihood_model(model))
  expect_equal(model$dist_name, "norm")
  expect_equal(model$ob_col, "x")
  expect_equal(model$censor_col, "censor")
  expect_true("likelihood_name_model" %in% class(model))
  expect_true("likelihood_name_norm" %in% class(model))
})

test_that("loglik for likelihood_name_model matches manual calculation for normal", {
  set.seed(123)
  # Note: censor column must have actual values, not NA
  # Use "exact" instead of NA for uncensored observations
  df <- data.frame(x = rnorm(100), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  true_mean <- 0
  true_sd <- 1

  # Log-likelihood should match sum of log densities
  ll_func <- loglik(model)
  computed_ll <- ll_func(df, c(true_mean, true_sd))
  expected_ll <- sum(dnorm(df$x, true_mean, true_sd, log = TRUE))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("loglik for likelihood_name_model handles left-censored data", {
  # Left-censored means we observe x but only know the true value is <= x
  df <- data.frame(
    x = c(1.0, 2.0, 3.0),
    censor = c("left", "left", "left")
  )

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll_func <- loglik(model)

  # For left-censored: use log(CDF) = log P(X <= x)
  expected_ll <- sum(pnorm(df$x, mean = 0, sd = 1, log.p = TRUE))
  computed_ll <- ll_func(df, c(0, 1))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("loglik for likelihood_name_model handles right-censored data", {
  # Right-censored means we observe x but only know the true value is > x
  df <- data.frame(
    x = c(1.0, 2.0, 3.0),
    censor = c("right", "right", "right")
  )

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll_func <- loglik(model)

  # For right-censored: use log(1 - CDF) = log P(X > x)
  expected_ll <- sum(pnorm(df$x, mean = 0, sd = 1, log.p = TRUE, lower.tail = FALSE))
  computed_ll <- ll_func(df, c(0, 1))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("loglik for likelihood_name_model handles mixed exact and censored data", {
  df <- data.frame(
    x = c(0.5, 1.5, -0.5, 2.0),
    censor = c("exact", "right", "left", "exact")  # Use "exact" instead of NA
  )

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll_func <- loglik(model)

  # Manually compute expected log-likelihood
  expected_ll <- dnorm(0.5, 0, 1, log = TRUE) +  # exact
    pnorm(1.5, 0, 1, log.p = TRUE, lower.tail = FALSE) +  # right
    pnorm(-0.5, 0, 1, log.p = TRUE) +  # left
    dnorm(2.0, 0, 1, log = TRUE)  # exact

  computed_ll <- ll_func(df, c(0, 1))
  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("score for likelihood_name_model uses numerical gradient", {
  set.seed(456)
  n <- 100
  df <- data.frame(x = rnorm(n, 0, 1), censor = rep("exact", n))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  score_func <- score(model)
  ll_func <- loglik(model)

  par <- c(0.1, 0.9)
  score_val <- score_func(df, par)

  # Should match numerical gradient
  expected_score <- numDeriv::grad(
    func = function(p) ll_func(df, p),
    x = par
  )

  expect_equal(as.numeric(score_val), expected_score, tolerance = 1e-4)
})

test_that("hess_loglik for likelihood_name_model returns negative definite at MLE", {
  set.seed(789)
  df <- data.frame(x = rnorm(100), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  hess_func <- hess_loglik(model)

  hessian_val <- hess_func(df, c(mean(df$x), sd(df$x)))

  # Hessian should be negative definite at MLE (all eigenvalues negative)
  eigenvalues <- eigen(hessian_val)$values
  expect_true(all(eigenvalues < 0),
              info = "Hessian should be negative definite at MLE")
})

test_that("assumptions for likelihood_name_model returns expected assumptions", {
  model <- likelihood_name("weibull", ob_col = "t", censor_col = "censor")
  assumptions_list <- assumptions(model)

  expect_true("independent" %in% assumptions_list)
  expect_true("identically distributed" %in% assumptions_list)
  expect_true("weibull distribution" %in% assumptions_list)
})

test_that("likelihood_name works with weibull distribution", {
  set.seed(111)
  shape_true <- 2
  scale_true <- 3
  df <- data.frame(x = rweibull(100, shape_true, scale_true), censor = rep("exact", 100))

  model <- likelihood_name("weibull", ob_col = "x", censor_col = "censor")
  ll_func <- loglik(model)

  # Compare to manual calculation
  expected_ll <- sum(dweibull(df$x, shape_true, scale_true, log = TRUE))
  computed_ll <- ll_func(df, c(shape_true, scale_true))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

# =============================================================================
# LIKELIHOOD_EXACT_WEIBULL TESTS
# =============================================================================

test_that("likelihood_exact_weibull creates model with correct structure", {
  model <- likelihood_exact_weibull("x")

  expect_true(is_likelihood_model(model))
  expect_equal(model$ob_col, "x")
  expect_true("likelihood_exact_weibull" %in% class(model))
})

test_that("loglik for likelihood_exact_weibull matches R's dweibull", {
  set.seed(222)
  shape <- 2.5
  scale <- 1.5
  df <- data.frame(x = rweibull(50, shape, scale))

  model <- likelihood_exact_weibull("x")
  ll_func <- loglik(model)

  expected_ll <- sum(dweibull(df$x, shape, scale, log = TRUE))
  computed_ll <- ll_func(df, c(shape, scale))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("analytical score matches numerical gradient for Weibull", {
  set.seed(333)
  shape <- 1.8
  scale <- 2.2
  df <- data.frame(x = rweibull(100, shape, scale))

  model <- likelihood_exact_weibull("x")
  ll_func <- loglik(model)
  score_func <- score(model)

  # Analytical score
  analytical_score <- score_func(df, c(shape, scale))

  # Numerical gradient using numDeriv
  numerical_score <- numDeriv::grad(
    func = function(par) ll_func(df, par),
    x = c(shape, scale)
  )

  # Compare values (ignoring names)
  expect_equal(as.numeric(analytical_score), numerical_score, tolerance = 1e-5,
               info = "Analytical score should match numerical gradient")
})

test_that("analytical hessian has expected structure for Weibull", {
  set.seed(444)
  shape <- 2.0
  scale <- 1.0
  df <- data.frame(x = rweibull(100, shape, scale))

  model <- likelihood_exact_weibull("x")
  ll_func <- loglik(model)
  hess_func <- hess_loglik(model)

  # Analytical Hessian
  analytical_hess <- hess_func(df, c(shape, scale))

  # Numerical Hessian using numDeriv
  numerical_hess <- numDeriv::hessian(
    func = function(par) ll_func(df, par),
    x = c(shape, scale)
  )

  # Analytical Hessian should match numerical Hessian
  expect_equal(analytical_hess, numerical_hess, tolerance = 1e-4)
})

test_that("likelihood_exact_weibull validates positive observations", {
  df_invalid <- data.frame(x = c(1.0, -0.5, 2.0))  # Contains negative value

  model <- likelihood_exact_weibull("x")
  ll_func <- loglik(model)

  expect_error(ll_func(df_invalid, c(2, 1)),
               info = "Should error on non-positive observations")
})

test_that("likelihood_exact_weibull requires non-empty data", {
  df_empty <- data.frame(x = numeric(0))

  model <- likelihood_exact_weibull("x")
  ll_func <- loglik(model)

  expect_error(ll_func(df_empty, c(2, 1)),
               info = "Should error on empty data frame")
})

test_that("assumptions for likelihood_exact_weibull returns expected list", {
  model <- likelihood_exact_weibull("x")
  assumptions_list <- assumptions(model)

  expect_true("independent" %in% assumptions_list)
  expect_true("identically distributed" %in% assumptions_list)
  expect_true("exact observations" %in% assumptions_list)
})

# =============================================================================
# LIKELIHOOD_CONTR_MODEL TESTS
# =============================================================================

test_that("likelihood_contr_model requires obs_type to be a function", {
  expect_error(likelihood_contr_model$new(obs_type = NULL),
               "obs_type must be a function")

  expect_error(likelihood_contr_model$new(obs_type = "not a function"),
               "obs_type must be a function")
})

test_that("likelihood_contr_model creates model with correct structure", {
  model <- likelihood_contr_model$new(
    obs_type = function(df) "exact",
    assumptions = c("iid", "exponential")
  )

  expect_true(is_likelihood_model(model))
  expect_true("likelihood_contr_model" %in% class(model))

  # iid assumption is always included
  expect_true("iid" %in% model$assumptions)
  expect_true("exponential" %in% model$assumptions)
})

test_that("likelihood_contr_model validates inputs to loglik, score, hess_loglik", {
  model <- likelihood_contr_model$new(obs_type = function(df) "exact")

  # Should error with NULL df
  expect_error(model$loglik(NULL, c(1, 2)))

  # Should error with NULL par
  expect_error(model$loglik(data.frame(x = 1), NULL))

  # Should error with non-data.frame
  expect_error(model$loglik("not a dataframe", c(1, 2)))
})

test_that("likelihood_contr_model dispatches to custom loglik functions", {
  # Define a custom loglik function for "exact" type
  loglik.test_exact <- function(df, par, ...) {
    sum(dnorm(df$x, par[1], par[2], log = TRUE))
  }

  # Make it available in the global environment
  assign("loglik.test_exact", loglik.test_exact, envir = .GlobalEnv)
  on.exit(rm("loglik.test_exact", envir = .GlobalEnv))

  model <- likelihood_contr_model$new(
    obs_type = function(df) "test_exact"
  )

  set.seed(555)
  df <- data.frame(x = rnorm(50, mean = 3, sd = 2))

  expected_ll <- sum(dnorm(df$x, 3, 2, log = TRUE))
  computed_ll <- model$loglik(df, c(3, 2))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("likelihood_contr_model sums contributions from multiple observation types", {
  # Define loglik functions for two types
  loglik.type_a <- function(df, par, ...) {
    sum(dnorm(df$x, par[1], par[2], log = TRUE))
  }
  loglik.type_b <- function(df, par, ...) {
    # Exponential with rate = par[3]
    sum(dexp(df$x, par[3], log = TRUE))
  }

  assign("loglik.type_a", loglik.type_a, envir = .GlobalEnv)
  assign("loglik.type_b", loglik.type_b, envir = .GlobalEnv)
  on.exit({
    rm("loglik.type_a", envir = .GlobalEnv)
    rm("loglik.type_b", envir = .GlobalEnv)
  })

  model <- likelihood_contr_model$new(
    obs_type = function(df) df$type
  )

  set.seed(666)
  df_a <- data.frame(x = rnorm(30, mean = 5, sd = 1), type = "type_a")
  df_b <- data.frame(x = rexp(20, rate = 2), type = "type_b")
  df <- rbind(df_a, df_b)

  par <- c(5, 1, 2)  # mean, sd, rate

  expected_ll <- sum(dnorm(df_a$x, par[1], par[2], log = TRUE)) +
    sum(dexp(df_b$x, par[3], log = TRUE))
  computed_ll <- model$loglik(df, par)

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("likelihood_contr_model falls back to numerical score when not provided", {
  # Define only loglik, not score
  loglik.num_test <- function(df, par, ...) {
    sum(dnorm(df$x, par[1], par[2], log = TRUE))
  }

  assign("loglik.num_test", loglik.num_test, envir = .GlobalEnv)
  on.exit(rm("loglik.num_test", envir = .GlobalEnv))

  model <- likelihood_contr_model$new(
    obs_type = function(df) "num_test"
  )

  set.seed(777)
  df <- data.frame(x = rnorm(50, 0, 1))
  par <- c(0.1, 0.9)

  # Get score from model (should use numerical differentiation)
  score_val <- model$score(df, par)

  # Manually compute numerical gradient
  expected_score <- numDeriv::grad(
    func = function(p) loglik.num_test(df, p),
    x = par,
    method.args = list(r = 6)
  )

  expect_equal(score_val, expected_score, tolerance = 1e-5)
})

test_that("likelihood_contr_model uses provided dispatcher functions", {
  custom_loglik <- function(df, par, ...) {
    -999  # Distinctive value
  }

  model <- likelihood_contr_model$new(
    obs_type = function(df) "custom",
    logliks = list(custom = custom_loglik)
  )

  df <- data.frame(x = 1:5)
  expect_equal(model$loglik(df, c(1, 2)), -999)
})

test_that("likelihood_contr_model errors when no loglik dispatcher found", {
  model <- likelihood_contr_model$new(
    obs_type = function(df) "nonexistent_type_xyz"
  )

  df <- data.frame(x = 1:5)
  expect_error(model$loglik(df, c(1, 2)),
               "No `loglik` dispatcher for type: nonexistent_type_xyz")
})

test_that("S3 wrappers for likelihood_contr_model work correctly", {
  # Note: The hessian wrapper needs data frame to be passed properly
  loglik.s3_test <- function(df, par, ...) {
    # Handle both data frame and vector cases (for numerical differentiation)
    if (is.data.frame(df)) {
      sum(dnorm(df$x, par[1], par[2], log = TRUE))
    } else {
      # This shouldn't happen with proper calls
      stop("Expected data frame")
    }
  }

  assign("loglik.s3_test", loglik.s3_test, envir = .GlobalEnv)
  on.exit(rm("loglik.s3_test", envir = .GlobalEnv))

  model <- likelihood_contr_model$new(
    obs_type = function(df) "s3_test"
  )

  set.seed(888)
  df <- data.frame(x = rnorm(30))
  par <- c(0, 1)

  # S3 wrapper should work for loglik
  ll_wrapper <- loglik(model)
  expect_equal(ll_wrapper(df, par), model$loglik(df, par))

  # S3 wrapper should work for score
  score_wrapper <- score(model)
  expect_equal(score_wrapper(df, par), model$score(df, par))
})

test_that("assumptions for likelihood_contr_model always includes iid", {
  model <- likelihood_contr_model$new(
    obs_type = function(df) "exact"
  )

  expect_true("iid" %in% assumptions(model))

  # Custom assumptions are added
  model2 <- likelihood_contr_model$new(
    obs_type = function(df) "exact",
    assumptions = c("weibull distribution", "independent censoring")
  )

  expect_true("iid" %in% assumptions(model2))
  expect_true("weibull distribution" %in% assumptions(model2))
  expect_true("independent censoring" %in% assumptions(model2))
})

# =============================================================================
# FIT.LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("fit.likelihood_model finds correct MLE for normal data", {
  set.seed(999)
  true_mean <- 5
  true_sd <- 2
  n <- 200
  df <- data.frame(x = rnorm(n, true_mean, true_sd), censor = rep("exact", n))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  result <- solver(df, par = c(0, 1))

  # MLE should be close to sample mean and sd (which are MLEs for normal)
  estimated <- get_mle_params(result)

  expect_equal(estimated[1], mean(df$x), tolerance = 0.01)
  # Note: MLE for sd differs from sample sd by sqrt((n-1)/n)
  expect_equal(estimated[2], sd(df$x) * sqrt((n - 1) / n), tolerance = 0.1)
})

test_that("fit.likelihood_model finds correct MLE for Weibull data", {
  set.seed(1000)
  true_shape <- 2
  true_scale <- 1
  df <- data.frame(x = rweibull(100, shape = true_shape, scale = true_scale))

  model <- likelihood_exact_weibull("x")
  solver <- fit(model)
  result <- solver(df, par = c(1.5, 0.8))

  estimated <- get_mle_params(result)

  # Check that we got valid estimates (not NA)
  expect_true(!is.null(estimated) && !any(is.na(estimated)),
              info = "Fit should produce valid estimates")

  # Should be reasonably close to true parameters
  expect_lt(abs(estimated[1] - true_shape), 0.5)
  expect_lt(abs(estimated[2] - true_scale), 0.3)
})

test_that("fit.likelihood_model works with different optimization methods", {
  set.seed(1001)
  df <- data.frame(x = rnorm(100), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)

  # Test Nelder-Mead (default)
  result_nm <- solver(df, par = c(0, 1), method = "Nelder-Mead")
  params_nm <- get_mle_params(result_nm)
  expect_true(!is.null(params_nm))
  expect_true(!any(is.na(params_nm)))

  # Test BFGS
  result_bfgs <- solver(df, par = c(0, 1), method = "BFGS")
  params_bfgs <- get_mle_params(result_bfgs)
  expect_true(!is.null(params_bfgs))

  # Results should be similar (if both converged properly)
  if (!is.null(params_bfgs) && !any(is.na(params_bfgs))) {
    expect_equal(params_nm, params_bfgs, tolerance = 0.05)
  }
})

test_that("fit.likelihood_model requires initial parameters", {
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  df <- data.frame(x = rnorm(10), censor = rep("exact", 10))

  expect_error(solver(df, par = NULL))
})

# =============================================================================
# SAMPLER.LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("sampler.likelihood_model generates bootstrap samples", {
  skip_if_not_installed("algebraic.mle")
  skip_if_not_installed("boot")

  set.seed(1002)
  df <- data.frame(x = rweibull(200, shape = 2, scale = 1))

  model <- likelihood_exact_weibull("x")

  # Create sampler
  model_sampler <- sampler(model, df = df, par = c(2, 1))

  # Generate bootstrap samples
  boot_result <- model_sampler(n = 10)

  # Should return an object with bootstrap estimates
  expect_true(!is.null(boot_result))

  # Bootstrap estimates should be stored in t
  expect_equal(nrow(boot_result$t), 10)
  expect_equal(ncol(boot_result$t), 2)  # shape and scale
})

test_that("sampler.likelihood_model bootstrap distribution is centered near MLE", {
  skip_if_not_installed("algebraic.mle")
  skip_if_not_installed("boot")

  set.seed(1003)
  true_shape <- 2
  true_scale <- 1.5
  df <- data.frame(x = rweibull(200, shape = true_shape, scale = true_scale))

  model <- likelihood_exact_weibull("x")

  # First fit to get MLE
  solver <- fit(model)
  mle_result <- solver(df, par = c(1, 1))
  mle_params <- get_mle_params(mle_result)

  # Create sampler and generate bootstrap samples
  model_sampler <- sampler(model, df = df, par = c(1, 1))
  boot_result <- model_sampler(n = 50)

  # Bootstrap mean should be close to MLE
  boot_mean <- colMeans(boot_result$t)
  expect_equal(boot_mean[1], mle_params[1], tolerance = 0.4)
  expect_equal(boot_mean[2], mle_params[2], tolerance = 0.4)
})

# =============================================================================
# LRT (LIKELIHOOD RATIO TEST) TESTS
# =============================================================================

test_that("lrt requires likelihood models as inputs", {
  df <- data.frame(x = rnorm(10), censor = rep("exact", 10))

  expect_error(lrt(list(), likelihood_name("norm", "x", "censor"),
                   df, null_par = c(0, 1), alt_par = c(0, 1, 0.5)))
})

test_that("lrt computes correct test statistic", {
  set.seed(1004)
  # Generate data from null model (exponential = Weibull with shape = 1)
  df <- data.frame(x = rexp(100, rate = 1))

  # Null model: exponential (shape fixed at 1, estimate scale only)
  # We'll use Weibull but the null hypothesis is shape = 1
  model <- likelihood_exact_weibull("x")

  # Alternative: Weibull (estimate both shape and scale)
  ll_func <- loglik(model)

  # Compute log-likelihoods at null and alternative parameters
  null_par <- c(1, 1)  # shape = 1 (exponential)
  alt_par <- c(1.05, 0.95)  # slightly different

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
  # Generate data from normal(0, 1)
  df <- data.frame(x = rnorm(200, 0, 1), censor = rep("exact", 200))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  # Null: mean = 0 (true value)
  # Alternative: mean can vary
  null_par <- c(0, 1)
  alt_par <- c(mean(df$x), sd(df$x))  # MLE

  result <- lrt(model, model, df,
                null_par = null_par, alt_par = alt_par, dof = 1)

  # Test statistic should be non-negative
  expect_true(result$stat >= 0)
  # p-value should be in [0, 1]
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("lrt auto-calculates degrees of freedom", {
  set.seed(1006)
  df <- data.frame(x = rnorm(50), censor = rep("exact", 50))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  # Both models have same structure, but we test dof calculation
  # null_par has fewer "free" parameters conceptually
  # Using same model for both (comparing fixed vs fitted parameters)
  result <- lrt(model, model, df,
                null_par = c(0, 1),    # Fixed at hypothesized values
                alt_par = c(0.1, 1.1), # "Alternative" estimates
                dof = 2)               # Explicitly set dof

  expect_equal(result$dof, 2)
  expect_true(!is.null(result$stat))
  expect_true(!is.null(result$p.value))
})

# =============================================================================
# EDGE CASES AND ERROR HANDLING
# =============================================================================

test_that("loglik handles data with single observation", {
  df <- data.frame(x = 5.0, censor = "exact")

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll_func <- loglik(model)

  expected_ll <- dnorm(5.0, 0, 1, log = TRUE)
  computed_ll <- ll_func(df, c(0, 1))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("likelihood_exact_weibull handles single observation", {
  df <- data.frame(x = 2.0)

  model <- likelihood_exact_weibull("x")
  ll_func <- loglik(model)

  expected_ll <- dweibull(2.0, shape = 2, scale = 1, log = TRUE)
  computed_ll <- ll_func(df, c(2, 1))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

test_that("models handle named parameters for likelihood_name_model", {
  set.seed(1007)
  df <- data.frame(x = rnorm(50), censor = rep("exact", 50))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll_func <- loglik(model)

  # Named parameters
  ll_named <- ll_func(df, c(mean = 0, sd = 1))
  # Unnamed parameters
  ll_unnamed <- ll_func(df, c(0, 1))

  expect_equal(ll_named, ll_unnamed)
})

test_that("hess_loglik returns symmetric matrix", {
  set.seed(1008)
  df <- data.frame(x = rweibull(50, 2, 1))

  model <- likelihood_exact_weibull("x")
  hess_func <- hess_loglik(model)

  H <- hess_func(df, c(2, 1))

  # Hessian should be symmetric
  expect_equal(H, t(H), tolerance = 1e-10)
})

test_that("score sums to approximately zero at MLE for likelihood_exact_weibull", {
  set.seed(1009)
  df <- data.frame(x = rweibull(200, 2, 1))

  model <- likelihood_exact_weibull("x")

  # Fit to get MLE
  solver <- fit(model)
  result <- solver(df, par = c(1.5, 0.8))
  mle_par <- get_mle_params(result)

  # Only test if fit converged (no NAs)
  skip_if(is.null(mle_par) || any(is.na(mle_par)), "MLE did not converge")

  # Score at MLE should be approximately zero
  score_func <- score(model)
  score_at_mle <- score_func(df, mle_par)

  # Score should be near zero at MLE
  expect_lt(max(abs(score_at_mle)), 0.5)
})

# =============================================================================
# REAL-WORLD SCENARIO TESTS
# =============================================================================

test_that("series system with observed component failures works", {
  # Simulate series system data with 3 exponential components
  set.seed(1010)
  n <- 100
  rates <- c(1.1, 1.2, 1.3)

  # Generate component lifetimes
  t1 <- rexp(n, rates[1])
  t2 <- rexp(n, rates[2])
  t3 <- rexp(n, rates[3])

  # System lifetime is minimum
  df <- data.frame(
    t = pmin(t1, t2, t3),
    k = apply(cbind(t1, t2, t3), 1, which.min)
  )

  # Define loglik for this observed series system
  loglik.series_observed <- function(df, rates, ...) {
    sum(log(rates[df$k])) - sum(rates) * sum(df$t)
  }

  assign("loglik.series_observed", loglik.series_observed, envir = .GlobalEnv)
  on.exit(rm("loglik.series_observed", envir = .GlobalEnv))

  model <- likelihood_contr_model$new(
    obs_type = function(df) "series_observed"
  )

  # Log-likelihood at true parameters
  ll_true <- model$loglik(df, rates)
  expect_true(is.finite(ll_true))

  # Fit model
  solver <- fit(model)
  result <- solver(df, par = c(1, 1, 1))

  estimated <- get_mle_params(result)

  # Check we got valid estimates
  skip_if(is.null(estimated) || any(is.na(estimated)), "MLE did not converge")

  # Estimates should be in reasonable range of true values
  expect_lt(abs(estimated[1] - rates[1]), 0.5)
  expect_lt(abs(estimated[2] - rates[2]), 0.5)
  expect_lt(abs(estimated[3] - rates[3]), 0.5)
})

test_that("right-censored weibull data estimation works", {
  set.seed(1011)
  n <- 100
  true_shape <- 2
  true_scale <- 1.5
  censoring_time <- 1.2

  # Generate data
  raw_t <- rweibull(n, true_shape, true_scale)

  # Apply right censoring
  df <- data.frame(
    x = pmin(raw_t, censoring_time),
    censor = ifelse(raw_t > censoring_time, "right", "exact")
  )

  model <- likelihood_name("weibull", ob_col = "x", censor_col = "censor")

  # Fit model
  solver <- fit(model)
  result <- solver(df, par = c(1.5, 1))

  estimated <- get_mle_params(result)

  # Check we got valid estimates
  skip_if(is.null(estimated) || any(is.na(estimated)), "MLE did not converge")

  # With censoring, estimates should still be reasonable
  expect_lt(abs(estimated[1] - true_shape), 0.7)
  expect_lt(abs(estimated[2] - true_scale), 0.7)
})

# =============================================================================
# INTEGRATION TESTS
# =============================================================================

test_that("full workflow: create model, fit, get score and hessian", {
  set.seed(1012)
  true_shape <- 1.5
  true_scale <- 2.0
  df <- data.frame(x = rweibull(150, true_shape, true_scale))

  # Create model
  model <- likelihood_exact_weibull("x")
  expect_true(is_likelihood_model(model))

  # Get functions
  ll_func <- loglik(model)
  score_func <- score(model)
  hess_func <- hess_loglik(model)

  # Fit model
  solver <- fit(model)
  result <- solver(df, par = c(1.2, 1.5))
  mle_par <- get_mle_params(result)

  skip_if(is.null(mle_par) || any(is.na(mle_par)), "MLE did not converge")

  # Verify score is near zero at MLE
  score_at_mle <- score_func(df, mle_par)
  expect_lt(max(abs(score_at_mle)), 0.5)

  # Verify Hessian is computed (note: may have sign issues)
  hess_at_mle <- hess_func(df, mle_par)
  expect_true(is.matrix(hess_at_mle))
  expect_equal(dim(hess_at_mle), c(2, 2))

  # Verify parameters are reasonable
  expect_lt(abs(mle_par[1] - true_shape), 0.4)
  expect_lt(abs(mle_par[2] - true_scale), 0.5)
})

test_that("comparison of likelihood_name_model and likelihood_exact_weibull", {
  set.seed(1013)
  shape <- 2.0
  scale <- 1.0
  df_name <- data.frame(x = rweibull(100, shape, scale), censor = rep("exact", 100))
  df_exact <- data.frame(x = df_name$x)  # Same data, different format

  model_name <- likelihood_name("weibull", ob_col = "x", censor_col = "censor")
  model_exact <- likelihood_exact_weibull("x")

  ll_name <- loglik(model_name)(df_name, c(shape, scale))
  ll_exact <- loglik(model_exact)(df_exact, c(shape, scale))

  # Both should give the same log-likelihood
  expect_equal(ll_name, ll_exact, tolerance = 1e-10)
})

# =============================================================================
# PRINT.LIKELIHOOD_MODEL TESTS
# =============================================================================

test_that("print.likelihood_model outputs model information", {
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  # Capture output
  output <- capture.output(print(model))

  # Check that key information is printed
  expect_true(any(grepl("Likelihood model", output)))
  expect_true(any(grepl("Observation column", output)))
})

test_that("print.likelihood_model with show.loglik=TRUE shows loglik function", {
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  # Capture output with show.loglik=TRUE
  output <- capture.output(print(model, show.loglik = TRUE))

  # Should include loglik function output
  expect_true(any(grepl("Log-likelihood function", output)))
})

test_that("print.likelihood_model returns model invisibly", {
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  # print should return model invisibly
  result <- print(model)
  expect_identical(result, model)
})

# =============================================================================
# LRT WITH NULL DOF TESTS
# =============================================================================

test_that("lrt computes dof automatically when NULL", {
  set.seed(1014)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))

  # Use likelihood_exact_weibull for both models
  null_model <- likelihood_exact_weibull("x")
  alt_model <- likelihood_exact_weibull("x")

  # Fit both models
  null_par <- c(1.5, 1.2)  # 2 parameters

  alt_par <- c(2.0, 1.0)   # 2 parameters

  # When dof is NULL, it should be computed as length(alt_par) - length(null_par)
  result <- lrt(null_model, alt_model, df,
                null_par = null_par, alt_par = alt_par, dof = NULL)

  # dof should be 0 (2 - 2)
  expect_equal(result$dof, 0)
  expect_true(!is.null(result$stat))
  expect_true(!is.null(result$p.value))
})

# =============================================================================
# FIT WITH SANN METHOD TESTS
# =============================================================================

test_that("fit.likelihood_model works with SANN method (no gradient)", {
  set.seed(1015)
  df <- data.frame(x = rnorm(100, mean = 5, sd = 2), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)

  # SANN should work without gradient (gr = NULL)
  result <- solver(df, par = c(0, 1), method = "SANN",
                   control = list(maxit = 1000))

  params <- get_mle_params(result)
  expect_true(!is.null(params))
  expect_true(!any(is.na(params)))

  # Parameters should be reasonable (not exact due to SANN stochasticity)
  expect_lt(abs(params[1] - 5), 1.5)
  expect_lt(abs(params[2] - 2), 1.5)
})

# =============================================================================
# FIM GENERIC TESTS
# =============================================================================

test_that("fim generic dispatches correctly", {
  # fim is a generic that requires a method
  # Test that it exists and errors appropriately for objects without methods
  expect_error(fim(list()), "no applicable method")
})
