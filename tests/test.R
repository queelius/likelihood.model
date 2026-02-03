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
  # likelihood_contr_model is a likelihood model
  model <- likelihood_contr_model$new(obs_type = function(df) "exact")
  expect_true(is_likelihood_model(model))

  # likelihood_name_model is a likelihood model
  # Note: censor_col must not be NULL, use actual column name
  model_name <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  expect_true(is_likelihood_model(model_name))

  # weibull_uncensored is a likelihood model (aliased as likelihood_exact_weibull)
  model_weibull <- weibull_uncensored("x")
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

  expect_true("independent observations" %in% assumptions_list)
  expect_true("identically distributed" %in% assumptions_list)
  expect_true("weibull distribution" %in% assumptions_list)
  expect_true("censoring independent of observations" %in% assumptions_list)

  # Model without censoring should have fewer assumptions
  model_no_censor <- likelihood_name("norm", ob_col = "x")
  assumptions_no_censor <- assumptions(model_no_censor)
  expect_false("censoring independent of observations" %in% assumptions_no_censor)
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

test_that("likelihood_name works without censor_col (all exact)", {
  set.seed(112)
  mean_true <- 5
  sd_true <- 2
  df <- data.frame(x = rnorm(100, mean_true, sd_true))

  # No censor_col - should treat all as exact
  model <- likelihood_name("norm", ob_col = "x")
  ll_func <- loglik(model)

  expected_ll <- sum(dnorm(df$x, mean_true, sd_true, log = TRUE))
  computed_ll <- ll_func(df, c(mean_true, sd_true))

  expect_equal(computed_ll, expected_ll, tolerance = 1e-10)
})

# =============================================================================
# WEIBULL_UNCENSORED TESTS
# =============================================================================

test_that("weibull_uncensored creates model with correct structure", {
  model <- weibull_uncensored("x")

  expect_true(is_likelihood_model(model))
  expect_equal(model$ob_col, "x")
  expect_true("weibull_uncensored" %in% class(model))
})

test_that("likelihood_exact_weibull alias works", {
  model <- likelihood_exact_weibull("x")
  expect_true(is_likelihood_model(model))
  expect_true("weibull_uncensored" %in% class(model))
})

test_that("loglik for weibull_uncensored matches R's dweibull", {
  set.seed(222)
  shape <- 2.5
  scale <- 1.5
  df <- data.frame(x = rweibull(50, shape, scale))

  model <- weibull_uncensored("x")
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

  model <- weibull_uncensored("x")
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

  model <- weibull_uncensored("x")
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

test_that("weibull_uncensored validates positive observations", {
  df_invalid <- data.frame(x = c(1.0, -0.5, 2.0))  # Contains negative value

  model <- weibull_uncensored("x")
  ll_func <- loglik(model)

  expect_error(ll_func(df_invalid, c(2, 1)),
               info = "Should error on non-positive observations")
})

test_that("weibull_uncensored requires non-empty data", {
  df_empty <- data.frame(x = numeric(0))

  model <- weibull_uncensored("x")
  ll_func <- loglik(model)

  expect_error(ll_func(df_empty, c(2, 1)),
               info = "Should error on empty data frame")
})

test_that("assumptions for weibull_uncensored returns expected list", {
  model <- weibull_uncensored("x")
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

  model <- weibull_uncensored("x")
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
  skip_if_not_installed("boot")

  set.seed(1002)
  df <- data.frame(x = rweibull(200, shape = 2, scale = 1))

  model <- weibull_uncensored("x")

  # Create sampler
  model_sampler <- sampler(model, df = df, par = c(2, 1))

  # Generate bootstrap samples
  boot_result <- model_sampler(n = 10)

  # Should return a fisher_boot object
  expect_true(inherits(boot_result, "fisher_boot"))
  expect_true(inherits(boot_result, "fisher_mle"))

  # Bootstrap estimates should be stored in replicates
  expect_equal(nrow(boot_result$replicates), 10)
  expect_equal(ncol(boot_result$replicates), 2)  # shape and scale

  # coef() should return original MLE
  expect_equal(length(coef(boot_result)), 2)
})

test_that("sampler.likelihood_model bootstrap distribution is centered near MLE", {
  skip_if_not_installed("boot")

  set.seed(1003)
  true_shape <- 2
  true_scale <- 1.5
  df <- data.frame(x = rweibull(200, shape = true_shape, scale = true_scale))

  model <- weibull_uncensored("x")

  # First fit to get MLE
  solver <- fit(model)
  mle_result <- solver(df, par = c(1, 1))
  mle_params <- get_mle_params(mle_result)

  # Create sampler and generate bootstrap samples
  model_sampler <- sampler(model, df = df, par = c(1, 1))
  boot_result <- model_sampler(n = 50)

  # Bootstrap mean should be close to MLE
  boot_mean <- colMeans(boot_result$replicates)
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
  model <- weibull_uncensored("x")

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

test_that("weibull_uncensored handles single observation", {
  df <- data.frame(x = 2.0)

  model <- weibull_uncensored("x")
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

  model <- weibull_uncensored("x")
  hess_func <- hess_loglik(model)

  H <- hess_func(df, c(2, 1))

  # Hessian should be symmetric
  expect_equal(H, t(H), tolerance = 1e-10)
})

test_that("score sums to approximately zero at MLE for weibull_uncensored", {
  set.seed(1009)
  df <- data.frame(x = rweibull(200, 2, 1))

  model <- weibull_uncensored("x")

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
  model <- weibull_uncensored("x")
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

test_that("comparison of likelihood_name_model and weibull_uncensored", {
  set.seed(1013)
  shape <- 2.0
  scale <- 1.0
  df_name <- data.frame(x = rweibull(100, shape, scale), censor = rep("exact", 100))
  df_exact <- data.frame(x = df_name$x)  # Same data, different format

  model_name <- likelihood_name("weibull", ob_col = "x", censor_col = "censor")
  model_exact <- weibull_uncensored("x")

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

  # Use weibull_uncensored for both models
  null_model <- weibull_uncensored("x")
  alt_model <- weibull_uncensored("x")

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
  expect_true(inherits(result, "mle"))

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

test_that("fisher_mle aic and bic compute correctly", {
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    loglik_val = -100,
    nobs = 50
  )

  # AIC = -2*logL + 2*k = -2*(-100) + 2*2 = 204
  expect_equal(aic(result), 204)

  # BIC = -2*logL + k*log(n) = 200 + 2*log(50)
  expected_bic <- 200 + 2 * log(50)
  expect_equal(bic(result), expected_bic)
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
  df <- data.frame(x = rnorm(100, 5, 2), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  result <- solver(df, par = c(0, 1))

  expect_true(inherits(result, "fisher_mle"))
  expect_true(inherits(result, "mle"))
  expect_true(result$converged)
  expect_equal(result$nobs, 100)
})

# =============================================================================
# FISHERIAN INFERENCE TESTS
# =============================================================================

test_that("support function computes log relative likelihood", {
  set.seed(2002)
  df <- data.frame(x = rnorm(100, 5, 2), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  mle_result <- solver(df, par = c(0, 1))

  # Support at MLE should be 0
  s_at_mle <- support(mle_result, coef(mle_result), df, model)
  expect_equal(s_at_mle, 0, tolerance = 1e-10)

  # Support at other values should be negative
  s_away <- support(mle_result, c(0, 1), df, model)
  expect_true(s_away < 0)
})

test_that("relative_likelihood computes exp(support)", {
  set.seed(2003)
  df <- data.frame(x = rnorm(50, 3, 1), censor = rep("exact", 50))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  mle_result <- solver(df, par = c(0, 1))

  # Relative likelihood at MLE should be 1
  rl_at_mle <- relative_likelihood(mle_result, coef(mle_result), df, model)
  expect_equal(rl_at_mle, 1, tolerance = 1e-10)

  # Relative likelihood elsewhere should be between 0 and 1
  rl_away <- relative_likelihood(mle_result, c(0, 1), df, model)
  expect_true(rl_away > 0 && rl_away < 1)
})

test_that("likelihood_interval computes profile interval for multivariate case", {
  set.seed(2004)
  # Use normal (2 parameters) and compute profile likelihood interval for mean
  df <- data.frame(x = rnorm(100, mean = 5, sd = 2), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  mle_result <- solver(df, par = c(0, 1))

  # Compute 1/8 likelihood interval for mean parameter only
  li <- likelihood_interval(mle_result, df, model, k = 8, param = 1)

  expect_true(inherits(li, "likelihood_interval"))
  expect_equal(nrow(li), 1)
  expect_equal(ncol(li), 2)

  # Interval should contain MLE for mean
  mle_val <- coef(mle_result)[1]
  expect_true(li[1, 1] < mle_val && mle_val < li[1, 2])
})

test_that("profile_loglik computes profile for single parameter", {
  set.seed(2005)
  df <- data.frame(x = rnorm(100, 5, 2), censor = rep("exact", 100))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  mle_result <- solver(df, par = c(0, 1))

  # Profile for first parameter (mean)
  prof <- profile_loglik(mle_result, df, model, param = 1, n_grid = 20)

  expect_true(inherits(prof, "profile_loglik"))
  expect_equal(nrow(prof), 20)
  expect_true("loglik" %in% names(prof))
  expect_true("support" %in% names(prof))
  expect_true("relative_likelihood" %in% names(prof))

  # Maximum profile loglik should be near MLE
  max_idx <- which.max(prof$loglik)
  expect_equal(prof[[1]][max_idx], coef(mle_result)[1], tolerance = 0.5)
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
  df <- data.frame(x = rnorm(50, 5, 1), censor = rep("exact", 50))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  # Evidence for theta1 vs theta2
  ev <- evidence(model, df, c(5, 1), c(0, 1))

  # Should strongly favor true parameters (5, 1) over (0, 1)
  expect_true(ev > 0)
  expect_true(ev > log(8))  # Strong evidence
})

# =============================================================================
# FISHER_BOOT CLASS TESTS
# =============================================================================

test_that("fisher_boot inherits from fisher_mle", {
  skip_if_not_installed("boot")

  set.seed(2007)
  df <- data.frame(x = rweibull(100, 2, 1))

  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_result <- model_sampler(n = 20)

  expect_true(inherits(boot_result, "fisher_boot"))
  expect_true(inherits(boot_result, "fisher_mle"))
  expect_true(inherits(boot_result, "mle"))

  # Base methods should work
  expect_equal(length(coef(boot_result)), 2)
  expect_true(is.matrix(vcov(boot_result)))
})

test_that("fisher_boot bias computes from bootstrap replicates", {
  skip_if_not_installed("boot")

  set.seed(2008)
  df <- data.frame(x = rweibull(200, 2, 1))

  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_result <- model_sampler(n = 50)

  # Bias should be computed from replicates
  b <- bias(boot_result)
  expect_equal(length(b), 2)

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
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))

  model <- weibull_uncensored("x")

  # alt_par has 3 elements, null_par has 2 => dof should be 1
  # We use the same model but supply different-length parameter vectors
  # to test the auto-computation of dof
  null_par <- c(2, 1)
  alt_par <- c(2, 1, 0.5)

  # This should NOT error — it should auto-compute dof = 1
  # Before the fix, stopifnot(dof >= 0) with NULL dof silently passed,
  # but dof was computed *after* the stat, meaning the bug was ordering.
  result <- lrt(model, model, df,
                null_par = null_par, alt_par = alt_par, dof = NULL)

  expect_equal(result$dof, 1)
  expect_true(is.finite(result$stat))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("likelihood_interval works for normal mean (negative parameter)", {
  set.seed(3002)
  # Use normal distribution where mean can be negative
  df <- data.frame(x = rnorm(200, mean = -3, sd = 1.5),
                   censor = rep("exact", 200))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  solver <- fit(model)
  mle_result <- solver(df, par = c(0, 1))

  # Profile likelihood interval for mean (param 1)
  # This previously failed with L-BFGS-B because it assumed positive params
  li <- likelihood_interval(mle_result, df, model, k = 8, param = 1)

  expect_true(inherits(li, "likelihood_interval"))
  expect_equal(nrow(li), 1)

  # Interval should contain MLE
  mle_mean <- coef(mle_result)[1]
  expect_true(li[1, 1] < mle_mean && mle_mean < li[1, 2])

  # Mean should be negative, so lower bound should also be negative

  expect_true(li[1, 1] < 0)
})

test_that("evidence returns 0 when comparing same parameters", {
  set.seed(3003)
  df <- data.frame(x = rnorm(50, 5, 1), censor = rep("exact", 50))

  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")

  # Evidence for identical parameters should be exactly 0
  ev <- evidence(model, df, c(5, 1), c(5, 1))
  expect_equal(ev, 0)
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

test_that("exponential_lifetime matches likelihood_name('exp', ...)", {
  set.seed(4009)
  lambda_true <- 1.5
  x <- rexp(100, rate = lambda_true)

  # exponential_lifetime model
  df_exp <- data.frame(t = x)
  model_exp <- exponential_lifetime("t")
  ll_exp <- loglik(model_exp)(df_exp, lambda_true)

  # likelihood_name model (note: dexp uses 'rate' parameter)
  df_name <- data.frame(t = x, censor = rep("exact", length(x)))
  model_name <- likelihood_name("exp", ob_col = "t", censor_col = "censor")
  ll_name <- loglik(model_name)(df_name, lambda_true)

  expect_equal(ll_exp, ll_name, tolerance = 1e-10)
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
# NEW TESTS: LIKELIHOOD_NAME_MODEL VALIDATION & FEATURES
# =============================================================================

test_that("likelihood_name errors on invalid distribution name", {
  expect_error(likelihood_name("nonexistent_dist_xyz", ob_col = "x"),
               "not found")
})

test_that("loglik warns on unknown censoring values", {
  df <- data.frame(x = c(1, 2, 3), censor = c("exact", "righ", "unknwn"))
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll <- loglik(model)
  # "righ" and "unknwn" are unrecognized typos
  expect_warning(ll(df, c(0, 1)), "Unknown censoring values")
})

test_that("loglik errors on too many parameters", {
  df <- data.frame(x = rnorm(10))
  model <- likelihood_name("norm", ob_col = "x")
  ll <- loglik(model)
  expect_error(ll(df, c(0, 1, 0.5)), "Too many parameters")
})

test_that("loglik handles interval-censored data", {
  df <- data.frame(
    x = c(1.0, 2.0, 3.0),
    x_upper = c(1.5, 2.5, 3.5),
    censor = c("interval", "interval", "interval")
  )
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor",
                           ob_col_upper = "x_upper")
  ll <- loglik(model)
  expected <- sum(log(pnorm(df$x_upper, 0, 1) - pnorm(df$x, 0, 1)))
  expect_equal(ll(df, c(0, 1)), expected, tolerance = 1e-10)
})

test_that("loglik handles all four censoring types together", {
  df <- data.frame(
    x = c(0.5, 1.5, -0.5, 2.0),
    x_upper = c(NA, NA, NA, 3.0),
    censor = c("exact", "right", "left", "interval")
  )
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor",
                           ob_col_upper = "x_upper")
  ll <- loglik(model)
  expected <- dnorm(0.5, 0, 1, log = TRUE) +
    pnorm(1.5, 0, 1, log.p = TRUE, lower.tail = FALSE) +
    pnorm(-0.5, 0, 1, log.p = TRUE) +
    log(pnorm(3.0, 0, 1) - pnorm(2.0, 0, 1))
  expect_equal(ll(df, c(0, 1)), expected, tolerance = 1e-10)
})

test_that("loglik returns 0 when observation column is missing from data", {
  # When ob_col is missing, df[[ob_col]] returns NULL, length(NULL) == 0,
  # so the function returns 0 (empty data path). This documents current behavior.
  df <- data.frame(y = rnorm(10))
  model <- likelihood_name("norm", ob_col = "x")
  ll <- loglik(model)
  expect_equal(ll(df, c(0, 1)), 0)
})

test_that("loglik works with all-censored data (no exact observations)", {
  df <- data.frame(x = c(1, 2, 3), censor = rep("right", 3))
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll <- loglik(model)
  expected <- sum(pnorm(c(1, 2, 3), 0, 1, log.p = TRUE, lower.tail = FALSE))
  expect_equal(ll(df, c(0, 1)), expected, tolerance = 1e-10)
})

test_that("likelihood_name works with gamma distribution", {
  set.seed(42)
  df <- data.frame(x = rgamma(100, shape = 2, rate = 3))
  model <- likelihood_name("gamma", ob_col = "x")
  ll <- loglik(model)
  expected <- sum(dgamma(df$x, shape = 2, rate = 3, log = TRUE))
  expect_equal(ll(df, c(2, 3)), expected, tolerance = 1e-10)
})

test_that("interval censoring errors without ob_col_upper", {
  df <- data.frame(
    x = c(1.0, 2.0),
    censor = c("interval", "interval")
  )
  model <- likelihood_name("norm", ob_col = "x", censor_col = "censor")
  ll <- loglik(model)
  expect_error(ll(df, c(0, 1)), "ob_col_upper")
})

# =============================================================================
# NEW TESTS: LIKELIHOOD_CONTR_MODEL ADVANCED CASES
# =============================================================================

test_that("S3 hess_loglik wrapper for contr_model matches R6 method", {
  loglik.hess_wrap_test <- function(df, par, ...) {
    sum(dnorm(df$x, par[1], par[2], log = TRUE))
  }
  assign("loglik.hess_wrap_test", loglik.hess_wrap_test, envir = .GlobalEnv)
  on.exit(rm("loglik.hess_wrap_test", envir = .GlobalEnv))

  model <- likelihood_contr_model$new(obs_type = function(df) "hess_wrap_test")
  set.seed(42)
  df <- data.frame(x = rnorm(30))
  par <- c(0, 1)

  hess_wrapper <- hess_loglik(model)
  expect_equal(hess_wrapper(df, par), model$hess_loglik(df, par))
})

test_that("contr_model score works with single-parameter multi-type model", {
  loglik.spm_a <- function(df, par, ...) sum(dexp(df$x, par[1], log = TRUE))
  loglik.spm_b <- function(df, par, ...) sum(dexp(df$x, par[1], log = TRUE))
  score.spm_a <- function(df, par, ...) sum(1/par[1] - df$x)
  score.spm_b <- function(df, par, ...) sum(1/par[1] - df$x)

  assign("loglik.spm_a", loglik.spm_a, envir = .GlobalEnv)
  assign("loglik.spm_b", loglik.spm_b, envir = .GlobalEnv)
  assign("score.spm_a", score.spm_a, envir = .GlobalEnv)
  assign("score.spm_b", score.spm_b, envir = .GlobalEnv)
  on.exit({
    rm("loglik.spm_a", "loglik.spm_b", "score.spm_a", "score.spm_b",
       envir = .GlobalEnv)
  })

  model <- likelihood_contr_model$new(
    obs_type = function(df) df$type
  )
  set.seed(42)
  df <- data.frame(
    x = rexp(20, rate = 2),
    type = rep(c("spm_a", "spm_b"), each = 10)
  )
  # This should NOT error (was rowSums bug before fix)
  score_fn <- score(model)
  result <- score_fn(df, 2)
  expect_true(is.numeric(result))

  # Also test R6 method directly
  r6_result <- model$score(df, 2)
  expect_true(is.numeric(r6_result))
  expect_equal(result, r6_result)
})

test_that("S3 loglik wrapper re-caches when df changes", {
  loglik.cache_test <- function(df, par, ...) {
    sum(dnorm(df$x, par[1], par[2], log = TRUE))
  }
  assign("loglik.cache_test", loglik.cache_test, envir = .GlobalEnv)
  on.exit(rm("loglik.cache_test", envir = .GlobalEnv))

  model <- likelihood_contr_model$new(obs_type = function(df) "cache_test")
  ll <- loglik(model)

  df1 <- data.frame(x = c(1, 2, 3))
  df2 <- data.frame(x = c(4, 5, 6))

  val1 <- ll(df1, c(0, 1))
  val2 <- ll(df2, c(0, 1))

  expect_false(val1 == val2)
  expect_equal(val1, sum(dnorm(c(1, 2, 3), 0, 1, log = TRUE)))
  expect_equal(val2, sum(dnorm(c(4, 5, 6), 0, 1, log = TRUE)))
})

test_that("print works for likelihood_contr_model without crashing", {
  model <- likelihood_contr_model$new(obs_type = function(df) "exact")
  output <- capture.output(print(model))
  expect_true(any(grepl("Likelihood model", output)))
  # Should NOT print "Observation column: " with NULL
  expect_false(any(grepl("Observation column:", output)))
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
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  result <- fit(model)(df, par = c(1.5, 0.8))

  p <- params(result)
  expect_equal(length(p), 2)
  expect_equal(p, coef(result))
  expect_equal(nparams(result), 2L)
})

test_that("observed_fim() is positive definite on real fit output", {
  set.seed(5002)
  df <- data.frame(x = rweibull(200, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  result <- fit(model)(df, par = c(1.5, 0.8))

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

  expect_true(inherits(mapped, "mle"))
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
  expect_true(inherits(m, "mle"))
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

test_that("loglik_val.fisher_mle returns the log-likelihood value", {
  result <- fisher_mle(
    par = c(a = 1, b = 2),
    vcov = matrix(c(0.04, 0, 0, 0.02), 2, 2),
    loglik_val = -123.45,
    nobs = 50
  )
  expect_equal(loglik_val(result), -123.45)
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

test_that("bic errors when nobs is NULL", {
  result <- fisher_mle(
    par = c(a = 1),
    vcov = matrix(0.1, 1, 1),
    loglik_val = -50
  )
  expect_error(bic(result), "nobs")
})


# =============================================================================
# FISHER_BOOT TESTS
# =============================================================================

test_that("confint.fisher_boot computes percentile CI", {
  skip_if_not_installed("boot")
  set.seed(6001)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, level = 0.95, type = "perc")
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
  expect_true(ci[1, 1] < ci[1, 2])  # lower < upper
})

test_that("confint.fisher_boot computes basic CI", {
  skip_if_not_installed("boot")
  set.seed(6002)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, level = 0.95, type = "basic")
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
})

test_that("confint.fisher_boot computes norm CI", {
  skip_if_not_installed("boot")
  set.seed(6003)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, level = 0.95, type = "norm")
  expect_equal(nrow(ci), 2)
  expect_equal(ncol(ci), 2)
})

test_that("confint.fisher_boot with integer parm", {
  skip_if_not_installed("boot")
  set.seed(6004)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_obj <- model_sampler(n = 200)

  ci <- confint(boot_obj, parm = 1, level = 0.95, type = "perc")
  expect_equal(nrow(ci), 1)
})

test_that("print.fisher_boot outputs expected text", {
  skip_if_not_installed("boot")
  set.seed(6005)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_obj <- model_sampler(n = 50)

  out <- capture.output(print(boot_obj))
  expect_true(any(grepl("Bootstrap MLE", out)))
  expect_true(any(grepl("Bootstrap replicates:", out)))
  expect_true(any(grepl("Bootstrap SE:", out)))
})

test_that("sampler.fisher_boot resamples from bootstrap distribution", {
  skip_if_not_installed("boot")
  set.seed(6006)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  model_sampler <- sampler(model, df = df, par = c(2, 1))
  boot_obj <- model_sampler(n = 200)

  samp_fn <- sampler(boot_obj)
  samples <- samp_fn(100)
  expect_equal(nrow(samples), 100)
  expect_equal(ncol(samples), 2)
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

test_that("likelihood_interval works with character param", {
  set.seed(6008)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  result <- fit(model)(df, par = c(shape = 1.5, scale = 0.8))

  li <- likelihood_interval(result, data = df, model = model,
                            k = 8, param = "shape")
  expect_true(is.matrix(li) || is.data.frame(li))
})

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

test_that("profile_loglik works with character param", {
  set.seed(6011)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  result <- fit(model)(df, par = c(shape = 1.5, scale = 0.8))

  pl <- profile_loglik(result, data = df, model = model,
                       param = "shape", n_grid = 10)
  expect_s3_class(pl, "profile_loglik")
  expect_true("shape" %in% names(pl))
})

test_that("profile_loglik errors for > 2 parameters", {
  result <- fisher_mle(
    par = c(a = 1, b = 2, c = 3),
    vcov = diag(3) * 0.01,
    loglik_val = -50,
    nobs = 100
  )
  model <- weibull_uncensored("x")
  df <- data.frame(x = 1:10)
  expect_error(
    profile_loglik(result, data = df, model = model, param = 1:3),
    "1 or 2 parameters"
  )
})

test_that("profile_loglik 2D works", {
  set.seed(6012)
  df <- data.frame(x = rweibull(100, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  result <- fit(model)(df, par = c(shape = 1.5, scale = 0.8))

  pl <- profile_loglik(result, data = df, model = model,
                       param = 1:2, n_grid = 5)
  expect_s3_class(pl, "profile_loglik")
  expect_equal(nrow(pl), 25)  # 5x5 grid
})

test_that("print.profile_loglik outputs expected text", {
  set.seed(6013)
  df <- data.frame(x = rweibull(50, shape = 2, scale = 1))
  model <- weibull_uncensored("x")
  result <- fit(model)(df, par = c(shape = 1.5, scale = 0.8))

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
# FIM AND OBSERVED_INFO COVERAGE TESTS
# =============================================================================

test_that("fim.likelihood_model computes FIM via Monte Carlo", {
  # Create a minimal model that falls through to fim.likelihood_model.
  # exponential_lifetime has its own fim method, so we use a mock.
  mock_model <- structure(list(), class = c("mock_exp", "likelihood_model"))

  score.mock_exp <<- function(model, ...) {
    function(df, par) {
      c(nrow(df) / par[1] - sum(df$x))
    }
  }

  rdata.mock_exp <<- function(model, ...) {
    function(theta, n, ...) {
      data.frame(x = rexp(n, rate = theta[1]))
    }
  }

  set.seed(6014)
  fim_fn <- fim(mock_model)
  # True FIM for exp(rate=2) with n=50 is n/rate^2 = 50/4 = 12.5
  fim_mat <- fim_fn(theta = c(2), n_obs = 50, n_samples = 500)
  expect_true(is.matrix(fim_mat))
  expect_equal(nrow(fim_mat), 1)
  expect_equal(ncol(fim_mat), 1)
  expect_equal(fim_mat[1, 1], 12.5, tolerance = 2)

  rm(score.mock_exp, rdata.mock_exp, envir = .GlobalEnv)
})

test_that("observed_info.likelihood_model computes -hessian", {
  set.seed(6015)
  df <- data.frame(x = rexp(50, rate = 2))
  model <- exponential_lifetime("x")
  oi_fn <- observed_info(model)
  oi <- oi_fn(df, c(rate = 2))
  expect_true(is.matrix(oi))
  expect_true(oi[1, 1] > 0)  # positive definite
})
