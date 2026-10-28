library(testthat)
library(likelihood.model)

# The nuisance optimum crosses zero although both fitted coordinates are positive.
loglik.correlated_domain <- function(model, ...) {
  function(df, par, ...) -0.5 * ((par[1] - 1)^2 + (par[2] - par[1])^2)
}
registerS3method("loglik", "correlated_domain", loglik.correlated_domain, envir = asNamespace("likelihood.model"))
model <- structure(list(), class = c("correlated_domain", "likelihood_model"))
fitted <- fisher_mle(c(a = 1, b = 1), vcov = matrix(c(1, 1, 1, 2), 2), loglik_val = 0)
test_that("positive nuisance estimates do not imply positive parameter domains", {
  interval <- likelihood_interval(fitted, data.frame(x = 0), model, param = "a")
  expected <- 1 + c(-1, 1) * sqrt(2 * log(8))
  expect_equal(as.numeric(interval), expected, tolerance = 1e-5)
  profile <- profile_loglik(fitted, data.frame(x = 0), model, param = "a", grid = expected)
  expect_equal(profile$support, rep(-log(8), 2), tolerance = 1e-5)
})

loglik.normal_domain <- function(model, ...) {
  function(df, par, ...) {
    if (par["sigma"] <= 0) return(-Inf)
    sum(dnorm(df$x, mean = par["mu"], sd = par["sigma"], log = TRUE))
  }
}
registerS3method("loglik", "normal_domain", loglik.normal_domain,
                 envir = asNamespace("likelihood.model"))

test_that("profile inference respects small parameter scales without sign constraints", {
  normal_model <- structure(list(), class = c("normal_domain", "likelihood_model"))
  for (scale in c(1, 1e-3, 1e-6)) {
    df <- data.frame(x = scale * c(-2, -1, 0, 1, 2))
    sigma2 <- mean(df$x^2)
    fitted <- fisher_mle(c(mu = 0, sigma = sqrt(sigma2)),
      vcov = diag(c(sigma2 / nrow(df), sigma2 / (2 * nrow(df)))),
      loglik_val = loglik(normal_model)(df, c(mu = 0, sigma = sqrt(sigma2))))
    # Analytic normal profile: -n/2 * log(1 + mu^2 / sigma_hat^2).
    expected <- c(-1, 1) * sqrt(sigma2 * (8^(2 / nrow(df)) - 1))
    interval <- likelihood_interval(fitted, df, normal_model, param = "mu")
    expect_equal(as.numeric(interval) / scale, expected / scale, tolerance = 1e-5)
    profile <- profile_loglik(fitted, df, normal_model, param = "mu", grid = expected)
    expect_equal(profile$support, rep(-log(8), 2), tolerance = 1e-5)
    fitted$vcov[2, 2] <- NA_real_
    profile_without_se <- profile_loglik(fitted, df, normal_model,
                                         param = "mu", grid = expected)
    expect_equal(profile_without_se$support, rep(-log(8), 2), tolerance = 1e-5)
  }
})
