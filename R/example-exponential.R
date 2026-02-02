#' @title Exponential lifetime model with right-censoring support
#' @description
#' A likelihood model for the Exponential(lambda) distribution with
#' optional right-censoring. This is a reference implementation
#' demonstrating:
#'
#' - **Closed-form MLE**: Overrides `fit()` to compute lambda_hat = d/T
#'   directly, bypassing `optim` entirely.
#' - **Analytical derivatives**: score, Hessian, and FIM in closed form.
#' - **Right-censoring**: Natural handling via the sufficient statistic
#'   (d, T) where d = number of exact observations and T = total time.
#' - **`rdata()` method**: For Monte Carlo validation and FIM estimation.
#'
#' The log-likelihood is:
#' \deqn{\ell(\lambda) = d \log \lambda - \lambda T}
#' where d is the number of exact (uncensored) observations and T is the
#' total observation time (sum of all times, whether censored or not).
#'
#' @param ob_col The name of the column containing observation times.
#' @param censor_col Optional column name indicating censoring status.
#'   When provided, values should be `"exact"` for uncensored observations
#'   and `"right"` for right-censored observations. When `NULL`, all
#'   observations are treated as exact.
#' @return An `exponential_lifetime` likelihood model object
#' @examples
#' # Uncensored exponential data
#' model <- exponential_lifetime("t")
#' df <- data.frame(t = rexp(100, rate = 2))
#' mle <- fit(model)(df)
#' coef(mle)  # should be close to 2
#'
#' # Right-censored data
#' model_c <- exponential_lifetime("t", censor_col = "status")
#' df_c <- data.frame(
#'   t = c(rexp(80, 2), rep(0.5, 20)),
#'   status = c(rep("exact", 80), rep("right", 20))
#' )
#' mle_c <- fit(model_c)(df_c)
#' coef(mle_c)
#' @export
exponential_lifetime <- function(ob_col, censor_col = NULL) {
  structure(
    list(ob_col = ob_col, censor_col = censor_col),
    class = c("exponential_lifetime", "likelihood_model")
  )
}


#' Assumptions for exponential_lifetime
#'
#' @param model An exponential_lifetime model
#' @param ... Additional arguments (ignored)
#' @return Character vector of assumptions
#' @export
assumptions.exponential_lifetime <- function(model, ...) {
  a <- c("independent", "identically distributed", "exponential distribution")
  if (!is.null(model$censor_col)) {
    a <- c(a, "non-informative right censoring")
  }
  a
}


#' Log-likelihood for exponential_lifetime
#'
#' Returns a function computing ell(lambda) = d * log(lambda) - lambda * T.
#'
#' @param model An exponential_lifetime model
#' @param ... Additional arguments (ignored)
#' @return A function(df, par, ...) computing the log-likelihood
#' @export
loglik.exponential_lifetime <- function(model, ...) {
  function(df, par, ...) {
    stats <- .exponential_sufficient_stats(df, model)
    lambda <- par[1]
    if (lambda <= 0) return(-.Machine$double.xmax / 2)
    stats$d * log(lambda) - lambda * stats$total_time
  }
}


#' Score for exponential_lifetime
#'
#' Returns a function computing d/lambda - T.
#'
#' @param model An exponential_lifetime model
#' @param ... Additional arguments (ignored)
#' @return A function(df, par, ...) computing the score vector
#' @export
score.exponential_lifetime <- function(model, ...) {
  function(df, par, ...) {
    stats <- .exponential_sufficient_stats(df, model)
    lambda <- par[1]
    stopifnot(lambda > 0)
    c(lambda = stats$d / lambda - stats$total_time)
  }
}


#' Hessian of the log-likelihood for exponential_lifetime
#'
#' Returns a function computing the 1x1 Hessian: -d/lambda^2.
#'
#' @param model An exponential_lifetime model
#' @param ... Additional arguments (ignored)
#' @return A function(df, par, ...) computing the Hessian matrix
#' @export
hess_loglik.exponential_lifetime <- function(model, ...) {
  function(df, par, ...) {
    stats <- .exponential_sufficient_stats(df, model)
    lambda <- par[1]
    stopifnot(lambda > 0)
    matrix(-stats$d / lambda^2, nrow = 1, ncol = 1,
           dimnames = list("lambda", "lambda"))
  }
}


#' Closed-form MLE for exponential_lifetime
#'
#' Computes the MLE directly as lambda_hat = d / T, bypassing `optim`.
#' This demonstrates that specialized models can provide exact solutions.
#'
#' @param object An exponential_lifetime model
#' @param ... Additional arguments (ignored)
#' @return A solver function that takes (df, par, ...) and returns a
#'   `fisher_mle` object. The `par` argument is accepted but ignored
#'   since the MLE is computed in closed form.
#' @export
fit.exponential_lifetime <- function(object, ...) {
  function(df, par = NULL, ...) {
    stats <- .exponential_sufficient_stats(df, object)

    if (stats$d == 0) {
      stop("Cannot compute MLE: no exact (uncensored) observations")
    }

    lambda_hat <- stats$d / stats$total_time

    fisher_mle(
      par = c(lambda = lambda_hat),
      loglik_val = stats$d * log(lambda_hat) - lambda_hat * stats$total_time,
      hessian = matrix(-stats$d / lambda_hat^2, 1, 1,
                       dimnames = list("lambda", "lambda")),
      score_val = c(lambda = 0),
      nobs = stats$n,
      converged = TRUE
    )
  }
}


#' Fisher information matrix for exponential_lifetime
#'
#' Returns the analytical FIM: n / lambda^2 (1x1 matrix).
#'
#' @param model An exponential_lifetime model
#' @param ... Additional arguments (ignored)
#' @return A function(theta, n_obs, ...) computing the FIM
#' @export
fim.exponential_lifetime <- function(model, ...) {
  function(theta, n_obs, ...) {
    lambda <- theta[1]
    stopifnot(lambda > 0)
    matrix(n_obs / lambda^2, nrow = 1, ncol = 1,
           dimnames = list("lambda", "lambda"))
  }
}


#' Random data generation for exponential_lifetime
#'
#' Generates exponential data, optionally with right-censoring at a fixed
#' censoring time.
#'
#' @param model An exponential_lifetime model
#' @param ... Additional arguments (ignored)
#' @return A function(theta, n, censor_time, ...) that returns a data frame
#' @importFrom stats rexp
#' @export
rdata.exponential_lifetime <- function(model, ...) {
  function(theta, n, censor_time = Inf, ...) {
    lambda <- theta[1]
    stopifnot(lambda > 0, n > 0)
    x <- rexp(n, rate = lambda)

    no_censoring <- is.null(model$censor_col) || is.infinite(censor_time)
    if (no_censoring) {
      df <- data.frame(x)
      names(df) <- model$ob_col
      return(df)
    }

    df <- data.frame(
      pmin(x, censor_time),
      ifelse(x > censor_time, "right", "exact")
    )
    names(df) <- c(model$ob_col, model$censor_col)
    df
  }
}


#' Compute sufficient statistics for the exponential model
#'
#' @param df Data frame
#' @param model An exponential_lifetime model
#' @return List with d (number of exact obs), total_time (sum of all times),
#'   n (total observations)
#' @keywords internal
.exponential_sufficient_stats <- function(df, model) {
  x <- df[[model$ob_col]]
  n <- length(x)
  stopifnot(n > 0, all(x >= 0))

  if (is.null(model$censor_col)) {
    d <- n
  } else {
    censor <- df[[model$censor_col]]
    exact <- is.na(censor) | censor == "exact"
    d <- sum(exact)
  }

  list(d = d, total_time = sum(x), n = n)
}
