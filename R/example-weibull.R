#' @title Weibull likelihood model (uncensored)
#' @description
#' A likelihood model for exact (uncensored) observations from a Weibull
#' distribution. The model assumes that the observations are independent
#' and identically distributed (i.i.d.).
#'
#' This is a reference implementation demonstrating how to satisfy the
#' likelihood_model concept with analytical derivatives. It provides:
#'
#' - `loglik.weibull_uncensored`: log-likelihood function
#' - `score.weibull_uncensored`: score (gradient) function
#' - `hess_loglik.weibull_uncensored`: Hessian of the log-likelihood
#'
#' Analytical derivatives are 10-100x faster than the default numerical
#' differentiation via numDeriv.
#'
#' This model may also be used as a contribution in
#' `likelihood_contr_model` for more complex models involving censoring.
#'
#' @param ob_col The name of the column in a data frame that
#' contains the observations.
#' @return A `weibull_uncensored` likelihood model object
#' @export
weibull_uncensored <- function(ob_col) {
  structure(
    list(ob_col = ob_col),
    class = c("weibull_uncensored", "likelihood_model")
  )
}

#' @rdname weibull_uncensored
#' @export
likelihood_exact_weibull <- weibull_uncensored

#' List the assumptions made by the model.
#'
#' Since this is a very simple model that assumes only one
#' observation type, the assumptions are simply "independent"
#' and "identically distributed".
#'
#' @param model The weibull likelihood model
#' @param ... Additional arguments (ignored)
#' @return Character vector of assumptions
#' @export
assumptions.weibull_uncensored <- function(model, ...) {
  c("independent", "identically distributed", "exact observations")
}


#' Log-likelihood function generator for the Weibull uncensored model.
#' @param model The weibull likelihood model
#' @param ... Additional arguments (not used)
#' @return A function to compute the log-likelihood given a data frame and
#' parameters.
#' @export
loglik.weibull_uncensored <- function(model, ...) {
  function(df, par, ...) {
    x <- as.data.frame(df)[[model$ob_col]]
    n <- length(x)
    stopifnot(n > 0, all(x > 0))
    k <- par[1] # shape
    l <- par[2] # scale
    n * log(k / l^k) + (k - 1) * sum(log(x)) - sum((x / l)^k)
  }
}

#' Score function generator for the Weibull uncensored model.
#' @param model The weibull likelihood model
#' @param ... Additional arguments (not used)
#' @return A function to compute the score given a data frame and parameters
#' @export
score.weibull_uncensored <- function(model, ...) {
  function(df, par, ...) {
    x <- as.data.frame(df)[[model$ob_col]]
    n <- length(x)
    stopifnot(n > 0, all(x > 0))
    k <- par[1] # shape
    l <- par[2] # scale
    x_l <- x / l
    x_l_k <- x_l^k
    c(
      shape = sum(log(x)) + n / k - n * log(l) - sum(x_l_k * log(x_l)),
      scale = k * sum(x_l_k) / l - n * k / l
    )
  }
}

#' Hessian of the log-likelihood function generator for the Weibull
#' uncensored model.
#' @param model The weibull likelihood model
#' @param ... Additional arguments (ignored)
#' @return A function to compute the Hessian of the log-likelihood given
#' a data frame and parameters
#' @export
hess_loglik.weibull_uncensored <- function(model, ...) {
  function(df, par, ...) {
    x <- as.data.frame(df)[[model$ob_col]]
    n <- length(x)
    stopifnot(n > 0, all(x > 0))
    k <- par[1] # shape
    l <- par[2] # scale
    x_l_k <- (x / l)^k
    log_x_l <- log(x / l)
    # d^2 ll / dk dl
    p_k_l <- -n / l + sum(x_l_k * (1 + k * log_x_l)) / l
    matrix(
      c(
        # d^2 ll / dk^2 = -n/k^2 - sum((x/l)^k * (log(x/l))^2)
        shape_shape = -n / k^2 - sum(x_l_k * log_x_l^2),
        shape_scale = p_k_l, shape_scale = p_k_l,
        # d^2 ll / dl^2 = nk/l^2 - k(k+1)/l^2 * sum((x/l)^k)
        scale_scale = n * k / l^2 - k * (k + 1) / l^2 * sum(x_l_k)
      ),
      nrow = 2
    )
  }
}
