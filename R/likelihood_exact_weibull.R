#' @title Exact Weibull likelihood model
#' @description
#' A likelihood model for observations that are exact
#' observations of a random variable that follows a Weibull
#' distribution. The model assumes that the observations are
#' independent and identically distributed (i.i.d.).
#'
#' This model defines all the relevant functions to satisfy
#' being a likelihood model. It explicitly provides the following
#' functions:
#'
#' - `loglik.likelihood_exact_weibull`: generates the log-likelihood
#' function for the model.
#' - `score.likelihood_exact_weibull`: generates the score function
#' for the model.
#' - `hess_loglik.likelihood_exact_weibull`: generates the Hessian
#' of the log-likelihood for the model.
#'
#' Note: This model may be used as a contribution in the
#' `likelihood_contr_model` if a more complicated likelihood model
#' is needed. For instance, if censoring is needed, then
#' you can also define a censoring likelihood model, or just
#' directly define it in the `likelihood_contr_model`
#' if it is simple enough.
#'
#' @param ob_col The name of the column in a data frame that
#' contains the observations.
#' @export
likelihood_exact_weibull <- function(ob_col) {
  res <- structure(list(ob_col = ob_col))
  class(res) <- c("likelihood_exact_weibull", "likelihood_model")
  res
}

#' List the assumptions made by the model.
#'
#' Since this is a very simple model that assumes only one
#' observation type, the assumptions are simply "independent"
#' and "identically distributed".
#'
#' @param model The weibull likelihood model
#' @export
assumptions.likelihood_exact_weibull <- function(model, ...) {
  c("independent", "identically distributed", "exact observations")
}


#' Log-likelihood function generator for the exact Weibull likelihood model.
#' @param model The weibull likelihood model
#' @param ... Additional arguments (not used)
#' @return A function to compute the log-likelihood given a data frame and
#' parameters.
#' @export
loglik.likelihood_exact_weibull <- function(model, ...) {
  function(df, par, ...) {
    x <- as.data.frame(df)[[model$ob_col]]
    n <- length(x)
    stopifnot(n > 0, all(x > 0))
    k <- par[1] # shape
    l <- par[2] # scale
    n * log(k / l^k) + (k - 1) * sum(log(x)) - sum((x / l)^k)
  }
}

#' Score function generator for the exact Weibull likelihood model.
#' @param model The weibull likelihood model
#' @param ... Additional arguments (not used)
#' @export
score.likelihood_exact_weibull <- function(model, ...) {
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

#' Hessian of the log-likelihood function generator for the weibull
#' likelihood model.
#' @param model The weibull likelihood model
#' @export
hess_loglik.likelihood_exact_weibull <- function(model, ...) {
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
