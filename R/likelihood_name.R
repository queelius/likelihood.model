#' This is a likelihood model generator based on the convention
#' used in R, where we have a name for a distribution,
#' like `normal`, `weibull`, `exp`, etc. and then
#' we have a set of functions that are named
#' `p<name>`, `d<name>`, `q<name>`, `r<name>`, etc.,
#' where the first letter is the first letter of the
#' distribution name. For instance, for the normal
#' distribution, we have `pnorm`, `dnorm`, `qnorm`,
#' `rnorm`, respectively for the cdf, pdf, quantile
#' function, and random variate generator.
#'
#' This function generates a likelihood model based on
#' a distribution name. The distribution name is
#' used to generate the log-likelihood function,
#' the score function, and the Hessian of the
#' log-likelihood function for either exact (pdf),
#' left-censored (cdf), right-censored (1 - cdf),
#' or interval-censored (cdf(b) - cdf(a)) observations.
#'
#' @param dist_name The name of the distribution.
#' @param ob_col The name of the column in a data frame that
#' contains the observations.
#' @param type_col The name of the column in a data frame that
#' contains the censoring type. The values in this column determine the type of
#' censoring.
#' - `exact` or NA: exact observation
#' - `left`: left-censored observation
#' - `right`: right-censored observation
#' - `interval`: interval-censored observation
#' @export
likelihood_name <- function(dist_name, ob_col, censor_col) {
  res <- structure(list(
    dist_name = dist_name,
    pdf = get(paste0("d", dist_name)),
    cdf = get(paste0("p", dist_name)),
    ob_col = ob_col,
    censor_col = censor_col
  ))

  class(res) <- c(
    paste0("likelihood_name_", dist_name),
    "likelihood_name_model",
    "likelihood_model"
  )
  res
}

#' List the assumptions made by the model.
#'
#' Since this is a very simple model that assumes
#' i.i.d. observations sampled from the named
#' distribution, a censoring mechanism (if any)
#' that is independent of the observations and parameters
#' of the model.
#' @param model The likelihood model
#' @param ... Additional arguments
#' @export
assumptions.likelihood_name_model <- function(model, ...) {
  c(
    "independent",
    "identically distributed",
    paste0(model$dist_name, " distribution"),
    "censoring independent of observations",
    "censoring independent of parameters"
  )
}

#' Log-likelihood function generator for the named likelihood model.
#' @param model The likelihood model
#' @param ... Additional arguments (not used)
#' @return A function to compute the log-likelihood given a data frame and
#' parameters.
#' @export
loglik.likelihood_name_model <- function(model, ...) {
  function(df, par, ...) {
    val <- 0
    censoring <- df[[model$censor_col]]
    xs <- df[[model$ob_col]]
    for (i in seq_along(xs)) {
      if (censoring[i] == "left") {
        args <- prepare_args(xs[i], par, model$cdf)
        val <- val + do.call(model$cdf, c(args, log.p = TRUE, ...))
      } else if (censoring[i] == "right") {
        args <- prepare_args(xs[i], par, model$cdf)
        val <- val + do.call(model$cdf,
                             c(args, log.p = TRUE, lower.tail = FALSE, ...))
      } else {
        args <- prepare_args(xs[i], par, model$pdf)
        val <- val + do.call(model$pdf, c(args, log = TRUE))
      }
    }
    val
  }
}
