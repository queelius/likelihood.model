#' This is a likelihood mode generator based on the convention
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
#' contains the censoring type. If `NULL`, then
#' no censoring is assumed (exact observations).
#' If `type_col` is not `NULL`, then values in this column
#' are used to determine the type of censoring.
#' - NA or `exact`: exact observation
#' - `left`: left-censored observation
#' - `right`: right-censored observation
#' - `interval`: interval-censored observation
#' @export
likelihood_name <- function(dist_name, ob_col, type_col = NULL) {
    if (is.null(type_col)) {
        res <- structure(list(dist_name = dist_name, ob_col = ob_col))
        class(res) <- c(
            paste0("likelihood_", dist_name, "_model"),
            "likelihood_model")
        res
    } else {
        res <- structure(list(dist_name = dist_name, ob_col = ob_col, type_col = type_col))
        class(res) <- c(
            paste0("likelihood_censored_", dist_name, "_model"),
            "likelihood_model")
        res
    }
}

#' List the assumptions made by the model.
#' 
#' Since this is a very simple model that assumes 
#' i.i.d. observations sampled from the named
#' distribution, a censoring mechanism (if any)
#' that is independent of the observations and parameters
#' of the model.
#' @param model The likelihood model
#' @export 
assumptions.likelihood_name <- function(model, ...) {
    c("independent", "identically distributed",
      paste0(model$dist_name, " distribution"), "censoring independent of observations",
      "censoring independent of parameters")
}

#' Log-likelihood function generator for the named likelihood model.
#' @param model The likelihood model
#' @param ... Additional arguments (not used)
#' @return A function to compute the log-likelihood given a data frame and
#' parameters.
#' @export
loglik.likelihood_name <- function(model, ...) {
    x <- as.data.frame(df)[[model$ob_col]]
    n <- length(x)
    stopifnot(n > 0)

    if (is.null(model$type_col)) {
        function(df, par, ...) {
            do.call(paste0("d", model$dist_name), c(list(x), par))
            sum(do.call(paste0("d", model$dist_name), c(list(x), par), log.p = TRUE))
        }
    } else {
        function(df, par, ...) {
            type <- as.character(df[[model$type_col]])
            ifelse(type == "left",
                   do.call(paste0("p", model$dist_name), c(list(x), par)),
                   ifelse(type == "right",
                          1 - do.call(paste0("p", model$dist_name), c(list(x), par)),
                          ifelse(type == "interval",
                                 do.call(paste0("p", model$dist_name), c(list(x[, 2]), par)) -
                                     do.call(paste0("p", model$dist_name), c(list(x[, 1]), par)),
                                 do.call(paste0("d", model$dist_name), c(list(x), par)))))
            sum(do.call(paste0("d", model$dist_name), c(list(x), par), log.p = TRUE))
        }
    }
}