#' Log-likelihood method
#'
#' This function returns the log-likelihood function of a model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the log-likelihood given a data frame and
#' parameters
#' @export
loglik <- function(model, ...) {
    UseMethod("loglik")
}

#' Score method
#'
#' This function returns the score function of a model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the score given a data frame and parameters
#' @export
score <- function(model, ...) {
    UseMethod("score")
}

#' Hessian of log-likelihood method
#'
#' This function returns the hessian of the log-likelihood function of a model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the hessian of the log-likelihood given a
#' data frame and parameters
#' @export
hess_loglik <- function(model, ...) {
    UseMethod("hess_loglik")
}

#' Compute the Fisher information matrix (FIM), which is
#' an expectation over the data-generating process (DGP) rather than a
#' particular observed sample.
#' 
#' @param model likelihood model
#' @param par true parameters
#' @param ... additional arguments
#' @return FIM
#' @export
fim <- function(model, par, ...) {
    UseMethod("fim")
}
