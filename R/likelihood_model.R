#' @name likelihood_model
#' @title Likelihood model
#' @description
#' A general likelihood model framework. When a likelihood model implements
#' this concept, it may be used to compute the log-likelihood, score, hessian
#' of the log-likelihood, Fisher information matrix, and other quantities for
#' the model. This framework is used to compute the standard errors of the
#' parameters of the model, maximum likelihood estimates, the asymptotic
#' sampling distribution of the MLE, and other quantities.
NULL

#' Log-likelihood method
#'
#' This function returns the log-likelihood function of a model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#' @return A log-likelihood function to compute the log-likelihood given a data frame
#' and parameters.
#' @export
loglik <- function(model, ...) {
    UseMethod("loglik")
}

#' Likelihood method
#' 
#' This function returns the likelihood function of a model.
#' 
#' @param model The likelihood model
#' @param ... Additional arguments to pass into the `loglik` function
#' @return A likelihood function to compute the likelihood given a data frame
likelihood <- function(model, ...) {
    ll <- loglik(model, ...)
    function(df, par) {
        exp(ll(df, par))
    }
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

#' Fisher information matrix method
#' 
#' Compute the Fisher information matrix (FIM), which is
#' an expectation over the data-generating process (DGP):
#' 
#'     FIM(theta|x) = E[-loglik_hessian(x)(data,theta)] / n
#'     FIM(theta|x) = E[score(x)(data,theta) %*% t(score(x)(data,theta))] / n
#'
#' where `x` is the likelihood model (say, exact failure times for
#' exponentially distributed lifetimes), `theta` is the parameters
#' of the likelihood model, `loglik_hessian` computes the hessian
#' of the log-likelihood function on the `data` (which is a random
#' quantity that we take the expectation over), and `score` is
#' the score function, the gradient of the log-likelhood function, and
#' `n` is the number of observations in each `data`. The FIM is
#' the negative of the expected hessian of the log-likelihood function
#' (or the expected outer product of the score function).
#' 
#' The FIM is a function of the parameters, and is used to compute
#' the standard errors of the parameters. The FIM is also used to
#' compute the covariance matrix of the parameters, which is used
#' to compute the standard errors of the parameters.
#' 
#' The FIM is also used to compute the Cramer-Rao lower bound (CRLB),
#' which is the inverse of the FIM. The CRLB is the minimum variance
#' of an unbiased estimator of the parameters. The CRLB is used to
#' compute the asymptotic relative efficiency (ARE) of an estimator
#' of the parameters, which is the ratio of the variance of the
#' estimator to the CRLB.
#' 
#' @param model likelihood model
#' @param par true parameters
#' @param ... additional arguments
#' @return FIM
#' @export
fim <- function(model, par, ...) {
    UseMethod("fim")
}
