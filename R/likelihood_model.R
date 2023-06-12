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


#' Retrieve the assumptions the likelihood model makes about the data.
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A list of assumptions
#' @export
assumptions <- function(model, ...) {
    UseMethod("assumptions")
}

#' Generic method for maximum likelihood estimation (MLE)
#' according to a likelihood model.
#' 
#' @param model The likelihood model
#' @param data The data
#' @param par0 The initial guess for the parameters
#' @param ... Additional arguments
#' @export
mle <- function(model, data, par0, ...) {
    UseMethod("mle")
}

#' Default MLE solver for subclasses of `likelihood_model`.
#' 
#' Note that `likelihood_model` is not a class, but a concept,
#' that other likelihood models implement. They should add
#' `likelihood_model` to their class definition, and then
#' they can use this function to compute the MLE.
#' 
#' This function uses the `optim` function to find the MLE of the
#' parameters of a likelihood model. It uses the `loglik` and `score`
#' methods to compute the log-likelihood and score function, respectively.
#' 
#' @param model The likelihood model
#' @param data The data
#' @param par0 The initial guess for the parameters
#' @param control Control parameters
#' @param ... Additional arguments to pass into the likelihood model's
#'            `loglik` and `score` methods
#' @return The MLE of the parameters
#' @importFrom stats optim
#' @importFrom algebraic.mle mle_numerical
#' @export
mle.likelihood_model <- function(
    model, data, par0, control = list(), ...) {

    stopifnot(!is.null(par0), !is.null(data))

    defaults <- list(
        lower = rep(-Inf, length(par0)),
        upper = rep(Inf, length(par0)),
        method = "Nelder-Mead",
        sann_sampler = NULL,
        hessian = TRUE)
    control <- modifyList(defaults, control)
    ll <- loglik(model)
    s <- NULL
    # which methods use gradient?
    if (control$method %in% c("BFGS", "CG", "L-BFGS-B")) {
        s <- score(model)
    }

    res <- optim(
        par = par0,
        fn = function(par) ll(data, par, ...),
        gr = if (is.null(s)) NULL else function(par) s(data, par, ...),
        method = control$method,
        lower = control$lower,
        upper = control$upper,
        hessian = control$hessian,
        control = list(fnscale=-1, reltol=1e-8, maxit=1000))

    if (is.null(res$hessian)) {
        res$hessian <- hess_loglik(model)(data, res$par, ...)
    }

    mle_numerical(res)
}


#' Likelihood method
#' 
#' This function returns the likelihood function of a model.
#' 
#' @param model The likelihood model
#' @param ... Additional arguments to pass into the `loglik` function
#' @return A likelihood function to compute the likelihood given a data frame
likelihood.likelihood_model <- function(model, ...) {
    ll <- loglik(model, ...)
    function(df, par, ...) {
        exp(ll(df, par, ...))
    }
}

#' Default implementation of the Fisher information matrix (FIM)
#' for a likelihood model
#' 
#' Since we do not assume a particular DGP, and often times it is
#' difficult to compute the FIM analytically, we use the empirical
#' FIM, which is the average of the outer product of the score
#' function over a empirical sample (`data`). The empirical FIM
#' is an estimate of the true FIM.
#' 
#' It has a number of uses, e.g., its inverse is the Cramer-Rao
#' lower bound (CRLB), which is the minimum variance of an unbiased
#' estimator of the parameters. The CRLB has a number of uses,
#' such as computing the asymptotic sampling distribution of the
#' maximum likelihood estimator (MLE) of the parameters for a specified
#' likelihood model, which is asymptotically normal with mean
#' equal to the true parameters and variance equal to the inverse
#' of the FIM. This may be used to compute standard errors of the
#' parameters, confidence intervals, and other quantities.
#' 
#' This is per-observation. If you want the total FIM, multiply by the
#' number of observations for a given data set.
#' 
#' If the true parameters are not known, then you can use an MLE
#' or some other value (say in a hypothesis test) for the parameters.
#' 
#' @param model The likelihood model
#' @param par The true parameters
#' @param data The data
#' @param ... Additional arguments
#' @return The empirical FIM
#' @examples
#' # generate data
#' df <- data.frame(x = rnorm(10000, 0, 1))
#' 
#' # define likelihood model
#' model <- likelihood_contr_model$new(
#'   logliks = list(
#'     normal = function(row, par) dnorm(row$x, par[1], par[2], log = TRUE)
#'   ),
#'   obs_type = function(row) "normal"
#' )
#' 
#' # compute FIM
#' fim(model, c(0, 1), df)
#' @export
fim.likelihood_model <- function(model, par, data, ...) {
    ## scores <- apply(data, 1, function(row) model$score(row, par, ...))
    ## tcrossprod(scores) / nrow(data)
    # fim_mc <- matrix(0, nrow = length(par), ncol = length(par))
    #R <- nrow(data)
    #for (i in 1:R) {
    #    s <- model$score(data[i, ], par, ...)
    #    fim_mc <- fim_mc + s %*% t(s)
    #}
    #fim_mc / R
    s <- score(model)
    Reduce("+",
        lapply(split(data, seq_len(nrow(data))),
            function(row) { s(row, par) }),
        function(s) tcrossprod(as.matrix(s))) / nrow(data)
}
