#' Likelihood model
#'
#' A general likelihood model framework. When an object `x` implements
#' the concept of `likelihood_model`, it may be used to compute the
#' standard errors of the parameters of the model, maximum likelihood
#' estimates, the asymptotic sampling distribution of the MLE, and
#' other quantities.
#' 
#' To satisfy the concept of `likelihood_model`, at a minimum the
#' object must implement ``loglik`. For optimal results, it may
#' also provide implementations for `score`, and `hess_loglik` methods.
#' @export
is_likelihood_model <- function(x) {
    inherits(x, "likelihood_model")
}

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
#' @param ... Additional arguments
#' @return A likelihood function to compute the likelihood given a data frame
#' and parameters.
#' @export
lik <- function(model, ...) {
    UseMethod("lik")
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

#' Default score method
#'
#' In case a `score` method is not provided, this function will be used.
#' It computes the score by numerical differentiation of the log-likelihood.
#' @param model The likelihood model
#' @param control Custom arguments to pass to `numDeriv::grad`.
#' @param ... Additional arguments (to pass into `loglik`)
#' @importFrom numDeriv grad
#' @export
score.likelihood_model <- function(
    model,
    control = list(),
    ...) {

    defaults <- list(
        method = "Richardson",
        side = NULL,
        method.args = list(r = 6))
    control <- modifyList(defaults, control)

    ll <- loglik(model, ...)
    function(df, par, ...) {
        stopifnot(!is.null(df), !is.null(par))

        grad(
            func = ll,
            x = par,
            method = control$method,
            side = control$side,
            method.args = control$method.args,
            ...)
    }
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

#' Default method to compute the hessian of the log-likelihood.
#' 
#' In case a `hess_loglik` method is not provided, this function will be
#' used.
#' It computes the hessian of the log-likelikhood using numerical
#' differentiation.
#' 
#' @param model The likelihood model
#' @param ... Additional arguments (to pass into `loglik`)
#' @param control Custom arguments to pass to `numDeriv::hessian`.
#' @importFrom numDeriv hessian
#' @export
hess_loglik.likelihood_model <- function(
    model,
    ...,
    control = list()) {

    defaults <- list(
        method = "Richardson",
        side = NULL,
        method.args = list(r = 6))
    control <- modifyList(defaults, control)

    if (exists("score.", mode = "function", envir = environment(model))) {
        ll <- loglik(model, ...)
    } else {
        ll <- lik(model, ...)
    }
    ll <- loglik(model, ...)
    function(df, par, ...) {
        stopifnot(!is.null(df), !is.null(par))
        hessian(
            func = ll,
            x = par,
            method = control$method,
            side = control$side,
            method.args = control$method.args,
            ...)
    }
}

#' Fisher information matrix method
#' 
#' This function calculates the Fisher information matrix (FIM), an expectation over
#' the data-generating process (DGP). The FIM is a crucial concept in statistics
#' because it provides information about the precision of estimates and the amount
#' of information that data carries about an unknown parameter.
#'
#' FIM is a function of the parameters, and is used to compute the standard errors 
#' of the parameters. It is also used to compute the covariance matrix of the parameters, 
#' which is in turn used to compute standard errors of the parameters.
#' 
#' Additionally, FIM is used to compute the Cramer-Rao lower bound (CRLB), 
#' the inverse of the FIM. CRLB represents the lower limit of the variance that 
#' an unbiased estimator can attain. This is used to compute the asymptotic relative 
#' efficiency (ARE) of an estimator of the parameters, which is the ratio of the 
#' variance of the estimator to the CRLB.
#'
#' The function computes FIM(x)(theta), the FIM of the likelihood model `x`, is based
#' on the following formulas:
#' 
#'     FIM(x)(theta) = E[-loglik_hessian(x)(ob,theta)]
#'     FIM(x)(theta) = E[score(x)(ob,theta) %*% t(score(x)(ob,theta))]
#' 
#' where the expectation is taken with respect to ob ~ DGP. The first formula
#' is the expected hessian of the log-likelihood function, and the second formula
#' is the expected outer product of the score function. The two formulas are
#' equivalent. 
#' 
#' @param model A likelihood model
#' @param ... Additional arguments
#' @return Function that computes the FIM given a sample and a parameter vector
#' @export
fim <- function(model, ...) {
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
#' @param ... Additional arguments
#' @export
fit <- function(model, ...) {
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
#' There are a few interesting options for the `control` argument:
#' 
#' * `method`: The optimization method to use. The default is `Nelder-Mead`,
#' which is a derivative-free method. Other options like include `BFGS`
#' are gradient-based methods, which may be preferable if you provide
#' a `score` function (rather than using the default finite-difference). There
#' is also the `SANN` method, which is a simulated annealing method. This
#' method is useful for multimodal likelihood functions, where the MLE
#' may be sensitive to the initial guess. The `SANN` method is a more global
#' method, but it is slower and may require some tweaking. Regardless,
#' if you do use `SANN`, you should follow it up with a local search
#' method like `Nelder-Mead` to refine the solution.
#' 
#' @param model The likelihood model
#' @param algo The algorithm to use for the MLE solver. It should
#' look like the `stats::optim` function, i.e., take similar arguments
#' (if you don't do anything with them in your algorithm, just accept them
#' as `...` and don't do anything with them) and, especially, return a list with
#' the same named values (it may also return more). Defaults to `stats::optim`.
#' @param ... Additional arguments to pass into the likelihood model's
#' `loglik`, `score`, and `hess_loglik` construction methods
#' @return An MLE solver (function) that returns an MLE object and accepts as
#' arguments:
#' - `df`: The data frame
#' - `par0`: The initial guess for the parameters
#' - `method`: additional optimization parameters to use, if available. `optim`
#'    makes use of this to choose among a variety of optimization algorithms.
#' - `lower`: Lower bounds for the parameters
#' - `upper`: Upper bounds for the parameters
#' - `hessian`: If `TRUE`, then the Hessian of the log-likelihood function
#'    is computed
#' - `control`: Control parameters for the optimization algorithm
#' - `...`: Additional arguments to pass into loglikelihood and score functions
#' 
#' @importFrom algebraic.mle mle_numerical
#' @importFrom stats optim
#' @export
fit.likelihood_model <- function(model, algo = stats::optim,
    use_numerical_hess = FALSE, ...) {

    ll <- loglik(model, ...)
    s <- score(model, ...)
    H <- hess_loglik(model, ...)


    function(
        df,
        par0,
        method = c("Nelder-Mead", "BFGS", "SANN", "CG", "L-BFGS-B", "Brent"),
        lower = -Inf, upper = Inf,
        hessian = TRUE,
        control = list(), ...) {

        stopifnot(!is.null(par0), is.data.frame(df))

        method <- match.arg(method)

        defaults <- list(
            fnscale = -1,
            maxit = 1000L)
        control <- modifyList(defaults, control)

        res <- algo(
            par = par0,
            fn = function(par, ...) ll(df, par, ...),
            gr = function(par, ...) s(df, par, ...),
            ...,
            method = method,
            lower = lower,
            upper = upper,
            hessian = hessian && use_numerical_hess,
            control = control)

        if (hessian && !use_numerical_hess) {
            res$hessian <- H(df, res$par, ...)
        }

        mle_numerical(res, superclasses = "mle_likelihood_model")
    }
}

#' Likelihood method
#' 
#' This function returns the likelihood function of a model.
#' 
#' @param model The likelihood model
#' @param ... Additional arguments to pass into the `loglik` function
#' @return A likelihood function to compute the likelihood given a data frame
lik.likelihood_model <- function(model, ...) {
    ll <- loglik(model, ...)
    function(df, par, ...) {
        exp(ll(df, par, ...))
    }
}

#' Default implementation of the Fisher information matrix (FIM)
#' for a likelihood model
#' 
#' This is an empirical FIM, which is the average of the outer product
#' of the score function over a empirical sample (`data`) or a large simulated
#' sample. The empirical FIM is an estimate of the true FIM.
#' 
#' @param model The likelihood model
#' @param ... Additional arguments to pass into the score and hess_loglik function
#' generators.
#' @return The empirical FIM function, which takes as arguments:
#' - `par`: The parameter vector
#' - `df`: a data frame of observations in the sample used to estimate the FIM. If
#' you can do a Monte-carlo simulation, just take a very large sample and that should
#' give you a very good estimate of the FIM.
#' @export
fim.likelihood_model <- function(model, ...) {

    s <- score(model, ...)
    H <- hess_loglik(model, ...)

    function(par, df, method = "vectorized", ...) {

        stopifnot(is.data.frame(df), !is.null(par))

        p <- length(par)
        R <- nrow(df)

        if (method == "iteration") {
            fim_mc <- matrix(rep(0, p*p), nrow = p, ncol = p)
            for (i in seq_len(R)) {
                si <- s(df[i, ], par, ...)
                fim_mc <- fim_mc + si %*% t(si)
            }
            return(fim_mc / R)
        } else { # if (method == "vectorized") {
            scores <- apply(df, 1, function(row) s(as.data.frame(t(row)), par, ...))
            return(tcrossprod(scores) / R)
        }
    }
}
