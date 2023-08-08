#' @title Exact Normal likelihood model
#' @description
#' A likelihood model for observations that are exact
#' observations of a random variable that follows a Normal
#' distribution. The model assumes that the observations are
#' independent and identically distributed (i.i.d.).
#' 
#' This model defines all the relevant functions to satisfy
#' being a likelihood model. It explicitly provides the following
#' functions:
#' 
#' - `loglik.likelihood_exact_normal`: generates the log-likelihood
#' function for the model.
#' - `score.likelihood_exact_normal`: generates the score function
#' for the model.
#' - `hess_loglik.likelihood_exact_normal`: generates the Hessian
#' of the log-likelihood for the model.
#' 
#' Note: This model may be used as a contribution in the
#' `likelihood_contr_model` if a more complicated likelihood model
#' is needed. We define this simple model mostly for illustrative
#' purposes.
#' 
#' @param ob_col The name of the column in a data frame that
#' contains the observations.
#' @export
likelihood_exact_normal <- function(ob_col) {
    res <- structure(list(ob_col = ob_col))
    class(res) <- c("likelihood_exact_normal", "likelihood_model")
    res
}

#' Customized MLE fit for the exact normal likelihood model.
#'
#' @param model normal likelihood model
#' @param keep_obs whether to store observations in `mle_exact_normal` object
#' @return an `mle` object of class `mle_exact_normal`
#' @importFrom stats var
#' @importFrom MASS ginv
#' @importFrom algebraic.mle mle
#' @importFrom generics fit
#' @export
fit.likelihood_exact_normal <- function(model, keep_obs = TRUE, ...) {
    
    function(df, keep_obs = keep_obs) {

        x <- as.data.frame(df)[[model$ob_col]]
        n <- length(x)
        stopifnot(n > 0)

        theta.hat <- c(
            mu = mean(x),
            var = (1 - 1 / n) * var(x))
        H <- hess_loglik(model, ...)(df, theta.hat)
        algebraic.mle::mle(
            theta.hat = theta.hat,
            loglike = loglik(model, ...)(df, theta.hat),
            score = score(model, ...)(df, theta.hat),
            sigma = ginv(H),
            info = H,
            obs = if (keep_obs) x else NULL,
            nobs = n,
            superclasses = c("mle_exact_normal")
        )
    }
}

#' Log-likelihood function generator for the exact Weibull likelihood model.
#' @param model The weibull likelihood model
#' @param ... Additional arguments (not used)
#' @return A function to compute the log-likelihood given a data frame and
#' parameters.
#' @export
loglik.likelihood_exact_normal <- function(model, ...) {
    function (df, par, ...) {
        stopifnot(!is.null(par), length(par) == 2)
        x <- as.data.frame(df)[[model$ob_col]]
        n <- length(x)
        -n / 2 * log(2 * pi * par[2]) - 1 / (2 * par[2]) *
            (sum(x^2) - 2 * par[1] * sum(x) + n * par[1]^2)
    }
}

#' Score function generator for the exact Weibull likelihood model.
#' @param model The weibull likelihood model
#' @param ... Additional arguments (not used)
#' @export
score.likelihood_exact_normal <- function(model, ...) {
    function(df, par, ...) {
        stopifnot(!is.null(par), length(par) == 2)
        x <- as.data.frame(df)[[model$ob_col]]
        n <- length(x)
        S <- sum(x)
        matrix(c(
            mu = 0.5 / par[2] * (S - n * par[1]),
            var = -0.5 * n / par[2] + 0.5 / par[2]^2 *
                (sum(x^2) - 2 * par[1] * S + n * par[1]^2)),
            nrow = 2)
    }
}

#' Fisher information matrix generator given data `x` for the normal
#' distribution
#'
#' @param x data (simple random sample)
#' @param observed whether to return the observed fisher information, default is
#'                 `TRUE`
#' @param ... additional arguments (not used)
#' @export
hess_loglik.likelihood_exact_normal <- function(model, ...)
    function(df, par, ...) {
        stopifnot(!is.null(par), length(par) == 2)
        x <- as.data.frame(df)[[model$ob_col]]
        n <- length(x)
        S <- sum(x)
        covar <- 0.5 / par[2]^2 * (S - n * par[1])
        matrix(c(
            mu_mu = -n / par[2],
            mu_var = covar, mu_var = covar,
            var_var = -0.5 * n / par[2]^2 + 1 / par[2]^3 *
                (sum(x^2) - 2 * par[1] * S + n * par[1]^2)),
            nrow = 2)
}

#' Computes the bias of an `mle_normal` object.
#' 
#' This function computes the bias of an `mle_normal` object.
#' It is a method for the generic function `bias`.
#' 
#' @param x An `mle_normal` object (subclass of `mle`) to compute the bias of.
#' @param par The true parameter value. If this is unknown (NULL), the bias is
#'            estimated.
#' @param ... Additional arguments (not unused).
#' @return A numeric vector of length 2, the bias of the `mle_normal` estimator.
#' @seealso \code{\link{bias}} for the generic function.
#' @importFrom algebraic.mle bias
#' @importFrom algebraic.dist params
#' @importFrom stats nobs
#' @export
bias.mle_exact_normal <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        c(mu = 0, var = -1 / nobs(x) * params(x)[2])
    } else {
        c(mu = 0, var = -1 / nobs(x) * par[2])
    }
}

#' FIM function generator for the normal distribution with
#' exact observations.
#' 
#' @param model The likelihood model
#' @param ... Additional arguments (not used)
#' @return The FIM function
#' @export
fim.likelihood_exact_normal_model <- function(model, ...) {

    function(par, ...) {
        stopifnot(!is.null(par), length(par) == 2)
        matrix(c(
            mu_mu = 1 / par[2],
            mu_var = 0, mu_var = 0,
            var_var = 1 / (2 * par[2]^2)),
            nrow = 2)
    }
}