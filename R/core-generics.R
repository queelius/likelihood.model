#' @title likelihood.model: Likelihood-Based Inference in the Fisherian Tradition
#' @description
#' The `likelihood.model` package provides a framework for likelihood-based
#' inference. The package is organized in layers:
#'
#' **Core Concept** (`core-generics.R`):
#' The `likelihood_model` "concept" -- an abstract interface that any model
#' can implement. At minimum, implement [loglik()]. Optionally provide
#' [score()] and [hess_loglik()] for analytical derivatives; defaults use
#' numerical differentiation via `numDeriv`.
#'
#' **Core Infrastructure**:
#' \itemize{
#'   \item [fisher_mle] / [fisher_boot]: Result objects from MLE fitting, with
#'     methods for [coef()], [vcov()], [confint()], [se()], [aic()], [bic()],
#'     [summary()].
#'   \item [fit()]: Default MLE solver using [optim()]. Models can override this
#'     with closed-form solutions (see [exponential_lifetime] for an example).
#'   \item Fisherian inference: [support()], [relative_likelihood()],
#'     [likelihood_interval()], [profile_loglik()], [evidence()] -- pure
#'     likelihood-based inference without probability statements.
#'   \item [lrt()]: Likelihood ratio test for nested models.
#' }
#'
#' **Model Builders**:
#' \itemize{
#'   \item [likelihood_contr_model]: R6 class for building models from
#'     heterogeneous observation types (exact, censored, etc.) with
#'     dynamic dispatch to type-specific functions.
#'   \item [likelihood_name()]: Wraps any standard R distribution (norm, weibull,
#'     exp, ...) with automatic censoring support.
#' }
#'
#' **Example Implementations**:
#' Reference implementations showing how to satisfy the `likelihood_model`
#' concept with hand-derived analytical solutions:
#' \itemize{
#'   \item [weibull_uncensored]: Weibull with exact observations only.
#'     Demonstrates analytical score and hessian (10-100x faster than numerical).
#'   \item [exponential_lifetime]: Exponential with right-censoring support.
#'     Demonstrates closed-form MLE (no optim needed), analytical FIM,
#'     and [rdata()] for Monte Carlo validation.
#' }
#'
"_PACKAGE"

#' Check if an object is a likelihood model
#'
#' Tests whether an object satisfies the `likelihood_model` concept by
#' checking if `"likelihood_model"` is in its class hierarchy. To be a
#' likelihood model, at a minimum the object must implement [loglik()].
#' For optimal results, it may also provide [score()] and [hess_loglik()].
#'
#' @param x An object to test
#' @return Logical indicating whether x is a likelihood model
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
#' @return A log-likelihood function to compute the log-likelihood given a data
#' frame and parameters.
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

#' Default score method
#'
#' In case a `score` method is not provided, this function will be used.
#' It computes the score by numerical differentiation of the log-likelihood.
#' @param model The likelihood model
#' @param control Custom arguments to pass to `numDeriv::grad`.
#' @param ... Additional arguments (to pass into `loglik`)
#' @return A function to compute the score given a data frame and parameters
#' @importFrom numDeriv grad
#' @export
score.likelihood_model <- function(model, control = list(), ...) {
  defaults <- list(
    method = "Richardson",
    side = NULL,
    method.args = list(r = 6)
  )
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
      df = df,
      ...
    )
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
#' It computes the hessian of the log-likelihood function using numerical
#' differentiation.
#'
#' @param model The likelihood model
#' @param ... Additional arguments (to pass into `loglik`)
#' @param control Custom arguments to pass to `numDeriv::hessian`.
#' @return A function(df, par, ...) that computes the Hessian matrix of the
#'   log-likelihood evaluated at `par` given data `df`.
#' @importFrom numDeriv hessian
#' @export
hess_loglik.likelihood_model <- function(model, control = list(), ...) {
  defaults <- list(
    method = "Richardson",
    side = NULL,
    method.args = list(r = 6)
  )
  control <- modifyList(defaults, control)
  ll <- loglik(model, ...)
  function(df, par, ...) {
    stopifnot(!is.null(df), !is.null(par))
    hessian(
      func = ll,
      x = par,
      method = control$method,
      side = control$side,
      method.args = control$method.args,
      df = df,
      ...
    )
  }
}

#' Fisher information matrix method
#'
#' This function calculates the Fisher information matrix (FIM), an expectation
#' over the data-generating process (DGP). The FIM is a crucial concept in
#' statistics because it provides information about the precision of estimates
#' and the amount of information that data carries about an unknown parameter.
#'
#' FIM is a function of the parameters, and is used to compute the standard
#' errors of the parameters. It is also used to compute the covariance matrix of
#' the parameters, which is in turn used to compute standard errors of the
#' parameters.
#'
#' Additionally, FIM is used to compute the Cramer-Rao lower bound (CRLB),
#' the inverse of the FIM. CRLB represents the lower limit of the variance that
#' an unbiased estimator can attain. This is used to compute the asymptotic
#' relative efficiency (ARE) of an estimator of the parameters, which is the
#' ratio of the variance of the estimator to the CRLB.
#'
#' The function computes FIM(x)(theta), the FIM of the likelihood model `x`, is
#' based on the following formulas:
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
#' @return Function that computes the FIM given a parameter vector and sample size
#' @export
fim <- function(model, ...) {
  UseMethod("fim")
}

#' Default FIM method using Monte Carlo estimation
#'
#' Computes the Fisher information matrix by Monte Carlo simulation.
#' Uses the outer product of scores.
#'
#' This default requires the model to implement `rdata` and `score` methods.
#'
#' @param model A likelihood model
#' @param ... Additional arguments passed to rdata
#' @return Function that takes (theta, n_obs, n_samples, ...) and returns FIM matrix
#' @export
fim.likelihood_model <- function(model, ...) {
  score_fn <- score(model)
  rdata_fn <- rdata(model)

  function(theta, n_obs, n_samples = 1000, ...) {
    p <- length(theta)
    fim_sum <- matrix(0, nrow = p, ncol = p)

    for (b in seq_len(n_samples)) {
      data_b <- rdata_fn(theta, n = n_obs, ...)
      s <- score_fn(data_b, theta)
      fim_sum <- fim_sum + outer(s, s)
    }
    fim_sum / n_samples
  }
}

#' Observed information matrix method
#'
#' Returns the observed information matrix, which is the negative Hessian
#' of the log-likelihood evaluated at the data and parameter values.
#'
#' At the MLE, the observed information is a consistent estimator of the
#' Fisher information matrix.
#'
#' @param model A likelihood model
#' @param ... Additional arguments
#' @return Function that computes -hess_loglik(df, par)
#' @export
observed_info <- function(model, ...) {
  UseMethod("observed_info")
}

#' Default observed information method
#'
#' Returns a function that computes -hess_loglik(df, par).
#'
#' @param model A likelihood model
#' @param ... Additional arguments passed to hess_loglik
#' @return Function that takes (df, par, ...) and returns observed information matrix
#' @export
observed_info.likelihood_model <- function(model, ...) {
  H_fn <- hess_loglik(model, ...)
  function(df, par, ...) {
    -H_fn(df, par, ...)
  }
}

#' Random data generation method
#'
#' Returns a function that generates random data from the model's
#' data-generating process (DGP) at a given parameter value.
#'
#' This is used by the default `fim` method for Monte Carlo estimation
#' of the Fisher information matrix.
#'
#' @param model A likelihood model
#' @param ... Additional arguments
#' @return Function that takes (theta, n, ...) and returns a data frame
#' @export
rdata <- function(model, ...) {
  UseMethod("rdata")
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

#' Print method for likelihood models
#'
#' @param x A likelihood model
#' @param show.loglik Logical; if TRUE, print the log-likelihood function
#' @param ... Additional arguments (ignored)
#' @return The likelihood model object, invisibly
#' @export
print.likelihood_model <- function(x, show.loglik = FALSE, ...) {
  cat("Likelihood model:", class(x)[1], "\n")
  cat("---------------\n")
  if (!is.null(x$ob_col)) cat("Observation column:", x$ob_col, "\n")
  if (!is.null(x$censor_col)) cat("Censoring column:", x$censor_col, "\n")
  cat("Assumptions:\n")
  for (assumption in assumptions(x)) {
    cat(" -", trimws(assumption), "\n")
  }
  if (show.loglik) {
    cat("---------------\n")
    cat("Log-likelihood function:\n")
    print(loglik(x))
  }
  invisible(x)
}
