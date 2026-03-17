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
#'     methods for [coef()], [vcov()], [confint()], [se()], [stats::AIC()], [stats::BIC()],
#'     [summary()].
#'   \item [fit()]: Default MLE solver using [optim()]. Models can override this
#'     with closed-form solutions (see [exponential_lifetime] for an example).
#'   \item Fisherian inference: [support()], [relative_likelihood()],
#'     [likelihood_interval()], [profile_loglik()], [evidence()] -- pure
#'     likelihood-based inference without probability statements.
#'   \item [lrt()]: Likelihood ratio test for nested models.
#' }
#'
#' **Example Implementation**:
#' [exponential_lifetime]: Exponential with right-censoring support.
#' Demonstrates closed-form MLE (no optim needed), analytical FIM,
#' and [rdata()] for Monte Carlo validation.
#'
#' For contribution-based models with heterogeneous observation types,
#' see the companion package \pkg{likelihood.contr}.
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

#' Default loglik method
#'
#' Provides a clear error when a model class does not implement `loglik`.
#'
#' @param model A likelihood model
#' @param ... Additional arguments
#' @return Never returns; always errors.
#' @export
loglik.likelihood_model <- function(model, ...) {
  stop(sprintf(
    "Class '%s' must implement loglik() to satisfy the likelihood_model concept",
    class(model)[1]
  ))
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
#' Computes the Fisher information matrix by Monte Carlo simulation using
#' the negative expected Hessian approach. For each of `n_samples` replicates,
#' generates a single-observation dataset via `rdata`, computes
#' `-hess_loglik(single_obs, theta)`, and averages. The result is
#' `n_obs * I_1(theta)` where `I_1(theta)` is the per-observation FIM.
#'
#' This default requires the model to implement `rdata` and `hess_loglik`
#' (or `loglik`, since `hess_loglik` falls back to numerical differentiation).
#'
#' @param model A likelihood model
#' @param ... Additional arguments passed to rdata
#' @return Function that takes (theta, n_obs, n_samples, ...) and returns FIM matrix
#' @export
fim.likelihood_model <- function(model, ...) {
  hess_fn <- hess_loglik(model, ...)
  rdata_fn <- rdata(model)

  function(theta, n_obs, n_samples = 1000, ...) {
    p <- length(theta)
    fim_sum <- matrix(0, nrow = p, ncol = p)
    n_ok <- 0L

    for (b in seq_len(n_samples)) {
      sample <- tryCatch({
        data_b <- rdata_fn(theta, n = 1, ...)
        H <- -hess_fn(data_b, theta, ...)
        if (any(!is.finite(H))) list(ok = FALSE) else list(ok = TRUE, val = H)
      }, error = function(e) list(ok = FALSE))

      if (sample$ok) {
        fim_sum <- fim_sum + sample$val
        n_ok <- n_ok + 1L
      }
    }

    if (n_ok == 0L) {
      stop("All Monte Carlo samples failed in FIM estimation")
    }
    if (n_ok < n_samples) {
      warning(sprintf("FIM estimation: %d of %d samples failed",
                      n_samples - n_ok, n_samples))
    }

    fim_mat <- n_obs * fim_sum / n_ok
    if (!is.null(names(theta))) {
      dimnames(fim_mat) <- list(names(theta), names(theta))
    }
    fim_mat
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

#' Default rdata method
#'
#' Provides a clear error when a model class does not implement `rdata`.
#' The `rdata` method is required by the default [fim.likelihood_model()]
#' for Monte Carlo FIM estimation.
#'
#' @param model A likelihood model
#' @param ... Additional arguments
#' @return Never returns; always errors.
#' @export
rdata.likelihood_model <- function(model, ...) {
  stop(sprintf(
    "Class '%s' does not implement rdata(); needed by fim() for Monte Carlo estimation",
    class(model)[1]
  ))
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
