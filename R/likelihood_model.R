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
#' @param object The `likelihood_model` object
#' @param ... Additional arguments to pass into the likelihood model's
#' `loglik`, `score`, and `hess_loglik` constructors.
#' @return An MLE solver (function) that returns an MLE object and accepts as
#' arguments:
#' - `df`: The data frame
#' - `par`: The initial guess for the parameters
#' - `control`: Control parameters for the optimization algorithm
#' - `...`: Additional arguments to pass into the likelihood model's
#'          constructed functions from `loglik`, `score`, and `hess_loglik`.
#' @importFrom algebraic.mle mle_numerical
#' @importFrom stats optim
#' @importFrom utils modifyList
#' @importFrom generics fit
#' @export fit
#' @export
fit.likelihood_model <- function(object, ...) {
  ll <- loglik(object, ...)
  s <- score(object, ...)
  H <- hess_loglik(object, ...)

  function(df, par,
           method = c("Nelder-Mead", "BFGS", "SANN", "CG", "L-BFGS-B", "Brent"),
           ..., control = list()) {
    stopifnot(!is.null(par))
    defaults <- list(fnscale = -1)
    control <- modifyList(defaults, control)
    method <- match.arg(method)

    sol <- optim(
      par = par,
      fn = function(par) ll(df, par, ...),
      gr = if (method == "SANN") {
        NULL
      } else {
        function(par) s(df, par, ...)
      },
      hessian = FALSE,
      method = method,
      control = control
    )

    sol$hessian <- H(df, sol$par, ...)
    mle_numerical(sol, superclasses = "mle_likelihood_model")
  }
}


#' Estimate the sampling distribution of the MLE for a likelihood model.
#'
#' We use the bootstrap method. In other words, we treat the data as an
#' empirical distribution and sample from it to get a new dataset, then
#' we fit the model to that dataset and return the MLE. We do this R
#' times and return the R MLEs.
#'
#' This is the default method, but if you want to use a different method,
#' you should define your own method for your likelihood model.
#'
#' @param model The likelihood model
#' @param ... Additional arguments to pass into the likelihood model
#' @param nthreads The number of threads to use for parallelization
#' @return A function that returns an bootstrapped sampling distribution of an
#' MLE.
#' @importFrom stats optim
#' @importFrom boot boot
#' @importFrom algebraic.mle mle_boot
#' @importFrom algebraic.dist sampler params
#' @export sampler
#' @export
sampler.likelihood_model <- function(model, df, par, ..., nthreads = 1L) {
  solver <- fit(model, ...)
  sol <- solver(df, par, ...)
  function(n, ...) {
    mle_boot(boot(
      data = df,
      statistic = function(df, ind) {
        params(solver(df[ind, , drop = FALSE], par = params(sol), ...))
      },
      R = n,
      parallel = "multicore",
      ncpus = nthreads, ...
    ))
  }
}

#' Print method for likelihood models
#' 
#' @param x A likelihood model
#' @param ... Additional arguments
#' @export
print.likelihood_model <- function(x, show.loglik = FALSE, ...) {
  cat("Likelihood model:", class(x)[1], "\n")
  cat("---------------\n")
  cat("Observation column:", x$ob_col, "\n")
  cat("Censoring column:", x$censor, "\n")
  cat("Assumptions:\n")
  # print the assumption list as a bulleted list
  # strip whitespace from the assumptions
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
