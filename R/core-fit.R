#' Fit a model
#'
#' Re-exported from \pkg{generics}. See [generics::fit()] for details.
#'
#' @importFrom generics fit
#' @export fit
#' @name fit
NULL

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
#' @importFrom stats optim
#' @importFrom utils modifyList
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

    if (length(par) == 1 && method == "Nelder-Mead") {
      method <- "BFGS"
    }

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

    hessian <- H(df, sol$par, ...)
    score_at_mle <- s(df, sol$par, ...)

    fisher_mle(
      par = sol$par,
      loglik_val = sol$value,
      hessian = hessian,
      score_val = score_at_mle,
      nobs = nrow(df),
      converged = (sol$convergence == 0),
      optim_result = sol
    )
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
#' @param x The likelihood model
#' @param df Data frame to bootstrap from
#' @param par Initial parameter values
#' @param ... Additional arguments to pass into the likelihood model
#' @param nthreads The number of threads to use for parallelization
#' @return A function that returns a bootstrapped sampling distribution of an
#' MLE (fisher_boot object).
#' @importFrom boot boot
#' @export
sampler.likelihood_model <- function(x, df, par, ..., nthreads = 1L) {
  solver <- fit(x, ...)
  sol <- solver(df, par, ...)
  function(n, ...) {
    boot_result <- boot(
      data = df,
      statistic = function(df, ind) {
        coef(solver(df[ind, , drop = FALSE], par = coef(sol), ...))
      },
      R = n,
      parallel = "multicore",
      ncpus = nthreads, ...
    )
    fisher_boot(boot_result, sol)
  }
}
