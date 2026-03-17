# ==========================================================================
# Import generics from algebraic.mle for S3 method registration.
# The generics themselves are available to users via Depends: algebraic.mle.
# ==========================================================================
#' @importFrom algebraic.mle se bias score_val observed_fim mse
#' @importFrom algebraic.dist sampler params nparams obs
NULL


#' Maximum Likelihood Estimate (Fisherian)
#'
#' Creates a `fisher_mle` object representing a maximum likelihood estimate
#' with methods for standard inference. This class emphasizes the Fisherian
#' approach to likelihood-based inference.
#'
#' @param par Numeric vector of parameter estimates (may be named)
#' @param vcov Variance-covariance matrix of the estimates
#' @param loglik_val Log-likelihood value at the MLE
#' @param hessian Hessian matrix of the log-likelihood at the MLE
#' @param score_val Optional score vector at the MLE (should be near zero)
#' @param nobs Number of observations used in estimation
#' @param converged Logical indicating if optimization converged
#' @param optim_result Raw result from optim() for diagnostics
#' @return An object of class `c("fisher_mle", "mle_fit")`
#' @importFrom stats coef confint cov printCoefmat qnorm rnorm
#' @export
fisher_mle <- function(par, vcov = NULL, loglik_val, hessian = NULL,
                       score_val = NULL, nobs = NULL, converged = TRUE,
                       optim_result = NULL) {

  if (is.null(vcov) && !is.null(hessian)) {
    info <- -hessian
    vcov <- tryCatch(
      solve(info),
      error = function(e) {
        warning("Hessian not invertible, vcov set to NULL")
        NULL
      }
    )
  }

  structure(
    list(
      par = par,
      vcov = vcov,
      loglik = loglik_val,
      hessian = hessian,
      score = score_val,
      nobs = nobs,
      converged = converged,
      optim = optim_result
    ),
    class = c("fisher_mle", "mle_fit")
  )
}


# --------------------------------------------------------------------------
# Base R generic implementations
# --------------------------------------------------------------------------

#' Extract coefficients from fisher_mle object
#' @param object A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return Named numeric vector of parameter estimates
#' @export
coef.fisher_mle <- function(object, ...) {
  object$par
}

#' Extract variance-covariance matrix from fisher_mle object
#' @param object A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return Variance-covariance matrix
#' @export
vcov.fisher_mle <- function(object, ...) {
  object$vcov
}

#' Extract log-likelihood from fisher_mle object
#' @param object A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return A logLik object
#' @export
logLik.fisher_mle <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- length(object$par)
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' Extract number of observations from fisher_mle object
#' @param object A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return Number of observations
#' @importFrom stats nobs
#' @export
nobs.fisher_mle <- function(object, ...) {
  object$nobs
}

#' Compute confidence intervals for fisher_mle parameters
#'
#' Computes asymptotic Wald confidence intervals based on the
#' estimated variance-covariance matrix.
#'
#' @param object A fisher_mle object
#' @param parm Parameter names or indices (NULL for all)
#' @param level Confidence level (default 0.95)
#' @param ... Additional arguments (ignored)
#' @return Matrix with columns for lower and upper bounds
#' @export
confint.fisher_mle <- function(object, parm = NULL, level = 0.95, ...) {
  cf <- coef(object)
  se <- se(object)

  if (is.null(parm)) {
    parm <- seq_along(cf)
  } else if (is.character(parm)) {
    parm <- match(parm, names(cf))
  }

  a <- (1 - level) / 2
  z <- qnorm(1 - a)

  ci <- cbind(
    cf[parm] - z * se[parm],
    cf[parm] + z * se[parm]
  )

  pct <- paste0(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%")
  colnames(ci) <- pct

  if (!is.null(names(cf))) {
    rownames(ci) <- names(cf)[parm]
  }

  ci
}

#' Print fisher_mle object
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return The fisher_mle object, invisibly
#' @export
print.fisher_mle <- function(x, ...) {
  cat("Maximum Likelihood Estimate (Fisherian)\n")
  cat("----------------------------------------\n")
  cat("Coefficients:\n")
  print(coef(x))
  cat("\nLog-likelihood:", x$loglik, "\n")
  if (!is.null(x$nobs)) {
    cat("Observations:", x$nobs, "\n")
  }
  if (!x$converged) {
    cat("WARNING: Optimization did not converge\n")
  }
  invisible(x)
}

#' Summarize fisher_mle object
#' @param object A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return A summary_fisher_mle object
#' @export
summary.fisher_mle <- function(object, ...) {
  cf <- coef(object)
  se <- se(object)
  ci <- confint(object)

  coef_table <- cbind(
    Estimate = cf,
    `Std. Error` = se,
    ci
  )

  structure(
    list(
      coefficients = coef_table,
      loglik = object$loglik,
      nobs = object$nobs,
      converged = object$converged,
      aic = stats::AIC(object)
    ),
    class = "summary_fisher_mle"
  )
}

#' Print summary of fisher_mle
#' @param x A summary_fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return The summary_fisher_mle object, invisibly
#' @export
print.summary_fisher_mle <- function(x, ...) {
  cat("Maximum Likelihood Estimate (Fisherian)\n")
  cat("----------------------------------------\n\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = 4, signif.stars = FALSE)
  cat("\n")
  cat("Log-likelihood:", format(x$loglik, digits = 4), "\n")
  cat("AIC:", format(x$aic, digits = 4), "\n")
  if (!is.null(x$nobs)) {
    cat("Number of observations:", x$nobs, "\n")
  }
  if (!x$converged) {
    cat("\nWARNING: Optimization did not converge\n")
  }
  invisible(x)
}


# --------------------------------------------------------------------------
# Custom accessors for fisher_mle
# --------------------------------------------------------------------------

#' Extract score vector from fisher_mle
#'
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return The score vector at the MLE (should be near zero)
#' @export
score_val.fisher_mle <- function(x, ...) {
  x$score
}

#' Extract standard errors from fisher_mle
#'
#' Computes standard errors as the square root of the diagonal
#' of the variance-covariance matrix.
#'
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of standard errors, or NA values if vcov is NULL
#' @export
se.fisher_mle <- function(x, ...) {
  if (is.null(x$vcov)) {
    return(rep(NA_real_, length(x$par)))
  }
  sqrt(diag(x$vcov))
}

#' Bias for fisher_mle
#'
#' Estimates the bias of the MLE. Without a model and true parameter value,
#' returns zeros (asymptotic bias is zero under regularity conditions).
#' When both `theta` and `model` are provided, performs a Monte Carlo
#' simulation study to estimate finite-sample bias.
#'
#' @param x A fisher_mle object
#' @param theta True parameter value (for simulation studies)
#' @param model A likelihood model with `rdata` and `fit` methods (optional)
#' @param n_sim Number of Monte Carlo replicates (default 1000)
#' @param ... Additional arguments (ignored)
#' @return Bias estimate vector
#' @export
bias.fisher_mle <- function(x, theta = NULL, model = NULL, n_sim = 1000, ...) {
  p <- length(x$par)

  if (is.null(theta) || is.null(model)) {
    return(rep(0, p))
  }

  if (is.null(x$nobs)) {
    stop("Cannot estimate MC bias: nobs not available in fisher_mle object")
  }

  rdata_fn <- rdata(model)
  fit_fn <- fit(model)
  estimates <- matrix(NA_real_, nrow = n_sim, ncol = p)
  n_ok <- 0L

  for (i in seq_len(n_sim)) {
    est <- tryCatch({
      sim_data <- rdata_fn(theta, n = x$nobs)
      coef(fit_fn(sim_data, par = theta))
    }, error = function(e) NULL)

    if (!is.null(est) && length(est) == p) {
      n_ok <- n_ok + 1L
      estimates[n_ok, ] <- est
    }
  }

  if (n_ok < 10L) {
    warning(sprintf("Only %d of %d MC replicates succeeded; bias estimate unreliable",
                    n_ok, n_sim))
    if (n_ok == 0L) return(rep(NA_real_, p))
  }

  colMeans(estimates[seq_len(n_ok), , drop = FALSE]) - theta
}


# --------------------------------------------------------------------------
# algebraic.mle interface compatibility
#
# fisher_mle inherits from "mle_fit" but uses different field names (e.g.,
# $par instead of $theta.hat, $hessian instead of $info). These methods
# ensure that algebraic.mle generics dispatch correctly to fisher_mle
# objects without falling through to *.mle_fit methods that access the
# wrong fields.
# --------------------------------------------------------------------------

#' Extract parameter estimates from fisher_mle
#'
#' @param x A fisher_mle object
#' @return Named numeric vector of parameter estimates
#' @export
params.fisher_mle <- function(x) {
  x$par
}

#' Number of parameters in fisher_mle
#'
#' @param x A fisher_mle object
#' @return Integer count of parameters
#' @export
nparams.fisher_mle <- function(x) {
  length(x$par)
}

#' Observed Fisher information matrix from fisher_mle
#'
#' Returns the negative Hessian of the log-likelihood evaluated at the
#' MLE, which estimates the Fisher information matrix.
#'
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return A matrix, or NULL if the Hessian was not computed
#' @export
observed_fim.fisher_mle <- function(x, ...) {
  if (is.null(x$hessian)) return(NULL)
  -x$hessian
}

#' Extract observed data from fisher_mle
#'
#' `fisher_mle` objects do not store observed data by design, so this
#' always returns NULL.
#'
#' @param x A fisher_mle object
#' @return Always `NULL`. `fisher_mle` objects do not store the observed data.
#' @export
obs.fisher_mle <- function(x) {
  x$obs  # NULL by design
}

#' Mean squared error for fisher_mle
#'
#' Computes MSE = Var + Bias^2 (scalar) or Vcov + bias %*% t(bias) (matrix).
#' Under regularity conditions, asymptotic bias is zero, so MSE equals the
#' variance-covariance matrix. When `model` is provided, uses Monte Carlo
#' bias estimation via [bias.fisher_mle()].
#'
#' @param x A fisher_mle object
#' @param theta True parameter value (for simulation studies)
#' @param model A likelihood model (optional, enables MC bias estimation)
#' @param n_sim Number of MC replicates for bias estimation (default 1000)
#' @param ... Additional arguments (ignored)
#' @return MSE matrix or scalar
#' @importFrom stats vcov
#' @export
mse.fisher_mle <- function(x, theta = NULL, ..., model = NULL, n_sim = 1000) {
  b <- bias(x, theta, model = model, n_sim = n_sim)
  V <- vcov(x)
  if (length(x$par) == 1L) V + b^2 else V + b %*% t(b)
}


# ==========================================================================
# Bootstrap MLE class: fisher_boot
# ==========================================================================

#' Bootstrap MLE Estimate
#'
#' Creates a `fisher_boot` object representing bootstrap-based inference
#' for maximum likelihood estimates.
#'
#' @param boot_result Result from boot::boot()
#' @param original_mle The original fisher_mle object
#' @return An object of class `c("fisher_boot", "fisher_mle", "mle_fit", "boot")`
#' @export
fisher_boot <- function(boot_result, original_mle) {
  structure(
    list(
      par = boot_result$t0,
      replicates = boot_result$t,
      vcov = cov(boot_result$t),
      loglik = original_mle$loglik,
      hessian = original_mle$hessian,
      nobs = original_mle$nobs,
      R = boot_result$R,
      boot = boot_result,
      converged = TRUE
    ),
    class = c("fisher_boot", "fisher_mle", "mle_fit", "boot")
  )
}

#' Confidence intervals from bootstrap
#'
#' Computes bootstrap confidence intervals using various methods.
#'
#' @param object A fisher_boot object
#' @param parm Parameter names or indices (NULL for all)
#' @param level Confidence level (default 0.95)
#' @param type Type of bootstrap CI: "perc", "bca", "norm", or "basic"
#' @param ... Additional arguments passed to boot::boot.ci
#' @return Matrix with columns for lower and upper bounds
#' @importFrom boot boot.ci
#' @export
confint.fisher_boot <- function(object, parm = NULL, level = 0.95,
                                 type = c("perc", "bca", "norm", "basic"),
                                 ...) {
  type <- match.arg(type)
  cf <- coef(object)

  if (is.null(parm)) {
    parm <- seq_along(cf)
  } else if (is.character(parm)) {
    parm <- match(parm, names(cf))
  }

  ci <- matrix(NA_real_, nrow = length(parm), ncol = 2)
  pct <- paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
                       trim = TRUE, digits = 3), "%")
  colnames(ci) <- pct

  for (i in seq_along(parm)) {
    idx <- parm[i]
    bci <- tryCatch(
      boot.ci(object$boot, conf = level, type = type, index = idx, ...),
      error = function(e) NULL
    )

    if (!is.null(bci)) {
      ci_vals <- switch(type,
        perc = bci$percent[4:5],
        bca = bci$bca[4:5],
        norm = bci$normal[2:3],
        basic = bci$basic[4:5]
      )
      ci[i, ] <- ci_vals
    }
  }

  if (!is.null(names(cf))) {
    rownames(ci) <- names(cf)[parm]
  }

  ci
}

#' Compute bootstrap bias estimate
#'
#' Estimates bias from bootstrap replicates:
#' bias = mean(bootstrap estimates) - original estimate
#'
#' @param x A fisher_boot object
#' @param theta Ignored (for compatibility)
#' @param ... Additional arguments (ignored)
#' @return Bias estimate vector
#' @export
bias.fisher_boot <- function(x, theta = NULL, ...) {
  colMeans(x$replicates) - x$par
}

#' Print fisher_boot object
#' @param x A fisher_boot object
#' @param ... Additional arguments (ignored)
#' @return The fisher_boot object, invisibly
#' @export
print.fisher_boot <- function(x, ...) {
  cat("Bootstrap MLE (Fisherian)\n")
  cat("-------------------------\n")
  cat("Bootstrap replicates:", x$R, "\n\n")
  cat("Coefficients:\n")
  print(coef(x))
  cat("\nBootstrap SE:", format(se(x), digits = 4), "\n")
  cat("Bootstrap bias:", format(bias(x), digits = 4), "\n")
  invisible(x)
}


#' Asymptotic sampler for fisher_mle
#'
#' Returns a function that samples from the asymptotic normal distribution
#' of the MLE: N(theta_hat, vcov).
#'
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return Function that takes n and returns n x p matrix of samples
#' @export
sampler.fisher_mle <- function(x, ...) {
  mu <- x$par
  sigma <- x$vcov

  if (is.null(sigma)) {
    stop("Cannot create sampler: vcov is NULL")
  }

  p <- length(mu)

  if (p > 1 && !requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for multivariate sampling. ",
         "Install it with: install.packages('mvtnorm')")
  }

  function(n) {
    if (p == 1) {
      matrix(rnorm(n, mean = mu, sd = sqrt(sigma)), ncol = 1)
    } else {
      mvtnorm::rmvnorm(n, mean = mu, sigma = sigma)
    }
  }
}

#' Bootstrap sampler for fisher_boot
#'
#' Returns a function that resamples from the bootstrap distribution.
#'
#' @param x A fisher_boot object
#' @param ... Additional arguments (ignored)
#' @return Function that takes n and returns n x p matrix of resampled estimates
#' @export
sampler.fisher_boot <- function(x, ...) {
  replicates <- x$replicates

  function(n) {
    idx <- sample.int(nrow(replicates), size = n, replace = TRUE)
    replicates[idx, , drop = FALSE]
  }
}
