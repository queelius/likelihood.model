# ==========================================================================
# Re-export generics from algebraic.mle
# ==========================================================================

#' Generics re-exported from algebraic.mle
#'
#' These generics are defined in \pkg{algebraic.mle} and re-exported
#' here so that users of \pkg{likelihood.model} can access them
#' without explicitly loading \pkg{algebraic.mle}.
#'
#' \describe{
#'   \item{\code{\link[algebraic.mle]{se}}}{Standard errors of parameter estimates}
#'   \item{\code{\link[algebraic.mle]{bias}}}{Bias of parameter estimates}
#'   \item{\code{\link[algebraic.mle]{aic}}}{Akaike Information Criterion}
#'   \item{\code{\link[algebraic.mle]{loglik_val}}}{Log-likelihood value at the MLE}
#'   \item{\code{\link[algebraic.mle]{score_val}}}{Score vector at the MLE}
#'   \item{\code{\link[algebraic.mle]{sampler}}}{Sampling distribution of the estimator}
#' }
#'
#' @importFrom algebraic.mle se bias aic loglik_val score_val sampler
#' @aliases se bias aic loglik_val score_val sampler
#' @export se
#' @export bias
#' @export aic
#' @export loglik_val
#' @export score_val
#' @export sampler
#' @name algebraic.mle-reexports
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
#' @return An object of class `c("fisher_mle", "mle")`
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
    class = c("fisher_mle", "mle")
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
      aic = aic(object)
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

#' Extract log-likelihood value from fisher_mle
#'
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return The log-likelihood value at the MLE
#' @export
loglik_val.fisher_mle <- function(x, ...) {
  x$loglik
}

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

#' Compute AIC for fisher_mle
#'
#' Returns the Akaike Information Criterion:
#' AIC = -2 * logL + 2 * k
#' where k is the number of parameters.
#'
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return AIC value
#' @export
aic.fisher_mle <- function(x, ...) {
  -2 * x$loglik + 2 * length(x$par)
}

#' Compute BIC
#'
#' Returns the Bayesian Information Criterion:
#' BIC = -2 * logL + k * log(n)
#' where k is the number of parameters and n is the sample size.
#'
#' @param x A fisher_mle object
#' @param ... Additional arguments (ignored)
#' @return BIC value
#' @export
bic <- function(x, ...) {
  UseMethod("bic")
}

#' @rdname bic
#' @export
bic.fisher_mle <- function(x, ...) {
  if (is.null(x$nobs)) {
    stop("Cannot compute BIC: nobs not available")
  }
  -2 * x$loglik + length(x$par) * log(x$nobs)
}

#' Bias for fisher_mle (asymptotic)
#'
#' Under regularity conditions, asymptotic bias of the MLE is zero.
#'
#' @param x A fisher_mle object
#' @param theta True parameter value (for simulation studies)
#' @param ... Additional arguments (ignored)
#' @return Bias estimate (vector of zeros)
#' @export
bias.fisher_mle <- function(x, theta = NULL, ...) {
  rep(0, length(x$par))
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
#' @return An object of class `c("fisher_boot", "fisher_mle", "mle", "boot")`
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
    class = c("fisher_boot", "fisher_mle", "mle", "boot")
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
#' @importFrom mvtnorm rmvnorm
#' @export
sampler.fisher_mle <- function(x, ...) {
  mu <- x$par
  sigma <- x$vcov

  if (is.null(sigma)) {
    stop("Cannot create sampler: vcov is NULL")
  }

  p <- length(mu)

  function(n) {
    if (p == 1) {
      matrix(rnorm(n, mean = mu, sd = sqrt(sigma)), ncol = 1)
    } else {
      rmvnorm(n, mean = mu, sigma = sigma)
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
