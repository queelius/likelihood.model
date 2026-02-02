#' Fisherian Likelihood Inference
#'
#' Functions for pure likelihood-based inference in the Fisherian tradition.
#' These emphasize the likelihood function itself as the primary object of
#' inference, without requiring probability statements about parameters.
#'
#' @name fisherian
NULL


#' Support Function (Log Relative Likelihood)
#'
#' Computes the support for parameter value theta relative to the MLE:
#' S(theta) = logL(theta) - logL(theta_hat)
#'
#' The support function is always <= 0, with maximum at the MLE.
#' Values of theta with support > -2 are considered well-supported.
#' Values with support > -log(8) ~ -2.08 correspond roughly to a
#' 95% likelihood interval.
#'
#' @param x A fisher_mle object
#' @param theta Parameter value(s) to evaluate
#' @param data Data frame used for likelihood computation
#' @param model The likelihood model used for fitting
#' @param ... Additional arguments passed to loglik
#' @return Support value(s): logL(theta) - logL(theta_hat)
#' @export
support <- function(x, ...) {
  UseMethod("support")
}

#' @rdname support
#' @export
support.fisher_mle <- function(x, theta, data, model, ...) {
  ll <- loglik(model, ...)
  ll(data, theta, ...) - x$loglik
}


#' Relative Likelihood
#'
#' Computes the relative likelihood (likelihood ratio) for theta:
#' R(theta) = L(theta) / L(theta_hat) = exp(S(theta))
#'
#' The relative likelihood is always between 0 and 1, with maximum 1 at the MLE.
#' Common cutoff values:
#' - R >= 0.15 (k=8): roughly equivalent to 95% confidence
#' - R >= 0.10 (k=10): more conservative
#' - R >= 0.05 (k=20): very conservative
#'
#' @param x A fisher_mle object
#' @param theta Parameter value(s) to evaluate
#' @param data Data frame used for likelihood computation
#' @param model The likelihood model used for fitting
#' @param ... Additional arguments passed to loglik
#' @return Relative likelihood value(s): L(theta)/L(theta_hat)
#' @export
relative_likelihood <- function(x, ...) {
  UseMethod("relative_likelihood")
}

#' @rdname relative_likelihood
#' @export
relative_likelihood.fisher_mle <- function(x, theta, data, model, ...) {
  exp(support(x, theta, data, model, ...))
}


#' Likelihood Interval
#'
#' Computes the likelihood interval for a parameter:
#' LI(k) = {theta : R(theta) >= 1/k} = {theta : S(theta) >= -log(k)}
#'
#' Unlike confidence intervals, likelihood intervals make no probability
#' statements about the parameter. They simply identify the set of
#' parameter values that are well-supported by the data.
#'
#' Common choices for k:
#' - k = 8: 1/8 likelihood interval (~95% CI equivalent)
#' - k = 15: 1/15 likelihood interval (~99% CI equivalent)
#' - k = 32: 1/32 likelihood interval (~99.9% CI equivalent)
#'
#' For multivariate parameters, specify `param` to get a profile likelihood
#' interval for that parameter (profiling over the others).
#'
#' @param x A fisher_mle object
#' @param data Data frame used for likelihood computation
#' @param model The likelihood model used for fitting
#' @param k Likelihood ratio cutoff (default 8, giving 1/8 interval)
#' @param param Index or name of parameter for profile interval (NULL for all)
#' @param ... Additional arguments passed to optimization
#' @return Matrix with lower and upper bounds for each parameter
#' @importFrom stats uniroot optim
#' @export
likelihood_interval <- function(x, ...) {
  UseMethod("likelihood_interval")
}

#' @rdname likelihood_interval
#' @export
likelihood_interval.fisher_mle <- function(x, data, model, k = 8,
                                            param = NULL, ...) {
  theta_hat <- coef(x)
  p <- length(theta_hat)
  target_support <- -log(k)
  ll <- loglik(model, ...)
  ll_max <- x$loglik

  if (is.null(param)) {
    param <- seq_len(p)
  } else if (is.character(param)) {
    param <- match(param, names(theta_hat))
  }

  # Get reasonable search bounds from Wald interval
  se_vals <- se(x)
  # Use wider search range than Wald CI
  search_mult <- sqrt(2 * k)

  result <- matrix(NA_real_, nrow = length(param), ncol = 2)
  colnames(result) <- c("lower", "upper")

  for (i in seq_along(param)) {
    idx <- param[i]
    se_i <- se_vals[idx]
    mle_i <- theta_hat[idx]

    if (p == 1) {
      # Univariate case: direct root finding
      support_fn <- function(theta_val) {
        ll(data, theta_val, ...) - ll_max - target_support
      }

      # Find lower bound
      lower_search <- mle_i - search_mult * se_i
      result[i, 1] <- tryCatch({
        uniroot(support_fn, c(lower_search, mle_i),
                extendInt = "yes", tol = 1e-6)$root
      }, error = function(e) NA_real_)

      # Find upper bound
      upper_search <- mle_i + search_mult * se_i
      result[i, 2] <- tryCatch({
        uniroot(support_fn, c(mle_i, upper_search),
                extendInt = "yes", tol = 1e-6)$root
      }, error = function(e) NA_real_)

    } else {
      # Multivariate case: profile likelihood
      # For parameter idx, optimize over all other parameters

      profile_support <- function(theta_i) {
        # Fix parameter idx at theta_i, optimize over others
        other_idx <- setdiff(seq_len(p), idx)
        init_other <- theta_hat[other_idx]

        if (length(other_idx) == 0) {
          # Only one parameter
          theta_full <- theta_i
          names(theta_full) <- names(theta_hat)
          return(ll(data, theta_full, ...) - ll_max)
        }

        # Reparameterize: use log-transform for positive MLE values
        # so the optimizer works in unconstrained space
        pos_mask <- init_other > 0
        init_unconstrained <- init_other
        init_unconstrained[pos_mask] <- log(init_other[pos_mask])

        opt_fn <- function(unconstrained) {
          other_params <- unconstrained
          other_params[pos_mask] <- exp(unconstrained[pos_mask])
          theta_full <- numeric(p)
          theta_full[idx] <- theta_i
          theta_full[other_idx] <- other_params
          names(theta_full) <- names(theta_hat)
          val <- tryCatch(
            suppressWarnings(-ll(data, theta_full, ...)),
            error = function(e) .Machine$double.xmax
          )
          if (!is.finite(val)) return(.Machine$double.xmax)
          val
        }

        opt <- tryCatch(
          suppressWarnings(
            optim(init_unconstrained, opt_fn, method = "BFGS",
                  control = list(maxit = 100))
          ),
          error = function(e) list(value = Inf)
        )

        -opt$value - ll_max
      }

      # Root finding for profile
      profile_root <- function(theta_i) {
        profile_support(theta_i) - target_support
      }

      # Find lower bound
      lower_search <- mle_i - search_mult * se_i
      result[i, 1] <- tryCatch({
        uniroot(profile_root, c(lower_search, mle_i),
                tol = 1e-5)$root
      }, error = function(e) NA_real_)

      # Find upper bound
      upper_search <- mle_i + search_mult * se_i
      result[i, 2] <- tryCatch({
        uniroot(profile_root, c(mle_i, upper_search),
                tol = 1e-5)$root
      }, error = function(e) NA_real_)
    }
  }

  if (!is.null(names(theta_hat))) {
    rownames(result) <- names(theta_hat)[param]
  }

  attr(result, "k") <- k
  attr(result, "relative_likelihood_cutoff") <- 1 / k
  class(result) <- c("likelihood_interval", "matrix")
  result
}

#' Print likelihood interval
#' @param x A likelihood_interval object
#' @param ... Additional arguments (ignored)
#' @return The likelihood_interval object, invisibly
#' @export
print.likelihood_interval <- function(x, ...) {
  k <- attr(x, "k")
  cutoff <- attr(x, "relative_likelihood_cutoff")
  cat(sprintf("1/%d Likelihood Interval (R >= %.3f)\n", k, cutoff))
  cat("-----------------------------------\n")
  print(unclass(x))
  invisible(x)
}


#' Profile Log-Likelihood
#'
#' Computes the profile log-likelihood for a subset of parameters.
#' For each value of the parameters of interest, the remaining
#' (nuisance) parameters are optimized out.
#'
#' The profile likelihood is useful for:
#' - Visualizing the likelihood surface
#' - Computing likelihood intervals
#' - Eliminating nuisance parameters
#'
#' @param x A fisher_mle object
#' @param data Data frame used for likelihood computation
#' @param model The likelihood model used for fitting
#' @param param Index or name of parameter(s) to profile
#' @param grid Optional grid of values to evaluate (vector or matrix)
#' @param n_grid Number of grid points if grid not specified (default 50)
#' @param range_mult Multiplier for grid range based on SE (default 4)
#' @param ... Additional arguments passed to loglik
#' @return A data frame with parameter values and profile log-likelihood
#' @importFrom stats optim
#' @export
profile_loglik <- function(x, ...) {
  UseMethod("profile_loglik")
}

#' @rdname profile_loglik
#' @export
profile_loglik.fisher_mle <- function(x, data, model, param,
                                       grid = NULL, n_grid = 50,
                                       range_mult = 4, ...) {
  theta_hat <- coef(x)
  p <- length(theta_hat)
  ll <- loglik(model, ...)

  if (is.character(param)) {
    param <- match(param, names(theta_hat))
  }

  if (length(param) > 2) {
    stop("Profile likelihood only supported for 1 or 2 parameters")
  }

  se_vals <- se(x)
  nuisance_idx <- setdiff(seq_len(p), param)

  # Evaluate profile log-likelihood at fixed values for the parameters of
  # interest, optimizing over nuisance parameters when present.
  eval_profile <- function(theta_i) {
    if (length(nuisance_idx) == 0) {
      theta_full <- theta_i
      names(theta_full) <- names(theta_hat)
      return(ll(data, theta_full, ...))
    }

    init_nuisance <- theta_hat[nuisance_idx]
    opt_fn <- function(nuisance) {
      theta_full <- numeric(p)
      theta_full[param] <- theta_i
      theta_full[nuisance_idx] <- nuisance
      names(theta_full) <- names(theta_hat)
      -ll(data, theta_full, ...)
    }
    opt <- optim(init_nuisance, opt_fn, method = "BFGS",
                 control = list(maxit = 100))
    -opt$value
  }

  if (length(param) == 1) {
    # 1D profile
    if (is.null(grid)) {
      center <- theta_hat[param]
      spread <- range_mult * se_vals[param]
      grid <- seq(center - spread, center + spread, length.out = n_grid)
    }

    profile_vals <- vapply(grid, eval_profile, numeric(1))

    result <- data.frame(
      param = grid,
      loglik = profile_vals,
      support = profile_vals - x$loglik,
      relative_likelihood = exp(profile_vals - x$loglik)
    )

    if (!is.null(names(theta_hat))) {
      names(result)[1] <- names(theta_hat)[param]
    }

  } else {
    # 2D profile
    if (is.null(grid)) {
      grid1 <- seq(theta_hat[param[1]] - range_mult * se_vals[param[1]],
                   theta_hat[param[1]] + range_mult * se_vals[param[1]],
                   length.out = n_grid)
      grid2 <- seq(theta_hat[param[2]] - range_mult * se_vals[param[2]],
                   theta_hat[param[2]] + range_mult * se_vals[param[2]],
                   length.out = n_grid)
      grid <- expand.grid(param1 = grid1, param2 = grid2)
    }

    profile_vals <- vapply(
      seq_len(nrow(grid)),
      function(i) eval_profile(as.numeric(grid[i, ])),
      numeric(1)
    )

    result <- data.frame(
      grid,
      loglik = profile_vals,
      support = profile_vals - x$loglik,
      relative_likelihood = exp(profile_vals - x$loglik)
    )

    if (!is.null(names(theta_hat))) {
      names(result)[1:2] <- names(theta_hat)[param]
    }
  }

  attr(result, "mle") <- x
  attr(result, "param") <- param
  class(result) <- c("profile_loglik", "data.frame")
  result
}

#' Print profile log-likelihood
#' @param x A profile_loglik object
#' @param ... Additional arguments (ignored)
#' @return The profile_loglik object, invisibly
#' @export
print.profile_loglik <- function(x, ...) {
  param <- attr(x, "param")
  mle <- attr(x, "mle")

  cat("Profile Log-Likelihood\n")
  cat("----------------------\n")
  cat("Parameter(s):", paste(names(coef(mle))[param], collapse = ", "), "\n")
  cat("MLE:", paste(format(coef(mle)[param], digits = 4), collapse = ", "), "\n")
  cat("Grid points:", nrow(x), "\n")
  cat("Max profile loglik:", format(max(x$loglik), digits = 4), "\n")
  invisible(x)
}


#' Deviance for likelihood models
#'
#' Computes the deviance, which is useful for model comparison.
#'
#' When called with a single model, returns -2 * logL (the deviance
#' relative to a saturated model).
#'
#' When called with two models, returns the deviance difference:
#' D = 2 * (logL_full - logL_reduced)
#'
#' Under the null hypothesis that the reduced model is correct,
#' D is asymptotically chi-squared with df = p_full - p_reduced.
#'
#' @param object A fisher_mle object
#' @param null_model Optional reduced/null model for comparison
#' @param ... Additional arguments (ignored)
#' @return Deviance value
#' @export
deviance.fisher_mle <- function(object, null_model = NULL, ...) {
  if (is.null(null_model)) {
    -2 * object$loglik
  } else {
    if (!inherits(null_model, "fisher_mle")) {
      stop("null_model must be a fisher_mle object")
    }
    2 * (object$loglik - null_model$loglik)
  }
}


#' Likelihood-Based Evidence
#'
#' Computes the strength of evidence for theta1 vs theta2:
#' E(theta1, theta2) = logL(theta1) - logL(theta2)
#'
#' Positive values favor theta1, negative values favor theta2.
#'
#' Conventional interpretation (following Royall):
#' - |E| > log(8) ~ 2.08: Strong evidence
#' - |E| > log(32) ~ 3.47: Very strong evidence
#'
#' @param model The likelihood model
#' @param data Data frame for likelihood computation
#' @param theta1 First parameter value
#' @param theta2 Second parameter value
#' @param ... Additional arguments passed to loglik
#' @return Evidence value (log likelihood ratio)
#' @export
evidence <- function(model, data, theta1, theta2, ...) {
  ll <- loglik(model, ...)
  ll(data, theta1, ...) - ll(data, theta2, ...)
}
