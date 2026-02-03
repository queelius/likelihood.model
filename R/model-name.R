#' Likelihood model generator for standard R distributions
#'
#' Creates a likelihood model based on R's distribution naming convention,
#' where distributions have functions named \code{d<name>} (PDF), \code{p<name>} (CDF),
#' \code{q<name>} (quantile), and \code{r<name>} (random generation).
#'
#' The model automatically handles exact and censored observations:
#' - Exact: uses PDF (\code{d<name>})
#' - Left-censored: uses CDF (\code{p<name>})
#' - Right-censored: uses survival function (1 - CDF)
#'
#' @param dist_name The name of the distribution (e.g., "norm", "weibull", "exp")
#' @param ob_col The name of the column containing observations (lower bound
#'   for interval-censored data)
#' @param censor_col The name of the column containing censoring type, or NULL
#'   for all exact observations. Valid values: "exact", "left", "right",
#'   "interval", or NA.
#' @param ob_col_upper The name of the column containing the upper bound for
#'   interval-censored observations, or NULL if no interval censoring is used.
#' @return A likelihood model object
#' @examples
#' # Simple case - all exact observations
#' model <- likelihood_name("norm", ob_col = "x")
#' df <- data.frame(x = rnorm(100))
#' ll <- loglik(model)
#' ll(df, c(mean = 0, sd = 1))
#'
#' # With censoring
#' model <- likelihood_name("weibull", ob_col = "time", censor_col = "status")
#' df <- data.frame(time = rweibull(100, 2, 1),
#'                  status = sample(c("exact", "right"), 100, replace = TRUE))
#'
#' # With interval censoring
#' model <- likelihood_name("norm", ob_col = "x", censor_col = "censor",
#'                          ob_col_upper = "x_upper")
#' @export
likelihood_name <- function(dist_name, ob_col, censor_col = NULL,
                            ob_col_upper = NULL) {
  stopifnot(is.character(dist_name), length(dist_name) == 1)
  stopifnot(is.character(ob_col), length(ob_col) == 1)

  pdf_name <- paste0("d", dist_name)
  cdf_name <- paste0("p", dist_name)
  if (!exists(pdf_name, mode = "function")) {
    stop(sprintf("Distribution '%s' not found: no function '%s'",
                 dist_name, pdf_name))
  }
  if (!exists(cdf_name, mode = "function")) {
    stop(sprintf("Distribution '%s' not found: no function '%s'",
                 dist_name, cdf_name))
  }

  res <- list(
    dist_name = dist_name,
    pdf = get(pdf_name, mode = "function"),
    cdf = get(cdf_name, mode = "function"),
    ob_col = ob_col,
    censor_col = censor_col,
    ob_col_upper = ob_col_upper
  )

  class(res) <- c(
    paste0("likelihood_name_", dist_name),
    "likelihood_name_model",
    "likelihood_model"
  )
  res
}

#' List the assumptions made by the model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#' @return Character vector of assumptions
#' @export
assumptions.likelihood_name_model <- function(model, ...) {
  assumptions <- c(
    "independent observations",
    "identically distributed",
    paste0(model$dist_name, " distribution")
  )
  if (!is.null(model$censor_col)) {
    assumptions <- c(assumptions,
      "censoring independent of observations",
      "censoring independent of parameters"
    )
  }
  assumptions
}

#' Log-likelihood function for named distribution models
#'
#' @param model The likelihood model
#' @param ... Additional arguments (not used)
#' @return A function(df, par, ...) that computes the log-likelihood
#' @export
loglik.likelihood_name_model <- function(model, ...) {
  function(df, par, ...) {
    xs <- df[[model$ob_col]]
    n <- length(xs)
    if (n == 0) return(0)

    if (is.null(model$censor_col)) {
      censor <- rep(NA, n)
    } else {
      censor <- df[[model$censor_col]]
    }

    exact <- is.na(censor) | censor == "exact"
    left <- !is.na(censor) & censor == "left"
    right <- !is.na(censor) & censor == "right"
    interval <- !is.na(censor) & censor == "interval"

    # Warn on unrecognized censoring values
    unhandled <- !exact & !left & !right & !interval
    if (any(unhandled)) {
      bad_vals <- unique(censor[unhandled])
      warning(sprintf(
        "Unknown censoring values ignored: %s. Use 'exact', 'left', 'right', 'interval', or NA.",
        paste(sQuote(bad_vals), collapse = ", ")
      ))
    }

    par_args <- prepare_args_list(par, model$pdf)

    val <- 0

    if (any(exact)) {
      val <- val + sum(do.call(model$pdf, c(list(xs[exact]), par_args,
                                             log = TRUE)))
    }

    if (any(left)) {
      val <- val + sum(do.call(model$cdf, c(list(xs[left]), par_args,
                                             log.p = TRUE)))
    }

    if (any(right)) {
      val <- val + sum(do.call(model$cdf, c(list(xs[right]), par_args,
                                             log.p = TRUE, lower.tail = FALSE)))
    }

    if (any(interval)) {
      if (is.null(model$ob_col_upper)) {
        stop("Interval censoring requires 'ob_col_upper' in the model")
      }
      xs_upper <- df[[model$ob_col_upper]]
      val <- val + sum(log(
        do.call(model$cdf, c(list(xs_upper[interval]), par_args)) -
        do.call(model$cdf, c(list(xs[interval]), par_args))
      ))
    }

    val
  }
}

#' Prepare parameter arguments for distribution function calls
#'
#' Converts a parameter vector to a named list suitable for do.call(),
#' matching parameter names to function formals.
#'
#' @param par Parameter vector (named or unnamed)
#' @param fn The distribution function to match against
#' @return Named list of parameters
#' @keywords internal
prepare_args_list <- function(par, fn) {
  # Get function formals, excluding first arg (x/q) and special args
  fn_args <- names(formals(fn))
  fn_args <- fn_args[!fn_args %in% c("x", "q", "p", "n", "log", "log.p",
                                      "lower.tail", "...")]

  if (length(par) > length(fn_args)) {
    stop(sprintf(
      "Too many parameters (%d) for distribution (expects at most %d: %s)",
      length(par), length(fn_args), paste(fn_args, collapse = ", ")
    ))
  }

  if (is.null(names(par))) {
    names(par) <- fn_args[seq_along(par)]
  }

  as.list(par)
}
