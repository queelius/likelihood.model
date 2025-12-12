#' Likelihood model generator for standard R distributions
#'
#' Creates a likelihood model based on R's distribution naming convention,
#' where distributions have functions named `d<name>` (PDF), `p<name>` (CDF),
#' `q<name>` (quantile), and `r<name>` (random generation).
#'
#' The model automatically handles exact and censored observations:
#' - Exact: uses PDF (d<name>)
#' - Left-censored: uses CDF (p<name>)
#' - Right-censored: uses survival function (1 - CDF)
#'
#' @param dist_name The name of the distribution (e.g., "norm", "weibull", "exp")
#' @param ob_col The name of the column containing observations
#' @param censor_col The name of the column containing censoring type, or NULL
#'   for all exact observations. Valid values: "exact", "left", "right", or NA.
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
#' @export
likelihood_name <- function(dist_name, ob_col, censor_col = NULL) {
  res <- list(
    dist_name = dist_name,
    pdf = get(paste0("d", dist_name)),
    cdf = get(paste0("p", dist_name)),
    ob_col = ob_col,
    censor_col = censor_col
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

    # Get censoring status (NULL censor_col means all exact)
    if (is.null(model$censor_col)) {
      censor <- rep(NA, n)
    } else {
      censor <- df[[model$censor_col]]
    }

    # Identify observation types
    exact <- is.na(censor) | censor == "exact"
    left <- !is.na(censor) & censor == "left"
    right <- !is.na(censor) & censor == "right"

    # Build parameter list (excluding the first argument which is x)
    par_args <- prepare_args_list(par, model$pdf)

    val <- 0

    # Exact observations: sum of log-PDF
    if (any(exact)) {
      val <- val + sum(do.call(model$pdf, c(list(xs[exact]), par_args, log = TRUE)))
    }

    # Left-censored: sum of log-CDF
    if (any(left)) {
      val <- val + sum(do.call(model$cdf, c(list(xs[left]), par_args, log.p = TRUE)))
    }

    # Right-censored: sum of log-survival
    if (any(right)) {
      val <- val + sum(do.call(model$cdf, c(list(xs[right]), par_args,
                                             log.p = TRUE, lower.tail = FALSE)))
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

  if (is.null(names(par))) {
    # Unnamed: match positionally
    names(par) <- fn_args[seq_along(par)]
  }

  as.list(par)
}
