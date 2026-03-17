#' Likelihood ratio test
#'
#' Computes the likelihood ratio test statistic and p-value for
#' nested models.
#'
#' @param null the likelihood model for the simpler (null) hypothesis nested
#'             within the alternative model
#' @param alt the likelihood model for the more complicated (alternative)
#'            hypothesis
#' @param data a data frame
#' @param null_par parameter values under the null model
#' @param alt_par parameter values under the alternative model
#' @param dof degrees of freedom (computed automatically if NULL)
#' @param ... additional arguments passed to loglik
#' @return An `lrt_result` object with components `stat`, `p.value`, and `dof`
#' @importFrom stats pchisq
#' @export
lrt <- function(null, alt, data, null_par, alt_par, dof = NULL, ...) {
  stopifnot(is_likelihood_model(null), is_likelihood_model(alt),
            is.data.frame(data))
  if (is.null(dof)) {
    dof <- length(alt_par) - length(null_par)
  }
  stopifnot(dof >= 0)
  stat <- -2 * (loglik(null)(df = data, par = null_par, ...) -
                  loglik(alt)(df = data, par = alt_par, ...))
  p.value <- pchisq(stat, df = dof, lower.tail = FALSE)
  structure(
    list(stat = stat, p.value = p.value, dof = dof),
    class = "lrt_result"
  )
}

#' Print method for likelihood ratio test results
#'
#' @param x An lrt_result object
#' @param ... Additional arguments (ignored)
#' @return The lrt_result object, invisibly
#' @export
print.lrt_result <- function(x, ...) {
  cat("Likelihood Ratio Test\n")
  cat("---------------------\n")
  cat("Test statistic:", format(x$stat, digits = 4), "\n")
  cat("Degrees of freedom:", x$dof, "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  invisible(x)
}
