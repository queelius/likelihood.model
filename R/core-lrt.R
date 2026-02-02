#' Likelihood ratio test
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
#' @return likelihood ratio test
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
  list(stat = stat, p.value = p.value, dof = dof)
}
