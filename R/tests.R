#' Likelihood ratio test
#'
#' @param null the likelihood value from a simpler likelihood model nested
#'             within the alternative model
#' @param alt the likelihood value from the more complicated model
#' @param data a data frame
#' @param dof degrees of freedom
#' @param ... additional arguments
#' @return likelihood ratio test
#' @importFrom stats pchisq
#' @export
lrt <- function(null, alt, data, null_par, alt_par, dof = NULL, ...) {
  stopifnot(is_likelihood_model(null), is_likelihood_model(alt),
            is.data.frame(data), dof >= 0)
  stat <- -2 * (loglik(null)(df = data, par = null_par, ...) -
                  loglik(alt)(df = data, par = alt_par, ...))
  if (is.null(dof)) {
    dof <- length(alt_par) - length(null_par)
  }
  p.value <- pchisq(stat, df = dof, lower.tail = FALSE)
  list(stat = stat, p.value = p.value, dof = dof)
}
