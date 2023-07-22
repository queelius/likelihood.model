#' Likelihood ratio test
#' 
#' @param model1 first likelihood model
#' @param model2 second likelihood model
#' @param data data
#' @param par1 parameters for first model
#' @param par2 parameters for second model
#' @param ... additional arguments
#' @return likelihood ratio test
#' @export
lrt <- function(ll1, ll2, data, par1, par2, df = NULL, ...) {
    stopifnot(is_likelihood_model(ll1), is_likelihood_model(ll2))
    # compute the likelihood ratio test
    lrt <- 2 * (ll2(df, par2, ...) - ll1(df, par1, ...))

    if (is.null(df)) {
        df <- length(par2) - length(par1)
    }
    # compute the p-value
    p.value <- pchisq(lrt, df = df, lower.tail = FALSE)
    list(lrt.value = lrt, p.value = p.value, df = df)
}