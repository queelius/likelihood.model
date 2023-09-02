#' Likelihood ratio test
#' 
#' @param null the null (simpler) likelihood model
#' @param alt the more complicated model
#' @param data the data which the likelihood models are fitted to
#' @param null_par parameters for null model
#' @param alt_par parameters for alternative model
#' @param dof degrees of freedom
#' @param ... additional arguments
#' @return likelihood ratio test
#' @export
lrt <- function(model1, model2, data, par1, par2, dof = NULL, ...) {
    stopifnot(is_likelihood_model(null), is_likelihood_model(alt))
    # compute the likelihood ratio test
    stat <- -2 * (loglik(null)(df = data, par = null_par, ...) - loglik(alt)(df = data, par = alt_par, ...))

    if (is.null(dof)) {
        dof <- length(alt_par) - length(null_par)
    }
    # compute the p-value
    p.value <- pchisq(stat, df = dof, lower.tail = FALSE)
    res <- list(stat = stat, p.value = p.value, dof = dof)
    class(res) <- c("likelihood_ratio_test", "hypothesis_test")
    res
}


#' Hypothesis test structure
#' 
#' @param stat test statistic
#' @param p.value p-value
#' @param dof degrees of freedom
#' @return hypothesis test
#' @export
#' @examples
#' # create a hypothesis test
#' test <- hypothesis_test(stat = 1.96, p.value = 0.05, dof = 1)
#' # print the test
#' test
#' # extract the p-value
#' pval(test)
#' # extract the degrees of freedom
#' dof(test)
#' # extract the test statistic
#' stat(test)
#' # check if the test is significant at the 5% level
#' is_significant_at(test, 0.05)
#' @importFrom stats pf pchisq qt qnorm
#' @export
hypothesis_test <- function(stat, p.value, dof, superclasses = NULL, ...) {
    res <- list(stat = stat, p.value = p.value, dof = dof, ...)
    class(res) <- unique(c(superclasses, "hypothesis_test"))
    res
}

#' Generic method for extracting the p-value from a hypothesis test
#' @param x a hypothesis test object
#' @param ... additional arguments to pass into the method
#' @return p-value
#' @export
pval <- function(x, ...) {
    UseMethod("pval")
}


#' p-value method for hypothesis tests
#' 
#' @param x a hypothesis test
#' @param ... additional arguments
#' @return p-value
#' @export
pval.hypothesis_test <- function(x, ...) {
    x$p.value
}

#' Generic method for extracting the degrees of freedom from a hypothesis test
#' @param x a hypothesis test object
#' @param ... additional arguments to pass into the method
#' @return degrees of freedom
#' @export
dof <- function(x, ...) {
    UseMethod("dof")
}

#' degrees of freedom method for hypothesis tests
#' 