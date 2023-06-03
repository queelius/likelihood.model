#' likelihood_model class
#'
#' This class encapsulates all necessary parts of a likelihood model.
#' It provides methods for computing the log-likelihood, score, and Hessian
#' of the model. It also allows for different functions depending on the 
#' observation type.
#' 
#' @field logliks list of functions for computing log-likelihoods
#' @field scores list of functions for computing scores
#' @field hess_logliks list of functions for computing Hessians
#' @field obs_type function that determines observation type
#' @field assumptions list of assumptions made by the model
#' @importFrom R6 R6Class
#' @export likelihood_model
likelihood_model <- R6::R6Class(
    "likelihood_model",
    public = list(
        logliks = NULL,
        scores = NULL,
        hess_logliks = NULL,
        obs_type = NULL,
        assumptions = NULL,

        #' @description
        #' Initializes a new instance of the class
        #' 
        #' @param logliks list of functions for computing log-likelihoods
        #' @param obs_type function that determines observation type
        #' @param scores list of functions for computing scores, default is
        #' NULL (finite difference method is used)
        #' @param hess_logliks list of functions for computing Hessians,
        #' default is NULL (finite difference method is used)
        #' @param assumptions list of assumptions made by the model, default is
        #' c("iid""), which assumes iid observations (this assumption is always
        #' made for this class, which is why we can sum the log-likelihood
        #' contributions)
        #' @return A new `likelihood_model` object
        initialize = function(logliks,
                              obs_type,
                              scores = NULL,
                              hess_logliks = NULL,
                              assumptions = c("iid")) {
            self$logliks <- logliks
            self$scores <- scores
            self$hess_logliks <- hess_logliks
            self$obs_type <- obs_type

            # iid assumption is always made
            self$assumptions <- unique(c("iid", assumptions))
        },

        #' @description
        #' Generic dispatcher for calculating total log-likelihood, score,
        #' and Hessian.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param methods methods to be dispatched
        #' @param ... additional arguments
        #' @return The sum of the dispatched methods
        #' @importFrom magrittr %>%
        dispatch = function(df, par, methods, ...) {
            if (is.null(methods)) {
                return(NULL)
            }
            sapply(1:nrow(df), function(i) {
                otype <- self$obs_type(df[i, ])
                methods[[otype]](df[i, ], par, ...)
            }) %>% sum()
        },

        #' @description
        #' Computes the log-likelihood of the model.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The total log-likelihood
        loglik = function(df, par, ...) {
            self$dispatch(df, par, self$logliks, ...)
        },

        #' @description
        #' Computes the score of the model.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The total score
        score = function(df, par, ...) {
            self$dispatch(df, par, self$scores, ...)
        },

        #' @description
        #' Computes the Hessian of the log-likelihood.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The Hessian of the log-likelihood
        hess_loglik = function(df, par, ...) {
            self$dispatch(df, par, self$hess_logliks, ...)
        }
    )
)

#' Log-likelihood method for likelihood_model
#'
#' This method returns the log-likelihood function for a likelihood_model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the log-likelihood given a data frame and 
#' parameters
#' @export
loglik.likelihood_model <- function(model, ...) {
    function(df, par) model$loglik(df, par, ...)
}

#' Score method for likelihood_model
#'
#' This method returns the score function for a likelihood_model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the score given a data frame and parameters
#' @export
score.likelihood_model <- function(model, ...) {
    function(df, par) model$score(df, par, ...)
}

#' Hessian of log-likelihood method for likelihood_model
#'
#' This method returns the hessian of the log-likelihood function for a 
#' likelihood_model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the hessian of the log-likelihood given a
#' data frame and parameters
#' @export
hess_loglik.likelihood_model <- function(model, ...) {
    function(df, par) model$hess_loglik(df, par, ...)
}

#' Compute an estimate of the Fisher information matrix (FIM) using MC
#' simulation.
#' 
#' This is per-observation. If you want the total FIM, multiply by the
#' number of observations in a sample from the DGP.
#' 
#' @param model likelihood model
#' @param par true parameters
#' @param data a sample from the assumed DGP (larger sample size is better)
#' @param ... additional arguments
#' @return MC estimate of the FIM
#' @export
#' @examples
#' # generate data
#' df <- data.frame(x = rnorm(10000, 0, 1))
#' 
#' # define likelihood model
#' model <- likelihood_model$new(
#'   logliks = list(
#'     normal = function(row, par) dnorm(row$x, par[1], par[2], log = TRUE)
#'   ),
#'   obs_type = function(row) "normal"
#' )
#' 
#' # compute FIM
#' fim(model, c(0, 1), df)
fim.likelihood_model <- function(model, par, data, ...) {

    # compute score. we only get a row of data at a time and apply
    # the score function to it
    fim_mc <- matrix(0, nrow = length(par), ncol = length(par))
    R <- nrow(data)
    for (i in 1:R) {
        # compute score
        s <- model$score(data[i, ], par, ...)

        # compute FIM
        fim_mc <- fim_mc + s %*% t(s)
    }

    fim_mc / R
}

fim <- function(model, par, data, ...) {
    # compute scores
    scores <- apply(data, 1, function(row) model$score(row, par, ...))

    # compute FIM
    tcrossprod(scores) / nrow(data)
}

