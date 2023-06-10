#' @title Likelihood contribution model
#' @description
#' This class encapsulates all necessary parts of a likelihood model.
#' A likelihood_model should provide the following methods:
#' 
#' - `loglik`: computes the log-likelihood of the model
#' - `score`: computes the score of the model
#' - `hess_loglik`: computes the Hessian of the log-likelihood given
#'  a data frame and parameters
#' - `fim`: computes the Fisher information matrix (expectation)
#' 
#' It provides methods for computing the log-likelihood, score, the
#' Hessian of log-likelihood, and the Fisher information matrix (FIM).
#' 
#' It also allows for different *likelihood contributions* depending
#' on the  observation type of a row in the data frame. For example,
#' if the data frame contains both exact and interval-censored
#' observations, then the log-likelihood contributions for exact
#' observations and interval-censored observations are computed
#' separately and then summed up (assumes i.i.d. samples).
#'  
#' @field logliks list of functions for computing log-likelihoods
#' @field scores list of functions for computing scores
#' @field hess_logliks list of functions for computing Hessians
#' @field obs_type function that determines observation type
#' @field assumptions list of assumptions made by the model
#' @importFrom R6 R6Class
#' @importFrom magrittr %>%
#' @export likelihood_contr_model
likelihood_contr_model <- R6::R6Class(
    "likelihood_contr_model",
    public = list(
        obs_type = NULL,
        logliks = NULL,
        scores = NULL,
        hess_logliks = NULL,
        assumptions = NULL,

        #' @description
        #' Initializes a new instance of the class
        #' 
        #' @param obs_type function that determines observation type
        #' @param logliks list of functions for computing log-likelihoods
        #' If an observation type does not have a log-likelihood dispatcher
        #' specified here, then we lazily construct a log-likelihood for said
        #' observation type by looking for a `loglik.type`, where `type` is the
        #' observation type. If we cannot find a `loglik.type`, then we throw
        #' an error.
        #' @param scores list of functions for computing scores.
        #' If an observation type does not have a score dispatcher specified
        #' here, then we lazily construct a score for said observation
        #' type using the following method:
        #'  (1) first, we look for a `score.type`, where `type` is the
        #' observation type.
        #'  (2) if (1) fails, then we use a finite difference method given the
        #' log-likelihood contribution for the observation type.
        #' @param hess_logliks list of functions for computing Hessians of the
        #' log-likelihood given the observed data. If an observation type does
        #' not have a Hessian dispatcher specified here, then we lazily
        #' construct a Hessian for said observation type using the following
        #' method:
        #' (1) first, we look for a `hess_loglik.type`, where `type` is the
        #' observation type.
        #' (2) if (1) fails, then we use a finite difference method given the
        #' log-likelihood contribution for the observation type.
        #' @param assumptions list of assumptions made by the model, default is
        #' c("iid""), which assumes iid observations (this assumption is always
        #' made for this class, which is why we can sum the log-likelihood
        #' contributions)
        #' @return A new `likelihood_contr_model` object
        initialize = function(obs_type,
                              logliks = NULL,
                              scores = NULL,
                              hess_logliks = NULL,
                              assumptions = c("iid")) {

            self$obs_type <- obs_type
            self$logliks <- logliks
            self$scores <- scores
            self$hess_logliks <- hess_logliks

            # iid assumption is always made
            self$assumptions <- unique(c("iid", assumptions))
        },

        #' @description
        #' Computes the log-likelihood of the model.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The total log-likelihood
        loglik = function(df, par, ...) {
            sapply(seq_len(nrow(df)), function(i) {
                self$get_loglik(self$obs_type(df[i, ]))(df[i, ], par, ...)
            }) |> sum()
        },

        #' @description
        #' Computes the score of the model.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The total score
        score = function(df, par, ...) {
            res <- sapply(seq_len(nrow(df)), function(i) {
                self$get_score(self$obs_type(df[i, ]))(df[i, ], par, ...)
            })
            if (ncol(df) == 1) {
                return(res |> sum())
            } else {
                return(res |> rowSums())
            }
        },

        #' @description
        #' Computes the Hessian of the log-likelihood.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The Hessian of the log-likelihood
        hess_loglik = function(df, par, ...) {
            res <- lapply(seq_len(nrow(df)), function(i) {
                self$get_hess_loglik(self$obs_type(df[i, ]))(df[i, ], par, ...)
            })
            return(Reduce("+", res))
        },

        #' @description
        #' Gets the loglikelihood contribution for an observation of `type`.
        #' 
        #' If the loglikelihood contribution for `type` does not have a
        #' corresponding method in the dispatcher, then we try to retrieve one
        #' from `loglik.type`, where `type` is the observation type. If this
        #' fails, then we throw an error.
        #' 
        #' @param type observation type
        #' @return The loglikelihood contribution for an observation
        get_loglik = function(type) {
            if (!(type %in% names(self$logliks))) {
                if (method_exists('loglik', type)) {
                    self$logliks[[type]] <- get(paste0('loglik.', type))
                } else {
                    stop(paste0("No `loglik` dispatcher for type: ", type))
                }
            }
            self$logliks[[type]]
        },

        #' @description
        #' Gets the score for an observation of `type`.
        #' 
        #' If the score for `type` does not have a corresponding method in the
        #' dispatcher, then we try to retrieve one from `score.type`, where
        #' `type` is the observation type. If this fails, then we use a
        #' finite difference method to compute the gradient of the
        #' log-likelihood function for the observation type.
        #' 
        #' @param type observation type
        #' @return The score for an observation
        #' @importFrom numDeriv grad
        #' @seealso \code{\link{numDeriv::grad}}
        #' @seealso \code{\link{likelihood_contr_model$get_loglik}}
        #' @export
        get_score = function(type) {
            if (!(type %in% names(self$scores))) {
                if (method_exists('score', type)) {
                    self$scores[[type]] <- get(paste0('score.', type))
                } else {
                    ll <- self$get_loglik(type)
                    s <- function(row, par) {
                        numDeriv::grad(func = function(par) { ll(row, par) },
                            x = par,
                            method.args = list(r = 6))
                    }
                    self$scores[[type]] <- s
                }
            }
            self$scores[[type]]
        },

        #' @description
        #' Gets the Hessian of the log-likelihood for an observation of `type`.
        #' 
        #' If the Hessian of the log-likelihood for `type` does not have a
        #' corresponding method in the dispatcher, then we try to retrieve one
        #' from `hess_loglik.type`, where `type` is the observation type. If
        #' this fails, then we use a finite difference method to compute the
        #' Hessian of the log-likelihood function for the observation type.
        #'
        #' @param type observation type
        #' @return The Hessian of the log-likelihood for an observation
        #' @importFrom numDeriv hessian
        #' @seealso \code{\link{numDeriv::hessian}}
        #' @seealso \code{\link{likelihood_contr_model$get_loglik}}
        #' @export
        get_hess_loglik = function(type) {
            if (!(type %in% names(self$hess_logliks))) {
                if (method_exists('hess_loglik', type)) {
                    self$hess_logliks[[type]] <-
                        get(paste0('hess_loglik.', type))
                } else {
                    ll <- self$get_loglik(type)
                    J <- function(row, par) {
                        numDeriv::hessian(func = function(par) { ll(row, par) },
                            x = par,
                            method.args = list(r = 6))
                    }
                    self$hess_logliks[[type]] <- J
                }
            }
            self$hess_logliks[[type]]
        }
    )
)

#' Log-likelihood method for likelihood_contr_model
#'
#' This method returns the log-likelihood function for a likelihood_contr_model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the log-likelihood given a data frame and 
#' parameters
#' @export
loglik.likelihood_contr_model <- function(model, ...) {
    function(df, par) model$loglik(df, par, ...)
}

#' Score method for likelihood_contr_model
#'
#' This method returns the score function for a likelihood_contr_model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the score given a data frame and parameters
#' @export
score.likelihood_contr_model <- function(model, ...) {
    function(df, par) model$score(df, par, ...)
}

#' Hessian of log-likelihood method for likelihood_contr_model
#'
#' This method returns the hessian of the log-likelihood function for a 
#' likelihood_contr_model
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the hessian of the log-likelihood given a
#' data frame and parameters
#' @export
hess_loglik.likelihood_contr_model <- function(model, ...) {
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
#' model <- likelihood_contr_model$new(
#'   logliks = list(
#'     normal = function(row, par) dnorm(row$x, par[1], par[2], log = TRUE)
#'   ),
#'   obs_type = function(row) "normal"
#' )
#' 
#' # compute FIM
#' fim(model, c(0, 1), df)
fim.likelihood_contr_model <- function(model, par, data, ...) {

    # compute scores
    #scores <- apply(data, 1, function(row) model$score(row, par, ...))

    # compute FIM
    #tcrossprod(scores) / nrow(data)


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


