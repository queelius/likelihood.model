#' @title Likelihood_contr_model
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
#' 
#' @importFrom R6 R6Class
#' @importFrom numDeriv grad hessian
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

            if (is.null(obs_type) | !is.function(obs_type)) {
                stop("obs_type must be a function")
            }
            self$obs_type <- obs_type
            self$logliks <- logliks
            self$scores <- scores
            self$hess_logliks <- hess_logliks

            # iid assumption is always made
            self$assumptions <- unique(c("iid", assumptions))
            class(self) <- c("likelihood_contr_model", "likelihood_model")
        },

        #' @description
        #' Computes the log-likelihood of the model.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The total log-likelihood
        loglik = function(df, par, ...) {
            private$validate(df, par)
            dfs <- split(df, self$obs_type(df))
            # Compute log-likelihood for each type and sum results
            sum(sapply(names(dfs), function(type) {
                self$get_loglik(type)(dfs[[type]], par, ...)
            }))
        },

        #' @description
        #' Computes the score of the model.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The total score
        score = function(df, par, ...) {
            private$validate(df, par)
            # Split data frame by observation type
            dfs <- split(df, self$obs_type(df))
            # Compute score for each type and sum results
            res <- sapply(names(dfs), function(type) {
                self$get_score(type)(dfs[[type]], par, ...)
            })
            res |> rowSums()
        },

        #' @description
        #' Computes the Hessian of the log-likelihood.
        #' 
        #' @param df dataframe for computation
        #' @param par parameters for computation
        #' @param ... additional arguments
        #' @return The Hessian of the log-likelihood
        hess_loglik = function(df, par, ...) {
            private$validate(df, par)
            # Find unique observation types
            types <- unique(self$obs_type(df))
            # Compute and sum Hessians for each type
            Reduce("+", lapply(types, function(type) {
                obs_df <- df[self$obs_type(df) == type, ]
                self$get_hess_loglik(type)(obs_df, par, ...)
            }))
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
                if (private$check_method('loglik', type)) {
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
                if (private$check_method('score', type)) {
                    self$scores[[type]] <- get(paste0('score.', type))
                } else {
                    ll <- self$get_loglik(type)
                    s <- function(df, par) {
                        numDeriv::grad(
                            func = function(par) {
                                ll(df, par)
                            },
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
                if (private$check_method('hess_loglik', type)) {
                    self$hess_logliks[[type]] <-
                        get(paste0('hess_loglik.', type))
                } else {
                    ll <- self$get_loglik(type)
                    J <- function(df, par) {
                        numDeriv::hessian(func = function(par) { ll(df, par) },
                            x = par,
                            method.args = list(r = 6))
                    }
                    self$hess_logliks[[type]] <- J
                }
            }
            self$hess_logliks[[type]]
        }
    ),

    private = list(
        validate = function(df, par) {
            if (is.null(df) || !is.data.frame(df) || is.null(par)) {
                stop("df and par must be provided")
            }
        },

        check_method = function(method, type) {
            exists(paste0(method, '.', type))
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

#' Retrieve the assumptions in the likelihood contributions model.
#'
#' @param model The likelihood contribution model
#' @param ... Additional arguments
#'
#' @return A list of assumptions
#' @export
assumptions.likelihood_contr_model <- function(model, ...) {
    model$assumptions
}


