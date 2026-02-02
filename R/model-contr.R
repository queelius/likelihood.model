#' @title Likelihood_contr_model
#' @description
#' This class encapsulates all necessary parts of a likelihood model.
#' A likelihood_model should provide the following methods:
#' 
#' - `loglik`: computes the log-likelihood of the model
#' - `score`: computes the score of the model
#' - `hess_loglik`: computes the Hessian of the log-likelihood given
#'  a data frame and parameters
#' 
#' It provides methods for computing the log-likelihood, score, and
#' the Hessian of log-likelihood.
#' 
#' It also allows for different likelihood contributions depending
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

            if (is.null(obs_type) || !is.function(obs_type)) {
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
            dfs <- split(df, self$obs_type(df))
            res <- sapply(names(dfs), function(type) {
                self$get_score(type)(dfs[[type]], par, ...)
            })
            if (is.matrix(res)) rowSums(res) else sum(res)
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
            dfs <- split(df, self$obs_type(df))
            Reduce("+", lapply(names(dfs), function(type) {
                self$get_hess_loglik(type)(dfs[[type]], par, ...)
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
#' Returns a log-likelihood function that caches the data frame split and
#' resolved dispatchers. During optimization, `optim` calls the returned
#' function hundreds of times with the same `df` but varying `par` --
#' caching eliminates repeated `split()` and `obs_type()` overhead.
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the log-likelihood given a data frame and
#' parameters
#' @export
loglik.likelihood_contr_model <- function(model, ...) {
    cache <- .new_contr_cache()
    function(df, par, ...) {
        c <- .ensure_cache(cache, df, model, "loglik")
        sum(vapply(names(c$dfs), function(type) {
            c$fns[[type]](c$dfs[[type]], par, ...)
        }, numeric(1)))
    }
}

#' Score method for likelihood_contr_model
#'
#' Returns a score function with the same caching strategy as
#' [loglik.likelihood_contr_model()].
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the score given a data frame and parameters
#' @export
score.likelihood_contr_model <- function(model, ...) {
    cache <- .new_contr_cache()
    function(df, par, ...) {
        c <- .ensure_cache(cache, df, model, "score")
        res <- sapply(names(c$dfs), function(type) {
            c$fns[[type]](c$dfs[[type]], par, ...)
        })
        if (is.matrix(res)) rowSums(res) else sum(res)
    }
}

#' Hessian of log-likelihood method for likelihood_contr_model
#'
#' Returns a hessian function with the same caching strategy as
#' [loglik.likelihood_contr_model()].
#'
#' @param model The likelihood model
#' @param ... Additional arguments
#'
#' @return A function to compute the hessian of the log-likelihood given a
#' data frame and parameters
#' @export
hess_loglik.likelihood_contr_model <- function(model, ...) {
    cache <- .new_contr_cache()
    function(df, par, ...) {
        c <- .ensure_cache(cache, df, model, "hess_loglik")
        Reduce("+", lapply(names(c$dfs), function(type) {
            c$fns[[type]](c$dfs[[type]], par, ...)
        }))
    }
}

#' Create a new cache environment for likelihood_contr_model wrappers.
#' @return An environment with df_ref, dfs, and fns slots
#' @keywords internal
.new_contr_cache <- function() {
    env <- new.env(parent = emptyenv())
    env$df_ref <- NULL
    env$dfs <- NULL
    env$fns <- NULL
    env
}

#' Ensure the cache is populated for the given data frame.
#'
#' If \code{df} is the same object as the cached one (checked via
#' \code{identical()}), returns the existing cache. Otherwise splits the
#' data frame and resolves dispatcher functions for each observation type.
#'
#' @param cache Cache environment from \code{.new_contr_cache()}
#' @param df Data frame
#' @param model A \code{likelihood_contr_model}
#' @param method One of \code{"loglik"}, \code{"score"}, \code{"hess_loglik"}
#' @return The cache environment (invisibly updated), with \code{$dfs} and
#'   \code{$fns} populated
#' @keywords internal
.ensure_cache <- function(cache, df, model, method) {
    if (!identical(df, cache$df_ref)) {
        cache$df_ref <- df
        cache$dfs <- split(df, model$obs_type(df))
        getter <- switch(method,
            loglik = model$get_loglik,
            score = model$get_score,
            hess_loglik = model$get_hess_loglik
        )
        cache$fns <- lapply(names(cache$dfs), getter)
        names(cache$fns) <- names(cache$dfs)
    }
    cache
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
