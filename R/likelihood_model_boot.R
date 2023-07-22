#' @name likelihood_model_boot
#' @title Likelihood model bootstrapping
#' @description
#' Bootstrap a likelihood model to get standard errors and confidence intervals
#' for the parameters when the asymptotic distribution may not kick in.
#' 
#' When you do the bootstrap on the model, it returns a boot object
#' (see `boot`), but specialized for maximum likelihood estimation, and so
#' we implement the `mle` API for it and assign it the class `mle_boot`.
#' 
#' For many methods, it just calls the boot method, e.g., CIs can be computed
#' with `boot.ci`, but we also provide methods for implementnig other kinds
#' of methods appropriate to `mle` objects, like `vcov` (variance-covariance
#' matrix) and such.
NULL

#' Bootstrap a model. The statistics that are bootstrapped are the
#' parameters of the model, or other salient characteristics that
#' are important to the model. In a likelihood model, for instance,
#' the parameters are bootstrapped.
#' 
#' @param model A likelihood model
#' @param ... Additional arguments to `statistic`
#' @return A function that returns an empirical sampling distribution of an
#' MLE. 
#' @export
bootstrap <- function(model, ...) {
    UseMethod("bootstrap")
}

#' Bootstrap a likelihood model
#' 
#' @param model The likelihood model
#' @param algo The algorithm to use for the MLE solver. It should
#' look like the `stats::optim` function, i.e., take similar arguments
#' (if you don't do anything with them in your algorithm, just accept them
#' as `...` and don't do anything with them) and, especially, return a list with
#' the same named values (it may also return more). Defaults to `stats::optim`.
#' @param ... Additional arguments to pass into the likelihood model's
#' `loglik`, `score`, and `hess_loglik` construction methods
#' @return A function that returns an empirical sampling distribution of an
#' MLE.
#' @importFrom boot boot
#' @importFrom stats optim
#' @export
bootstrap.likelihood_model <- function(model, algo = stats::optim, ...) {

    ll <- loglik(model, ...)
    s <- score(model, ...)
    #H <- hess_loglik(model, ...)
    #mle_solver <- mle(model, ...)

    function(
        df,
        par,
        R = 999,
        method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
        lower = -Inf,
        upper = Inf,
        sim = c("ordinary", "parametric",
                 "balanced", "permutation", "antithetic"),
        control = list(), ...) {

        stopifnot(!is.null(par), is.data.frame(df))
        sim <- match.arg(sim)
        method <- match.arg(method)

        mle.boot <- function(df, indices, ...) {
            df.b <- df[indices, ]
            res.b <- optim(
                par = par,
                fn = function(par, ...) ll(df.b, par, ...),
                gr = function(par, ...) s(df.b, par, ...),
                ...,
                method = method,
                lower = lower,
                upper = upper,
                hessian = FALSE,
                control = control)
            # make sure par.boot converged?
            if (res.b$convergence != 0) {
                warning("MLE did not converge in bootstrap")
            }
            print(res.b$par)
            res.b$par
        }

        # we should have algebraic.mle::mle_boot also accept a true param
        # or MLE
        boot(data = df, statistic = mle.boot, R = R, sim = sim, ...)
    }
}



#' Bootstrap a likelihood model
#' 
#' In other words, we treat the data as an empirical distribution and
#' sample from it to get a new dataset, then we fit the model to that dataset
#' and return the MLE. We do this R times and return the R MLEs.
#' 
#' @param model The likelihood model
#' @param ... Additional arguments to pass into the likelihood model's
#' @param par The starting parameter values
#' @param R The number of bootstrap samples to take
#' @param nthreads The number of threads to use for parallelization
#' @param method_seq The sequence of optimization methods to use. By default,
#' we use SANN and then Nelder-Mead.
#' @return A function that returns an bootstrapped sampling distribution of an
#' MLE.
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom stats optim
#' @export
bootstrap.likelihood_model <- function(
    model,
    par = NULL,
    R = 999,
    algo = stats::optim,
    method_seq = c("SANN", "Nelder-Mead"),
    nthreads = NULL,
    ...) {

    ll <- loglik(model, ...)
    s <- score(model, ...)
    if (is.null(nthreads)) {
        nthreads <- detectCores() - 1L
    }

    function(df, par = par, R = R, method_seq = method_seq, ...) {

        stopifnot(!is.null(par), is.data.frame(df))

        n <- nrow(df)
        sol <- function(ind) {
            df.b <- df[ind, ]
            est <- par
            for (method in method_seq) {
                est <- optim(
                    par = est,
                    fn = function(par) ll(df.b, par),
                    gr = if (method == "SANN") {
                        NULL
                    } else {
                        function(par) s(df.b, par)
                    },
                    method = method,
                    control = list(fnscale = -1))
                est <- est$par
            }
            est$par
        }
        registerDoParallel(cores = nthreads)
        boots <- foreach(i = seq_len(R), .combine = rbind) %dopar% sol(
            ind = sample.int(
                n = n,
                size = n,
                replace = TRUE)
        )
        stopImplicitCluster()

        boots
    }
}

