library(tidyverse)
library(algebraic.mle)
devtools::load_all()
#library(likelihood.model)


#################################################################
# simulation data
####################################################################
n <- 20000
rates <- c(1.5, 1.25, 2.0)
clear.latent <- TRUE
set.seed(12325)

dgp <- function(rates, n, clear.latent = TRUE) {
    df <- data.frame(t1 = rexp(n, rates[1]),
                     t2 = rexp(n, rates[2]),
                     t3 = rexp(n, rates[3]))
    df$t <- apply(df, 1, min)
    df$k <- apply(df, 1, which.min)

    if (clear.latent) {
        for (i in 1:n) {
            for (p in 1:3) {
                if (p == df$k[i]) next
                df[i, p] <- NA
            }
        }
    }

    df
}

df <- dgp(rates, n, clear.latent)
head(df)

####################################################################
# score, loglike, hess loglik functions for observed
################################################
score.observed <- function(df, rates, ...) {
    s <- rep(-sum(df$t), length(rates))
    for (i in seq_len(nrow(df))) {
        s[df$k[i]] <- s[df$k[i]] + 1 / rates[df$k[i]]
    }
    s
}

loglik.observed <- function(df, rates, ...) {
    sum(log(rates[df$k])) - sum(rates) * sum(df$t)
}

hess_loglik.observed <- function(df, rates, ...) {
    p <- length(rates)
    H <- matrix(0, p, p)
    for (i in seq_len(nrow(df))) {
        j <- df$k[i]
        H[j, j] <- H[j, j] - 1 / rates[j]^2
    }
    H
}

####################################################################
# likelihood model
################################################
model.observed1 <- likelihood_contr_model$new(
    obs_type = function(df) { "observed" },
)

model.observed2 <- likelihood_contr_model$new(
    obs_type = function(df) { "what" },
    loglik = list(what = loglik.observed),
)

model.observed3 <- likelihood_contr_model$new(
    obs_type = function(df) { "observed" },
)



####################################################################
# MLE
################################################
# let's call use make optim control list so that the tolerance
# is extremely small, say 1e-30
# this is because we want to compare the loglike values
# and the default tolerance is 1e-8
# and we want to make sure that the loglike values are the same
# for all methods

control1 <- list(maxit = 100000, pgtol = 1e-50, factr = 1e-50)
control2 <- list(maxit = 100000, abstol = 1e-50, reltol = 1e-50)
control3 <- control1

mle.solver.observed1 <- likelihood.model::mle(model.observed1)
mle.solver.observed2 <- likelihood.model::mle(model.observed2)
mle.solver.observed3 <- likelihood.model::mle(model.observed3)
mle.observed1 <- mle.solver.observed1(df, par0 = rates, control = control1)
mle.observed2 <- mle.solver.observed2(df, par0 = rates, control = control2)
mle.observed3 <- mle.solver.observed3(df, par0 = rates, lower = 0,
    control = control3, method = "L-BFGS-B")

confint(mle.observed1)
confint(mle.observed2)
confint(mle.observed3)

lls <- list("1" = loglike(mle.observed1),
            "2" = loglike(mle.observed2),
            "3" = loglike(mle.observed3))
print(lls)

####################################################################
# other
################################################

R <- 200000
df2 <- dgp(rates, n=R, clear.latent = FALSE)
mle.observed1 <- mle.solver.observed1(df2, par0 = rates, control = control1)
#mle.observed2 <- mle.solver.observed2(df2, par0 = rates, control = control2)

fim.ob.observed1 <- algebraic.mle::fim(mle.observed1)
#fim.ob.observed2 <- algebraic.mle::fim(mle.observed2)
fim.exp.observed1 <- fim(model.observed1)(df2, rates)
#fim.exp.observed2 <- fim.likelihood_model(model.observed2)(df2, rates)

sqrt(diag(fim.exp.observed1))
#round(diag(fim.exp.observed2), digits=5)
sqrt(diag(fim.ob.observed1/R))
#round(diag(fim.ob.observed2/R), digits=5)

F1 <- algebraic.mle::fim(mle.observed1)/R
sigma <- solve(fim.exp.observed1)
sd <- sqrt(diag(sigma))


# plot sd of component 1 vs sample size
N <- 100
plot(x=15:N,sd[1]/sqrt(15:N),
    type="l",
    col="red",
    ylab="Standard error",
    xlab="Sample size",
    main="Standard error of the MLE")
# plot sd of component 2 vs sample size
lines(x=1:N,sd[2]/sqrt(1:N),col="blue")
# plot sd of component 3 vs sample size
lines(x=1:N,sd[3]/sqrt(1:N),col="green")



##############################################################
# let's sample MLE a bunch of times using MC simulation
############################################
R <- 200000

score.ob <- function(df, par, ...) {
    s <- rep(-sum(df$t), length(par))
    for (i in seq_len(nrow(df))) {
        s[df$k[i]] <- s[df$k[i]] + 1 / par[df$k[i]]
    }
    s
}

loglik.ob <- function(df, par, ...) {
    sum(log(par[df$k])) - sum(par) * sum(df$t)
}

hess_loglik.ob <- function(df, par, ...) {
    p <- length(par)
    H <- matrix(0, p, p)
    for (i in seq_len(nrow(df))) {
        j <- df$k[i]
        H[j, j] <- H[j, j] - 1 / par[j]^2
    }
    H
}

model.ob <- likelihood_contr_model$new(
    obs_type = function(df) { "ob" },
)

control <- list(maxit = 10000)
mle.solver <- mle(model.ob)
par <- rates

emp.data <- dgp(par, n = R, clear.latent = FALSE)
fim.exp <- fim(model.ob)(df, par)
sigma.exp <- solve(fim.exp)
sd.exp <- sqrt(diag(sigma.exp))

N <- seq(10, 200000, 10)
it <- 1
coverage <- c(0, 0, 0)
coverage_matrix <- matrix(0, nrow = length(N), ncol = 3)
for (n in N) {
    df <- dgp(par, n = n, clear.latent = FALSE)
    mle.ob <- mle.solver(df, par0 = par, control = control)
    
    CI <- confint(mle.ob)
    # check to make sure par is in the 95% CI
    if (par[1] >= CI[1, 1] && par[1] <= CI[1, 2]) {
        coverage[1] <- coverage[1] + 1
    }
    if (par[2] >= CI[2, 1] && par[2] <= CI[2, 2]) {
        coverage[2] <- coverage[2] + 1
    }
    if (par[3] >= CI[3, 1] && par[3] <= CI[3, 2]) {
        coverage[3] <- coverage[3] + 1
    }

    if (it %% 100 == 0) {
        cat("n =", n, "\n")
        cat("coverages:", coverage / it, "\n")
    }

    coverage_matrix[it, ] <- coverage / it

    fim.ob <- algebraic.mle::fim(mle.ob)
    sigma.ob <- solve(fim.ob)
    sd.ob <- sqrt(diag(sigma.ob))

    it <- it + 1
}

print(coverate_matrix)