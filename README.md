Likelihood contribution model
================

# `likelihood.model`

`likelihood.model` is an R package that allows you to…

## Installation

You can install the development version of `likelihood.model` from
[GitHub](https://github.com/queelius/likelihood.model) with:

``` r
#if (!require(devtools)) {
#    install.packages("devtools")
#}
#devtools::install_github("queelius/likelihood.model")

library(algebraic.mle)
library(likelihood.model)
```

    ## 
    ## Attaching package: 'likelihood.model'

    ## The following objects are masked from 'package:algebraic.mle':
    ## 
    ##     fim, score

## Example: Exponential series system

Consider a series system with
![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m
"m") components, where the
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p
"p")th component in the
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i
"i")-th system has lifetime ![T\_{i
p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_%7Bi%20p%7D
"T_{i p}") with rate
![\\lambda\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_p
"\\lambda_p"), and ![T\_{i
j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_%7Bi%20j%7D
"T_{i j}") for all
![i,j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i%2Cj
"i,j") are independent. The system lifetime is the minimum of the
component lifetimes, i.e., ![T\_i = \\min\\{T\_{i 1}, \\ldots, T\_{i
m}\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_i%20%3D%20%5Cmin%5C%7BT_%7Bi%201%7D%2C%20%5Cldots%2C%20T_%7Bi%20m%7D%5C%7D
"T_i = \\min\\{T_{i 1}, \\ldots, T_{i m}\\}").

Let’s draw a sample for this model for
![m=3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m%3D3
"m=3") and ![\\lambda =
(1, 1.25, 1.5)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda%20%3D%20%281%2C%201.25%2C%201.5%29
"\\lambda = (1, 1.25, 1.5)") for a sample size of
![n=75](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D75
"n=75").

``` r
n <- 70
rates <- c(1.1, 1.25, 1.5)
set.seed(123)

df <- data.frame(t1 = rexp(n, rates[1]),
                 t2 = rexp(n, rates[2]),
                 t3 = rexp(n, rates[3]))
df$t <- apply(df, 1, min)

# map each observation to the corresponding index (component) that was minimum
df$k <- apply(df, 1, which.min)

head(df)
```

    ##           t1        t2         t3          t k
    ## 1 0.76677933 1.3139236 0.44257008 0.44257008 3
    ## 2 0.52419116 1.2959302 1.07292810 0.52419116 1
    ## 3 1.20823170 2.0289164 0.40455806 0.40455806 3
    ## 4 0.02870669 1.2172397 0.06243488 0.02870669 1
    ## 5 0.05110089 0.3040114 0.21642352 0.05110089 1
    ## 6 0.28772838 0.1908104 1.17325487 0.19081037 2

Column
![t\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_p
"t_p") is the failure time of the
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p
"p")-th component, column
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t
"t") is the column for the series system lifetime, and column
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k
"k") is the component that causes the system failure.

Of course, normally we don’t know the component lifetimes, but only the
series system lifetime (or a censored version of it). However,
sometimes, we may know more or less. We will consider three cases:

1.  We know the component lifetimes exactly. (Columns ![t\_1, \\ldots,
    t\_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t_1%2C%20%5Cldots%2C%20t_m
    "t_1, \\ldots, t_m") are observed.)
2.  We know the system lifetime and the component cause of failure.
    (Columns
    ![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t
    "t") and
    ![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k
    "k") are observed.)
3.  We know the system lifetime, and some subset of the components that
    contains the component that caused the system failure. (New columns
    will be introduced to indicate which components are observed.)

### Case 1: Complete knowledge of component lifetimes

Suppose we have a series system, but we have complete knowledge of the
component lifetimes. The component lifetimes are exponentially
distributed with rate
![\\lambda\_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_k
"\\lambda_k") for component
![k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k
"k").

The system lifetime is the minimum of the component lifetimes, but we
actually know the component lifetimes somehow, so we just estimate each
![\\lambda\_k](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_k
"\\lambda_k") separately
(![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m
"m") independent MLE problems), which we know to have the MLE
![\\hat{\\lambda}\_p
= 1/{\\bar{t}\_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Clambda%7D_p%20%3D%201%2F%7B%5Cbar%7Bt%7D_p%7D
"\\hat{\\lambda}_p = 1/{\\bar{t}_p}") where
![\\bar{t}\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbar%7Bt%7D_p
"\\bar{t}_p") is the sample mean of the component lifetimes for
component
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p
"p").

``` r
(lambda.hat <- c(1/mean(df$t1), 1/mean(df$t2), 1/mean(df$t3)))
```

    ## [1] 1.070043 1.250383 1.533236

#### Likelihood model

Let’s derive the full likelihood model for this problem anyway. The
likelihood function is given by so we use those instead. The likelihood
function is given by ![L(\\lambda\_1, \\dots, \\lambda\_p) =
\\prod\_{i=1}^n \\prod\_{j=1}^m \\lambda\_j \\exp(-\\lambda\_j t\_{i
j})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;L%28%5Clambda_1%2C%20%5Cdots%2C%20%5Clambda_p%29%20%3D%20%5Cprod_%7Bi%3D1%7D%5En%20%5Cprod_%7Bj%3D1%7D%5Em%20%5Clambda_j%20%5Cexp%28-%5Clambda_j%20t_%7Bi%20j%7D%29
"L(\\lambda_1, \\dots, \\lambda_p) = \\prod_{i=1}^n \\prod_{j=1}^m \\lambda_j \\exp(-\\lambda_j t_{i j})"),
whose log-likelihood function is given by ![\\ell (\\lambda\_1, \\dots,
\\lambda\_p) = \\sum\_{i=1}^n \\sum\_{j=1}^m \\log \\lambda\_p -
\\lambda\_p t\_{i
j}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell%20%28%5Clambda_1%2C%20%5Cdots%2C%20%5Clambda_p%29%20%3D%20%5Csum_%7Bi%3D1%7D%5En%20%5Csum_%7Bj%3D1%7D%5Em%20%5Clog%20%5Clambda_p%20-%20%5Clambda_p%20t_%7Bi%20j%7D
"\\ell (\\lambda_1, \\dots, \\lambda_p) = \\sum_{i=1}^n \\sum_{j=1}^m \\log \\lambda_p - \\lambda_p t_{i j}"),
whose score function is given by ![s(\\lambda\_1, \\dots, \\lambda\_p) =
(\\partial \\ell / \\partial \\lambda\_1 \\cdots \\partial \\ell /
\\partial
\\lambda\_m)^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%28%5Clambda_1%2C%20%5Cdots%2C%20%5Clambda_p%29%20%3D%20%28%5Cpartial%20%5Cell%20%2F%20%5Cpartial%20%5Clambda_1%20%5Ccdots%20%20%5Cpartial%20%5Cell%20%2F%20%5Cpartial%20%5Clambda_m%29%5ET
"s(\\lambda_1, \\dots, \\lambda_p) = (\\partial \\ell / \\partial \\lambda_1 \\cdots  \\partial \\ell / \\partial \\lambda_m)^T")
where ![\\partial \\ell / \\partial \\lambda\_p = n/\\lambda\_p -
\\sum\_{i=1}^n t\_{i
p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpartial%20%5Cell%20%2F%20%5Cpartial%20%20%5Clambda_p%20%3D%20n%2F%5Clambda_p%20-%20%5Csum_%7Bi%3D1%7D%5En%20t_%7Bi%20p%7D
"\\partial \\ell / \\partial  \\lambda_p = n/\\lambda_p - \\sum_{i=1}^n t_{i p}"),
and whose Hessian matrix is a diagonal matrix whose
![(p,p)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28p%2Cp%29
"(p,p)")-th element is given by ![\\partial^2 \\ell / \\partial
\\lambda\_p^2 =
-\\frac{n}{\\lambda\_p^2}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cpartial%5E2%20%5Cell%20%2F%20%5Cpartial%20%5Clambda_p%5E2%20%3D%20-%5Cfrac%7Bn%7D%7B%5Clambda_p%5E2%7D
"\\partial^2 \\ell / \\partial \\lambda_p^2 = -\\frac{n}{\\lambda_p^2}").

To solve the MLE, we just solve the score equations ![s(\\lambda\_1,
\\dots, \\lambda\_p)
= 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;s%28%5Clambda_1%2C%20%5Cdots%2C%20%5Clambda_p%29%20%3D%200
"s(\\lambda_1, \\dots, \\lambda_p) = 0"), which has the unique solution
![\\hat{\\lambda}\_p
= 1/{\\bar{t}\_p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7B%5Clambda%7D_p%20%3D%201%2F%7B%5Cbar%7Bt%7D_p%7D
"\\hat{\\lambda}_p = 1/{\\bar{t}_p}").

So, the sampling distribution of
![\\hat\\lambda\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Clambda_p
"\\hat\\lambda_p") is given by   
![&#10;\\hat\\lambda\_p \\sim \\text{Gamma}(n, 1 /
\\lambda\_p),&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Chat%5Clambda_p%20%5Csim%20%5Ctext%7BGamma%7D%28n%2C%201%20%20%2F%20%5Clambda_p%29%2C%0A
"
\\hat\\lambda_p \\sim \\text{Gamma}(n, 1  / \\lambda_p),
")  
and asymptotically, by the CLT, it is given by   
![&#10;\\hat\\lambda\_p \\sim \\text{Normal}(\\lambda\_p, \\lambda\_p^2
/
n),&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Chat%5Clambda_p%20%5Csim%20%5Ctext%7BNormal%7D%28%5Clambda_p%2C%20%5Clambda_p%5E2%20%2F%20n%29%2C%0A
"
\\hat\\lambda_p \\sim \\text{Normal}(\\lambda_p, \\lambda_p^2 / n),
")  
which we may in either case replace
![\\lambda\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_p
"\\lambda_p") with
![\\hat\\lambda\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Clambda_p
"\\hat\\lambda_p") to get respectively the approximate sampling and
approximate asymptotic sampling distribution of
![\\hat\\lambda\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Clambda_p
"\\hat\\lambda_p").

``` r
model.exp.complete <- likelihood_model$new(
    logliks = list(
        observed = function(row, par) {
            dexp(row$t1, par[1], log = TRUE) +
            dexp(row$t2, par[2], log = TRUE) +
            dexp(row$t3, par[3], log = TRUE)
        }
    ),
    scores = list(
        observed = function(row, par) {
            c(1/par[1] - row$t1,
              1/par[2] - row$t2,
              1/par[3] - row$t3)
        }
    ),
    hess_logliks = list(
        observed = function(row, par) {
            -matrix(c(1/par[1]^2, 0         , 0,
                      0         , 1/par[2]^2, 0,
                      0         , 0         , 1/par[3]^2),
                nrow=length(par))
        }
    ),
    obs_type = function(row) {
        "observed"
    }
)
ll.exp.complete <- loglik(model.exp.complete)
s.exp.complete <- score(model.exp.complete)
H.exp.complete <- hess_loglik(model.exp.complete)
```

The MLE is just:

``` r
res.exp.complete <- optim(rates, fn = ll.exp.complete, gr = s.exp.complete,
    df = df, hessian = TRUE, control = list(fnscale = -1))
print(res.exp.complete$par)
```

    ## [1] 1.070296 1.250357 1.533323

``` r
print(lambda.hat)
```

    ## [1] 1.070043 1.250383 1.533236

``` r
ll.exp.complete(df, res.exp.complete$par)
```

    ## [1] -159.703

``` r
print(res.exp.complete$value)
```

    ## [1] -159.703

``` r
print(res.exp.complete$hessian)
```

    ##           [,1]      [,2]      [,3]
    ## [1,] -61.10696   0.00000   0.00000
    ## [2,]   0.00000 -44.77448   0.00000
    ## [3,]   0.00000   0.00000 -29.77357

``` r
print(H.exp.complete(df, res.exp.complete$par))
```

    ##          [,1]      [,2]      [,3]
    ## [1,] -61.1069   0.00000   0.00000
    ## [2,]   0.0000 -44.77445   0.00000
    ## [3,]   0.0000   0.00000 -29.77356

We see that the two MLEs, the analytical one and the numerical one, are
the same.

The sampling distribution of the MLE is given by:

``` r
mle.samp <- mle_numerical(res.exp.complete)
summary(mle.samp)
```

    ## Maximum likelihood estimator of type mle_numerical is normally distributed.
    ## The estimates of the parameters are given by:
    ## [1] 1.070296 1.250357 1.533323
    ## The standard error is  0.1279248 0.1494461 0.1832671 .
    ## The asymptotic 95% confidence interval of the parameters are given by:
    ##             2.5%    97.5%
    ## param1 0.8598784 1.280713
    ## param2 1.0045396 1.496174
    ## param3 1.2318754 1.834771
    ## The MSE of the estimator is  0.07228573 .
    ## The log-likelihood is  -159.703 .
    ## The AIC is  325.406 .

### Case 2: Known component cause, but masked component lifetimes

In this case, we assume we can only observe the system lifetime and the
component that caused the system failure.

We might be tempted to think that we can just consider each observation,
where we know the system lifetime and the componente cause, and then
just estimate the component lifetime for that component based on its
failure time. We have less information about the parameters, because
each observation only gives us information about one component, but
that’s okay.

However, this is incorrect, because when we condition on the component
that caused the system failure, we are not observing a random sample of
that component’s lifetime, but a conditional sample. For instance, if
![k
= 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;k%20%3D%201
"k = 1"), then we are observing ![T\_{i 1} | T\_{i 1} \< T\_{i 2}
\\text{ and } T\_{i 1} \<
T\_{i 3}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_%7Bi%201%7D%20%7C%20T_%7Bi%201%7D%20%3C%20T_%7Bi%202%7D%20%5Ctext%7B%20and%20%7D%20T_%7Bi%201%7D%20%3C%20T_%7Bi%203%7D
"T_{i 1} | T_{i 1} \< T_{i 2} \\text{ and } T_{i 1} \< T_{i 3}").

What is this distribution? We can derive it as:

  
![&#10; f\_{T\_i|K\_i}(t\_i | k\_i) = f\_{T\_i,K\_i}(t\_i,k\_i) /
f\_{K\_i}(k\_i) =&#10; f\_{K\_i|T\_i}(k\_i|t\_i) f\_{T\_i}(t\_i) /
f\_{K\_i}(k\_i).&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%20%20%20f_%7BT_i%7CK_i%7D%28t_i%20%7C%20k_i%29%20%3D%20f_%7BT_i%2CK_i%7D%28t_i%2Ck_i%29%20%2F%20f_%7BK_i%7D%28k_i%29%20%3D%0A%20%20%20%20%20%20%20%20f_%7BK_i%7CT_i%7D%28k_i%7Ct_i%29%20f_%7BT_i%7D%28t_i%29%20%2F%20f_%7BK_i%7D%28k_i%29.%0A
"
    f_{T_i|K_i}(t_i | k_i) = f_{T_i,K_i}(t_i,k_i) / f_{K_i}(k_i) =
        f_{K_i|T_i}(k_i|t_i) f_{T_i}(t_i) / f_{K_i}(k_i).
")  
By the memoryless property of the exponential distribution, we have   
![&#10; f\_{K\_i|T\_i}(k\_i|t\_i) =
f\_{K\_i}(k\_i),&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%20%20%20f_%7BK_i%7CT_i%7D%28k_i%7Ct_i%29%20%3D%20f_%7BK_i%7D%28k_i%29%2C%0A
"
    f_{K_i|T_i}(k_i|t_i) = f_{K_i}(k_i),
")  
since the failure rates are constant (and thus independent of the time),
and thus the probability that a componet is the cause of a system
failure is independent of the time of the system failure.

That means that   
![&#10; f\_{T\_i|K\_i}(t\_i | k\_i) =
f\_{T\_i}(t\_i),&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%20%20%20f_%7BT_i%7CK_i%7D%28t_i%20%7C%20k_i%29%20%3D%20f_%7BT_i%7D%28t_i%29%2C%0A
"
    f_{T_i|K_i}(t_i | k_i) = f_{T_i}(t_i),
")  
which is the density of the series system. Since the minimum of
exponentially distributed random variables is exponentially distributed
with a falure rate equal to the sum of the failure rates of the
components, this estimator is an estimator of the sum of the failure
rates (or the failure rate of the series system). Let’s do this
computation:

``` r
model.exp.series <- likelihood_model$new(
    logliks = list(
        observed = function(row, rate) {
            dexp(row$t, rate[row$k], log = TRUE)
        }
    ),
    obs_type = function(row) {
        "observed"
    }
)
ll.exp.series <- loglik(model.exp.series)
res.exp.series <- optim(rates, fn = ll.exp.series, df = df, hessian = TRUE,
    control = list(fnscale = -1))

print(res.exp.series$par)
```

    ## [1] 2.691628 4.273178 3.720431

Each of these are more or less reasonable estimates of the true failure
rate of the system, which is ![(1.1 + 1.25 + 1.5)
= 3.85](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%281.1%20%2B%201.25%20%2B%201.5%29%20%3D%203.85
"(1.1 + 1.25 + 1.5) = 3.85").

However, we can combine each of these estimates into a weighted mean to
get a better estimate. The weights are the number of observations for
each component. We actually went ahead and analytically computed the MLE
for this, also:

``` r
# print some summary statistics
# first, show how many rows have k = j
N <- rep(0, 3)
for (j in 1:3) {
    N[j] <- sum(df$k == j)
    cat("counts(", j, ") = ", N[j], "\n")
}
```

    ## counts( 1 ) =  17 
    ## counts( 2 ) =  24 
    ## counts( 3 ) =  29

``` r
t.sum <- sum(df$t)
cat("total system lifetimes = ", t.sum, "\n")
```

    ## total system lifetimes =  19.72759

``` r
# sum of system lifetimes for k = j 
t.sum.k <- rep(0, 3)
for (j in 1:3) {
    t.sum.k[j] <- sum(df$t[df$k == j])
    cat("sum of system lifetimes for k = ", j, " = ", t.sum.k[j], "\n")
}
```

    ## sum of system lifetimes for k =  1  =  6.31582 
    ## sum of system lifetimes for k =  2  =  5.616639 
    ## sum of system lifetimes for k =  3  =  7.795127

``` r
# MLE (analytical solution)
# lambda_j = N_j / t_j
(lambda.hat2 <- N / t.sum.k)
```

    ## [1] 2.691654 4.273018 3.720273

``` r
res.exp.series$par
```

    ## [1] 2.691628 4.273178 3.720431

``` r
res.exp.series$par[1] * N[1] / nrow(df) +
    res.exp.series$par[2] * N[2] / nrow(df) +
    res.exp.series$par[3] * N[3] / nrow(df)
```

    ## [1] 3.660092

``` r
print(sum(rates))
```

    ## [1] 3.85

We see that this is a pretty good estimator of the true failure rate of
the system, but it tells us nothing about the failure rate of the
components.

#### Joint distribution of ![T\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_i "T_i") and ![K\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_i "K_i")

Instead of maximizing the conditional likelihood based on ![T\_i |
K\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_i%20%7C%20K_i
"T_i | K_i"), we should be maximizing the joint likelihood of
![T\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T_i
"T_i") and
![K\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_i
"K_i"):

``` r
model.exp.known <- likelihood_model$new(
    logliks = list(
        observed = function(row, rate) {
            log(rate[row$k]) - sum(rate) * row$t
        }
    ),
    scores = list(
        observed = function(row, rate) {
            v <- rep(-row$t, length(rate))
            v[row$k] <- v[row$k] + 1 / rate[row$k]
            v
        }
    ),
    hess_logliks = list(
        observed = function(row, rate) {
            p <- length(rate)
            H <- matrix(0, p, p)
            H[row$k, row$k] <- -1 / rate[row$k]^2
        }
    ),
    obs_type = function(row) {
        "observed"
    }
)
ll.exp.known <- loglik(model.exp.known)
s.exp.known <- score(model.exp.known)
H.exp.known <- hess_loglik(model.exp.known)
```

``` r
rate.known <- optim(rates, fn = ll.exp.known, df = df,
    hessian = TRUE, control = list(fnscale = -1))

summary(mle_numerical(rate.known))
```

    ## Maximum likelihood estimator of type mle_numerical is normally distributed.
    ## The estimates of the parameters are given by:
    ## [1] 0.8619292 1.2164711 1.4701796
    ## The standard error is  0.2090483 0.248311 0.2730054 .
    ## The asymptotic 95% confidence interval of the parameters are given by:
    ##             2.5%    97.5%
    ## param1 0.5180754 1.205783
    ## param2 0.8080359 1.624906
    ## param3 1.0211257 1.919233
    ## The MSE of the estimator is  0.1798914 .
    ## The log-likelihood is  -56.65176 .
    ## The AIC is  119.3035 .
