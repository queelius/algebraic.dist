#' Concstructor for `dfr_dist` objects
#' 
#' We assume that the hazard rate is a function of time and any other
#' predictors. We also assume that intergate(rate(t), 0, Inf) = infinity
#' and that the support is (0, Inf).
#' 
#' @param rate A function that computes the hazard rate at time `t`.
#' @param par The parameters of the distribution. Defaults to `NULL`,
#'            which means that the parameters are unknown.
#' @param eps The epsilon update for numerical integration. Defaults to 0.01.
#' @return A `dfr_dist` object.
#' @export
dfr_dist <- function(rate, par = NULL, eps = 0.01) {
    structure(
        list(rate = rate,
             par = NULL,
             eps = eps),
    class = c("dfr_dist", "univariate_dist", "dist"))
}

#' Function for determining whether an object is a `dfr_dist` object.
#' 
#' @param x The object to test.
#' @return A logical value indicating whether `x` is a `dfr_dist` 
#' object.
#' @export
is_dfr_dist <- function(x) {
    inherits(x, "dfr_dist")
}

#' Method for obtaining the hazard function of a `dfr_dist` object.
#'
#' @param x The object to obtain the hazard function of.
#' @return A function that computes the hazard function of the distribution.
#' It accepts `t`, the time at which to compute the hazard function, and
#' `par`, the parameters of the distribution. If `par` is `NULL`, then
#' the parameters of the `dfr_dist` object `x` are used. It also
#' accepts a `log` argument that determines whether to compute the log of
#' the hazard function. Finally, it passes any additional arguments to the
#' `rate` function of the `dfr_dist` object `x`.
#' @export
hazard.dfr_dist <- function(x) {
    function(t, par = NULL, ...) {
        par <- params(x, par)
        x$rate(t, par, ...)
    }
}

#' Method for obtaining the quantile (inverse cdf) of an object.
#'
#' @param x The object to obtain the inverse cdf of.
#' @param ... Additional arguments to pass into `cdf` constructor.
#' @return A function that computes the quantile of the distribution.
#' It accepts `p`, the probability at which to compute the quantile, 
#' `par`, the parameters of the distribution, and `...`, any additional
#' arguments to pass into the constructed cdf.
#' @export
quantile.dfr_dist <- function(x, ...) {

    F <- cdf(x, ...)
    function(p, par = NULL, ...) {
        par <- params(x, par)
        # F(t) = p, so we want to find t such that F(t) = p.
        uniroot_args <- list(
            f = function(t) {
                F(t, par, ...) - p
            },
            interval = c(0, 1e3),
            extendInt = "upX")
        do.call(uniroot, uniroot_args)$root
    }
}

#' checks to see if any of the parameters are `NA` or `NULL`
#' and if so, replaces them with the default values
#' in `x` (`x$par`)
#'
#' We may think of this function as a "constructor" for the parameters
#' of a `dfr_dist` object. If we don't know some parameters, we can
#' pass `NA` or `NULL` for them, and this function will replace them
#' with the default values in `x$par`.
#'
#' @export
params.dfr_dist <- function(x, par = NULL) {
    if (is.null(par)) {
        stopifnot(!is.null(x$par))
        return(x$par)
    }

    # do component-wise replacement of `par` with `x$par`
    # whenever `par` is `NA` or `NULL`
    for (i in 1:length(x$par)) {
        if (is.null(par[[i]]) || is.na(par[[i]])) {
            par[[i]] <- x$par[[i]]
        }
    }
    par
}

#' Sampling function for `dfr_dist` objects.
#' 
#' Since S(t,par) = exp(-cum_hz(t,par)), we can sample from the
#' distribution by letting t = 0 (or some other positive number if
#' we want to condition on T > t_min), sampling from an exponential
#' distribution with `lambda = rate(t, par)`, and then rejecting
#' the sample if `runif(1) > S(t, par)`. If accepted, add that
#' observation to the sample, otherwise reject it, let `t = t + eps`
#' where `eps` is some small number, and repeat. We continue this
#' process until we have `n` observations for the sample.
#'
#' @param x The object to obtain the sampler of.
#' @param ... Additional arguments to pass into the survival function
#' @return A function that samples from the distribution. It accepts
#' `n`, the number of samples to take, `t` is the time at which to start
#' sampling, `par` are the parameters of the distribution, and `eps` is
#' the  update for numerical integration. Finally, we pass additional
#' arguments `...` into the hazard function.
#' @export
sampler.dfr_dist <- function(x, ...) {
    S <- surv(x, ...)
    function(n, t = 0, par = NULL, eps = x$eps, ...) {
        par <- params(x, par)
        replicate(n, {
            while (runif(1) > S(t, par)) {
                t <- t + eps
            }
            t
        })
    }
}

#' Method for obtaining the cdf of a `dfr_dist` object.
#'
#' @param x The object to obtain the cdf of.
#' @param ... Additional arguments to pass into the `cum_haz`
#' constructor.
#' @return A function that computes the cdf of the distribution.
#' It accepts `t`, the time at which to compute the cdf, `par`,
#' the parameters of the distribution, `log.p` argument that
#' determines whether to compute the log of the cdf, `lower.limit`,
#' whether to compute the lower limit (F(t)) or upper limit
#' (S(t) = 1-F(t)). Finally, it passes any additional arguments `...`
#' to the `rate` function of the `dfr_dist` object `x`.
#' @export
cdf.dfr_dist <- function(x, ...) {
    H <- cum_haz(x, ...)
    function (t, par = NULL, log.p, lower.limit = TRUE, ...) {
        par <- params(x, par)
        haz <- H(t, par, ...)
        if (lower.limit) {
            p <- 1 - exp(-haz)
            return(ifelse(log.p, log(p), p))
        }
        else {            
            return(ifelse(log.p, -haz, exp(-haz)))
        }
    }
}

#' Method for obtaining the pdf of a `dfr_dist` object.
#' 
#' @param x The object to obtain the pdf of.
#' @param ... Additional arguments to pass.
#' @return A function that computes the pdf of the distribution.
#' It accepts `y`, the value at which to compute the pdf, `t`, the time
#' at which to compute the pdf, `par` is the parameters of the
#' distribution, and `log` determines whether to compute the log of
#' the pdf. Finally, it passes any additional arguments `...` to
#' the `rate` function of the `dfr_dist` object `x`.
#' @export
pdf.dfr_dist <- function(x, ...) {
    H <- cum_haz(x, ...)
    function(t, par = NULL, log = FALSE, ...) {
        par <- params(x, par)
        if (log) {
            log(x$rate(t, par, ...)) - H(t, par, ...)
        } else {
            x$rate(t, par, ...) * exp(-H(t, par, ...))
        }
    }
}


#' Method for retrieving the support of an object `x`.
#'
#' @param x The object to obtain the support of.
#' @param ... Additional arguments to pass.
#' @return A support object for `x`, an interval (0,Inf).
#' @export
sup.dfr_dist <- function(x) {
    interval$new(0, Inf, FALSE, FALSE)
}

#' Method for obtaining the hazard function of a `dfr_dist` object.
#' @param x The object to obtain the hazard function of.
#' @param ... Additional arguments to pass into the `integrate` function.
#' @return A function that computes the hazard function of the distribution.
#' It accepts `t`, the time at which to compute the hazard function, and
#' `par`, the parameters of the distribution. If `par` is `NULL`, then the
#' parameters of the `dfr_dist` object `x` are used. Finally, it passes any
#' additional arguments `...` to the `rate` function of the `dfr_dist`
#' object `x`.
#' @importFrom stats integrate
#' @importFrom utils modifyList
#' @export
cum_haz.dfr_dist <- function(x, ...) {

    integrator_defaults <- list(
        lower = 0, subdivisions = 1000L, abs_tol = 1e-3)
    integrator <- modifyList(integrator_defaults, list(...))

    function(t, par = NULL, ...) {
        par <- params(x, par)
        # we call `integrate` with the options in `int_options` and
        # we integrate `rate` over `t` with respect to the parameters
        # `par` and any additional arguments `...`
        res <- do.call(integrate,
            modifyList(integrator, list(
                f = function(t) x$rate(t, par, ...))))
        if (res$message != "OK") {
            warning(res$message)
        }
        if (res$abs.error > abs_tol) {
            warning("Absolute error in cumulative hazard is greater than 1e-3")
        }
        res$value
    }
}


#' Method for obtaining the survival function of a `dfr_dist` object.
#' @param x The object to obtain the survival function of.
#' @param ... Additional arguments to pass into the `cum_haz`
#' constructor.
#' @return A function that computes the survival function of the
#' distribution.
#' It accepts `t`, the time at which to compute the survival, `par`,
#' the parameters of the distribution, `log.p` argument that
#' determines whether to compute the log of the survival, and
#' it passes any additional arguments into the `rate` function of
#' the `dfr_dist` object `x`.
#' @export
surv.dfr_dist <- function(x, ...) {
    H <- cum_haz(x, ...)
    function (t = 0, par = NULL, log.p, ...) {
        par <- params(x, par)
        haz <- -H(t, par, ...)
        ifelse(log.p, -haz, exp(-haz))
    }
}


#' Print method for `dfr_dist` objects.
#' @param x The `dfr_dist` object to print.
#' @param ... Additional arguments (not used)
#' @export
print.dfr_dist <- function(x, ...) {
  cat("Dynamic failure rate (DFR) distribution with failure rate:\n")
  print(x$rate)

  cat("It has a survival function given by:\n")
  cat("    S(t|rate) = exp(-H(t,...))\n")
  cat("where H(t,...) is the cumulative hazard function.\n")
}





#' Log-likelihood contribution of a `dfr_dist` object.
#' @param x The `dfr_dist` object to compute the log-likelihood contribution
#' for a vector of exact failure times `t`.
#' @importFrom likelihood.model loglik
loglik.dfr_dist <- function(x, ...) {
    function(t, par = NULL, ...) {
        par <- params(x, par)
        log(x$rate(t, par, ...))
    }
}