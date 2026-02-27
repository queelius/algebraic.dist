#' Convert an object to a probability distribution.
#'
#' Generic method for converting objects (such as fitted models) into
#' distribution objects from the \code{algebraic.dist} package.
#'
#' @param x The object to convert to a distribution.
#' @param ... Additional arguments to pass to methods.
#' @return A \code{dist} object.
#' @examples
#' # Identity for existing distributions
#' d <- normal(0, 1)
#' identical(as_dist(d), d)
#' @export
as_dist <- function(x, ...) {
    UseMethod("as_dist", x)
}

#' @rdname as_dist
#' @export
as_dist.dist <- function(x, ...) {
    x
}

#' Generic method for obtaining the hazard function of an object.
#'
#' @param x The object to obtain the hazard function of.
#' @param ... Additional arguments to pass.
#' @return A function computing the hazard rate at given points.
#' @examples
#' x <- exponential(2)
#' h <- hazard(x)
#' h(1)  # hazard rate at t = 1 (constant for exponential)
#' @export
hazard <- function(x, ...) {
    UseMethod("hazard", x)
}

#' Generic method for obtaining the cdf of an object.
#'
#' @param x The object to obtain the cdf of.
#' @param ... Additional arguments to pass.
#' @return A function computing the cumulative distribution function.
#' @examples
#' x <- normal(0, 1)
#' F <- cdf(x)
#' F(0)    # 0.5 (median of standard normal)
#' F(1.96) # approximately 0.975
#' @export
cdf <- function(x, ...) {
    UseMethod("cdf", x)
}

#' Generic method for obtaining the survival function of an object.
#'
#' @param x The object to obtain the survival function of.
#' @param ... Additional arguments to pass.
#' @return A function computing the survival function \eqn{S(t) = P(X > t)}.
#' @examples
#' x <- exponential(1)
#' S <- surv(x)
#' S(0)  # 1 (survival at time 0)
#' S(1)  # exp(-1), approximately 0.368
#' @export
surv <- function(x, ...) {
    UseMethod("surv", x)
}

#' Generic method for obtaining the quantile (inverse cdf) of an object.
#'
#' @param x The object to obtain the quantile of.
#' @param ... Additional arguments to pass.
#' @return A function computing the quantile (inverse CDF).
#' @examples
#' x <- normal(0, 1)
#' Q <- inv_cdf(x)
#' Q(0.5)   # 0 (median of standard normal)
#' Q(0.975) # approximately 1.96
#' @export
inv_cdf <- function(x, ...) {
    UseMethod("inv_cdf", x)
}

#' Generic method for sampling from distribution-like objects.
#'
#' It creates a sampler for the `x` object. It returns a function
#' that accepts a parameter `n` denoting the number of samples
#' to draw from the `x` object and also any additional parameters
#' `...` are passed to the generated function.
#'
#' @param x the `x` object to create a sampler for
#' @param ... additional arguments to pass
#' @return A function that takes \code{n} and returns \code{n} samples.
#' @examples
#' x <- normal(0, 1)
#' samp <- sampler(x)
#' set.seed(42)
#' samp(5)  # draw 5 samples from standard normal
#' @export
sampler <- function(x, ...) {
    UseMethod("sampler", x)
}

#' Generic method for obtaining the parameters of an object.
#'
#' @param x The object to obtain the parameters of.
#' @return A named vector (or list) of distribution parameters.
#' @examples
#' x <- normal(5, 2)
#' params(x)  # mu = 5, var = 2
#'
#' y <- exponential(3)
#' params(y)  # rate = 3
#' @export
params <- function(x) {
    UseMethod("params", x)
}

#' Generic method for obtaining the number of parameters of 
#' distribution-like object `x`.
#'
#' @param x the object to obtain the number of parameters for
#' @return Integer; the number of parameters.
#' @examples
#' d <- empirical_dist(matrix(rnorm(30), ncol = 3))
#' nparams(d)  # 0 (non-parametric)
#' @export
nparams <- function(x) {
    UseMethod("nparams", x)
}

#' Generic method for obtaining the expectation of `f` with respect to
#' `x`.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass into `g`.
#' @return The expected value of \code{g(x)}.
#' @examples
#' \donttest{
#' x <- exponential(1)
#' # E[X] for Exp(1) is 1
#' expectation(x, function(t) t)
#' }
#' @export
expectation <- function(x, g, ...) {
    UseMethod("expectation", x)
}

#' Generic method for obtaining the marginal distribution of a distribution
#' object `x` over components `indices`.
#' @param x The distribution object.
#' @param indices The indices of the marginal distribution to obtain.
#' @return A distribution object for the marginal over \code{indices}.
#' @examples
#' x <- mvn(c(0, 0), diag(2))
#' m <- marginal(x, 1)  # marginal over first component
#' mean(m)               # 0
#' @export
marginal <- function(x, indices) {
    UseMethod("marginal", x)
}


#' Generic method for obtaining the conditional distribution of a distribution
#' object `x` given condition `P`.
#' @param x The empirical distribution object.
#' @param P The predicate function to condition `x` on
#' @param ... additional arguments to pass into `P`
#' @return A distribution object for the conditional distribution.
#' @examples
#' \donttest{
#' d <- empirical_dist(1:100)
#' # condition on values greater than 50
#' d_gt50 <- conditional(d, function(x) x > 50)
#' mean(d_gt50)
#' }
#' @export
conditional <- function(x, P, ...) {
    UseMethod("conditional", x)
}

#' Generic method for applying a map `f` to distribution object `x`.
#' @param x The distribution object.
#' @param g The function to apply.
#' @param ... Additional arguments to pass into `g`.
#' @return A distribution representing the push-forward of \code{x} through \code{g}.
#' @examples
#' \donttest{
#' d <- empirical_dist(1:20)
#' d_sq <- rmap(d, function(x) x^2)
#' mean(d_sq)  # E[X^2] for uniform 1..20
#' }
#' @export
rmap <- function(x, g, ...) {
    UseMethod("rmap", x)
}

#' Generic method for retrieving the support of a (dist) object `x`.
#' 
#' The returned value should have the following operations:
#'  - `min`: a vector, the minimum value of the support for each component.
#'  - `max`: a vector, the maximum value of the support for each component.
#'  - `call`: a predicate function, which returns TRUE if the value is in
#'    the support, and FALSE otherwise.
#'  - `sample`: a function, which returns a sample from the support. Note that
#'    the returned value is not guaranteed to be in the support of `x`. You may need
#'    to call `call` to check.
#' @param x The object to obtain the support of.
#' @return A support object for `x`.
#' @examples
#' x <- normal(0, 1)
#' S <- sup(x)
#' infimum(S)   # -Inf
#' supremum(S)  # Inf
#'
#' y <- exponential(1)
#' S2 <- sup(y)
#' infimum(S2)  # 0
#' @export
sup <- function(x) {
    UseMethod("sup", x)
}

#' Retrieve the observations used to construct a distribution-like object. This is
#' useful for obtaining the data used to construct an empirical distribution, but
#' it is also useful for, say, retrieving the sample that was used by a fitted
#' object, like an maximum likelihood estimate.
#'
#' @param x the object to retrieve the observations from
#' @return The data (matrix or vector) used to construct \code{x}.
#' @examples
#' d <- empirical_dist(1:10)
#' obs(d)  # returns the vector 1:10
#' @export
obs <- function(x) {
    UseMethod("obs", x)
}
