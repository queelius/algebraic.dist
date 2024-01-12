#' Takes an expression `e` and a list `vars` and returns a
#' lazy `edist` (expression distribution object), that is a subclass
#' of `dist` that can be used in place of a `dist` object.
#' 
#' @param e the expression to evaluate against the arguments.
#' @param vars the list of distributions (with variable names)
#'             to evaluate the expression `e` against.
#' @examples
#' edist(expression(x + y),
#'       list("x" = normal(mu = 0, var = 1),
#'            "y" = exponential(rate = 1)))
#' @return An `edist` object.
#' @export
edist <- function(e, vars) {
  structure(list(e = e, vars = vars),
            class = c(paste0("edist_", e, "_", names(vars)), "edist", "dist"))
}

#' Function to determine whether an object `x` is an `edist` object.
#' @param x The object to test
#' @export
is_edist <- function(x) {
  inherits(x, "edist")
}

#' Method for obtaining the parameters of an `edist` object.
#' @param x The object to obtain the parameters of
#' @return A named vector of parameters
#' @export
#' @examples
#' edist(expression(x + y),
#'      list("x" = normal(mu = 0, var = 1),
#'          "y" = exponential(rate = 1)))
#' params(x)
#' #> x.mu x.var y.rate
#' #>     0     1      1
params.edist <- function(x) {
  c(sapply(x$vars, params), names(x$vars))
}

#' Method for obtaining the variance-covariance matrix (or scalar)
#' 
#' @param x The `edist` object to retrieve the variance-covariance matrix from
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `edist` object
#' @export
vcov.edist <- function(x) {
  cov(sampler(x)(1000))
}

#' Method for obtaining the mean of an `edist` object.
#' @param x The `edist` object to retrieve the mean from
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `edist` object
#' @export
mean.edist <- function(x) {
  colMeans(sampler(x)(1000))
}

#' Method for printing an `edist` object.
#' @param x The object to print
#' @param ... Additional arguments to pass (not used)
#' 
#' @export
#' @examples
#' edist(expression(x + y),
#'      list("x" = normal(mu = 0, var = 1),
#'           "y" = exponential(rate = 1)))
#' #> Distribution x + y
#' #>   x ~ normal(mu = 0, var = 1)
#' #>   y ~ exponential(rate = 1)
print.edist <- function(x, ...) {
  cat("Distribution", x$e, "\n")
    for (i in 1:length(x$vars)) {
        cat("  ", names(x$vars)[i], "~", x$vars[[i]], "\n")
    }
}

#' Method for obtaining the distribution objects of an `edist` object.
#' @param x The `edist` object to obtain the distribution objects of
#' @export
dists <- function(x) {
  if (is_edist(x)) {
    x$vars
  } else {
    stop("x is not an edist object")
  }
} 

#' Method for obtaining the sampler of an `edist` object.
#' @param x The `edist` object to obtain the sampler of.
#' @param ... Additional arguments to pass into each of the `sampler`
#'            function generators.
#' @return A function that takes a number of samples `n`, `...`
#'         which is passed into the expression `x$e` and returns
#'         the result of applying the expression `x$e` to the
#'         sampled values.
#' @export
sampler.edist <- function(x, ...) {

  # we use `sampler` to obtain the sampler
  # for each `dist` object in `x$vars`
  # so, here is what `x$vars` looks like:
  # `x$vars <- list("x" = dist1, "y" = dist2)`
  # and here is an example of what `x$e` looks like:
  # `x$e <- expression(x + y)`
  # so, we want to sample from `dist1` and `dist2`
  # and then evaluate `x + y` against the sampled values
  samplers <- lapply(x$vars, sampler, ...)
  cnames <- names(x$vars)

  # todo: accept a param and override the params in the e$vars
  # distributions. each distribution should have a nparam method
  # we can call (including `edist`!), so use it
  function(n = 1, ...) {
    # we sample from each `dist` object
    samples <- sapply(samplers, function(sampler) sampler(n, ...))
    colnames(samples) <- cnames

    # we evaluate the expression `x$e` using `eval` in environment
    # x$vars against the sampled values
    esamples <- list()
    for (i in 1:ncol(samples)) {
      es <- eval(x$e, envir = samples[[i]])
      esamples[[i]] <- es
    }
    esamples
  }
}