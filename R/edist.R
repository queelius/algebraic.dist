#' Takes an expression `e` and a list `vars` and returns a
#' lazy `edist` (expression distribution object), that is a subclass
#' of `dist` that can be used in place of a `dist` object.
#' 
#' @param e the expression to evaluate against the arguments.
#' @param vars the list of distributions (with variable names)
#'             to evaluate the expression `e` against.
#' @examples
#' edist(expression(x + y),
#'       list(x = normal(mu = 0, var = 1),
#'            y = exponential(rate = 1)))
#' @return An `edist` object.
#' @export
edist <- function(e, vars) {
  # retrieve the class of each of the distributions
  # to include that in the class of the edist object
  # we are creating
  classes <- sapply(vars, class)
  structure(list(e = e, vars = vars),
            class = c(paste0("edist_", e, "_", classes, "edist", "dist")))
}

#' Generic method for simplifying distributions.
#' 
#' @param x The distribution to simplify
#' @return The simplified distribution
#' @export
simplify <- function(x) UseMethod("simplify")


#' Default Method for simplifming a `dist` object.
#' @param x The `dist` object to simplify
#' @return The `dist` object
#' @export
simplify.dist <- function(x) {
  x
}

#' Method for simplifying an `edist` object.
#' 
#' This is a complicated function that walks the expression tree
#' and tries to simplify it. Since sometimes a simplification
#' made at some level of the expression tree can lead to a
#' simplification at a higher level of the expression tree,
#' we need to walk the expression tree from the bottom up.
#' 
#' Also, since some simplifications can lead to a change in
#' the class of the distribution, we need to be careful to
#' update the class of the distribution as we simplify it.
#' 
#' Finally, the simplifications we initially choose to do can
#' prevent us from doing other simplifications that may ultimately
#' be more beneficial. So, we need to try all valid simplifications,
#' creating new `edist` objects for each simplification, and then
#' choose the one that is the most simplified.
#' 
#' @param x The `edist` object to simplify
#' @return The simplified object
#' @export
#' @examples
#' e <- edist(expression(x + y),
#'            list(x = exponential(rate = 2),
#'                 y = exponential(rate = 1)))
#' print(e)
#' #> Distribution x + y
#' #>   x ~ exponential(rate = 2)
#' #>   y ~ exponential(rate = 1)
#' e.simp <- simplify(e)
#' print(e.simp)
#' #> Exponential distribution with rate 3
simplify.edist <- function(x) {

  

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
#' @param n The number of samples to take (default: 1000)
#' @param ... Additional arguments to pass (not used)
#' @return The variance-covariance matrix of the `edist` object
#' @export
vcov.edist <- function(x, n = 1000, ...) {
  samp <- sampler(x)(n)
  # if we have a matrix, we want the covariance matrix
  if (is.matrix(samp)) cov(samp)
  # otherwise, we just want the variance
  else var(samp)
}

#' Method for obtaining the mean of an `edist` object.
#' @param x The `edist` object to retrieve the mean from
#' @param n The number of samples to take (default: 1000)
#' @param ... Additional arguments to pass (not used)
#' @return The mean of the `edist` object
#' @export
mean.edist <- function(x, n = 1000, ...) {
  samp <- sampler(x)(n)
  # if we have a matrix, we want the column means
  if (is.matrix(samp)) colMeans(samp)
  # otherwise, we just want the mean
  else mean(samp)
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
  cat("Distribution", deparse(x$e), "\n")
  for (i in 1:length(x$vars)) {
    cat("  ", names(x$vars)[i], " ~ ")
    print(x$vars[[i]])
    cat("\n")
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
    esamples <- list()
    for (i in 1:nrow(samples)) {
      es <- eval(x$e, envir = as.list(samples[i,]))
      esamples[[i]] <- es
    }


    # convert esamples to vector or matrix, depending on
    # dimensionality of the first element in the list
    if (is.matrix(esamples[[1]])) {
      esamples <- do.call(rbind, esamples)
    } else {
      esamples <- unlist(esamples)
    }
    esamples
  }
}