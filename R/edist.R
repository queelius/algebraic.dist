#' Takes an expression `e` and a list `vars` and returns a
#' lazy `edist` (expression distribution object), that is a subclass
#' of `dist` that can be used in place of a `dist` object.
#' 
#' @param e the expression to evaluate against the arguments.
#' @param vars the list of distributions (with variable names)
#'             to evaluate the expression `e` against.
#' @export
edist <- function(e, vars) {
    structure(list(
        e = e,
        vars = vars),
        class = c("edist", "dist"))  # c(paste0("edist(", e, ")", "edist")))
}

#' Function to determine whether an object `x` is an `expr_dist` object.
#' @param x The object to test
#' @export
is_edist <- function(x) {
    inherits(x, "edist")
}

#' Method for obtaining the sampler of an `expr_dist` object.
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

