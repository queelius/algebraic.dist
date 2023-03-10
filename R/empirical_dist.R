#' \code{empirical_dist} models the concept of an empirical distribution.
#' We coerce the data into a data farme (tibble), so that whether a vector,
#' matrix, or data frame is given, it works the same. We do not use a matrix
#' because each of the components may be of a different type.
#' @export
NULL

#' @param data data to construct empirical distribution from. if matrix or
#'             data frame, each row is a joint observation, if a vector, each
#'             element is an observation. whatever data is, it must be
#'             convertible to a tibble.
#' @export
empirical_dist <- function(data)
{
  df <- as.tibble(data)
  structure(df,class=unique(c("empirical_dist",class(df))))
}

nobs.empirical_dist <- function(x)
{
  nrow(x)
}

obs.empirical_dist <- function(x)
{
  x
}

#' Method for obtaining the pdf of a \code{empirical_dist} object.
#'
#' @param x The object to obtain the pdf of.
#' @note sort tibble lexographically and do a binary search to find upper
#'       and lower bound in \code{log(nobs(x))} time.
#' @export
pdf.empirical_dist <- function(x)
{
  function(t)
  {
    count <- 0L
    for (o in x)
      if (t == o) count <- count + 1L
    count / nrow(x)
  }
}

#' Method for obtaining the sampler for a \code{empirical_dist} object.
#'
#' @param x The object to obtain the sampler of.
#' @export
sampler.empirical_dist <- function(x,...)
{
    function(n=1) x[sample(1:nrow(x),size=n,replace=T)]
}

#' Method for obtaining the expectation of \code{f} with respect to a
#' \code{empirical_dist} object \code{x}.
#'
#' @param x The distribution object.
#' @param g The function to take the expectation of.
#' @param ... Additional arguments to pass (to function \code{g}).
#' @export
expectation.empirical_dist <- function(x,g,...)
{
  res <- g(x[1,],...)
  for (i in 2:nrow(x))
    res <- res + g(x[i,],...)
  res / nrow(x)
}

#' Method for obtaining the mean of \code{empirical_dist} object \code{x}.
#'
#' @param x The distribution object.
#' @export
mean.empirical_dist <- function(x)
{
    expectation(x,function(t) t,...)
}

#' Method for obtaining the variance of \code{univariate_dist} object \code{x}.
#'
#' @param x The distribution object.
#' @export
variance.empirical_dist <- function(x)
{
    mu <- mean(x,...)
    expectation(x,function(t) (t-mu)^2)
}
