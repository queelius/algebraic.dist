% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic_dist.R
\name{sup}
\alias{sup}
\title{Generic method for retrieving the support of a (dist) object \code{x}.}
\usage{
sup(x)
}
\arguments{
\item{x}{The object to obtain the support of.}
}
\value{
A support object for \code{x}.
}
\description{
The returned value should have the following operations:
\itemize{
\item \code{min}: a vector, the minimum value of the support for each component.
\item \code{max}: a vector, the maximum value of the support for each component.
\item \code{call}: a predicate function, which returns TRUE if the value is in
the support, and FALSE otherwise.
\item \code{sample}: a function, which returns a sample from the support. Note that
the returned value is not guaranteed to be in the support of \code{x}. You may need
to call \code{call} to check.
}
}
