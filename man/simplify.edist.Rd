% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algebra.R
\name{simplify.edist}
\alias{simplify.edist}
\title{Method for simplifying an \code{edist} object.}
\usage{
\method{simplify}{edist}(x, ...)
}
\arguments{
\item{x}{The \code{edist} object to simplify}

\item{...}{Additional arguments to pass (not used)}
}
\value{
The simplified object
}
\description{
This is a complicated function that walks the expression tree
and tries to simplify it. Since sometimes a simplification
made at some level of the expression tree can lead to a
simplification at a higher level of the expression tree,
we need to walk the expression tree from the bottom up.
}
\details{
Also, since some simplifications can lead to a change in
the class of the distribution, we need to be careful to
update the class of the distribution as we simplify it.

Finally, the simplifications we initially choose to do can
prevent us from doing other simplifications that may ultimately
be more beneficial. So, we need to try all valid simplifications,
creating new \code{edist} objects for each simplification, and then
choose the one that is the most simplified.
}
