% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exponential.R
\name{density.exponential}
\alias{density.exponential}
\title{Method to obtain the pdf of an \code{exponential} object.}
\usage{
\method{density}{exponential}(x, ...)
}
\arguments{
\item{x}{The object to obtain the pdf of}

\item{...}{Additional arguments (not used)}
}
\value{
A function that computes the pdf of the exponential distribution
at a given point \code{t}. Also accepts a \code{rate} argument that
determines the failure rate of the exponential distribution (defaults
to the failure rate of object \code{x}) and a \code{log} argument that determines
whether to compute the log of the pdf.
}
\description{
Method to obtain the pdf of an \code{exponential} object.
}
