% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal.R
\name{cdf.normal}
\alias{cdf.normal}
\title{Method for obtaining the cdf of an \code{normal} object.}
\usage{
\method{cdf}{normal}(x, ...)
}
\arguments{
\item{x}{The object to obtain the cdf of}

\item{...}{Additional arguments to pass (not used)}
}
\value{
A function that computes the cdf of the normal distribution.
It accepts as input a parameter vector \code{q}, a mean vector \code{mu},
a variance \code{var}, and a \code{log} argument
determining whether to compute the log of the cdf. By default,
\code{mu} and \code{var} are the mean and variance of object \code{x} and \code{log}
is \code{FALSE}. Finally, it accepts additional arguments \code{...} to
pass to the \code{pnorm} function.
}
\description{
Method for obtaining the cdf of an \code{normal} object.
}
