% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal.R
\name{density.normal}
\alias{density.normal}
\title{Method for obtaining the pdf of an \code{normal} object.}
\usage{
\method{density}{normal}(x, ...)
}
\arguments{
\item{x}{The object to obtain the pdf of}

\item{...}{Additional arguments to pass (not used)}
}
\value{
A function that computes the pdf of the normal distribution.
It accepts as input a parameter vector \code{x}, a mean vector \code{mu},
a variance-covariance matrix \code{var}, and a \code{log} argument
determining whether to compute the log of the pdf. By default,
\code{mu} and \code{var} are the mean and variance of object \code{x}.
}
\description{
Method for obtaining the pdf of an \code{normal} object.
}
