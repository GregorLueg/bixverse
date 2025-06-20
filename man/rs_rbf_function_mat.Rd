% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extendr-wrappers.R
\name{rs_rbf_function_mat}
\alias{rs_rbf_function_mat}
\title{Apply a Radial Basis Function (to a matrix)}
\usage{
rs_rbf_function_mat(x, epsilon, rbf_type)
}
\arguments{
\item{x}{Numeric Matrix. The distances you wish to apply the Gaussian kernel
onto.}

\item{epsilon}{Float. Epsilon parameter for the RBF.}

\item{rbf_type}{String. Needs to be from \code{c("gaussian", "bump", "inverse_quadratic")}.}
}
\value{
The affinities after the Kernel was applied.
}
\description{
Applies a radial basis function (RBF) to a given distance
matrix. Has at the option to apply a Gaussian, Bump or Inverse Quadratic
RBF.
}
