% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_coexp_cor.R
\name{cor_module_tom}
\alias{cor_module_tom}
\title{Update the correlation matrix to a TOM}
\usage{
cor_module_tom(object, signed = TRUE, version = c("v2", "v1"), .verbose = TRUE)
}
\arguments{
\item{object}{The class, see \code{\link[=bulk_coexp]{bulk_coexp()}}. You need to have
applied \code{\link[=cor_module_processing]{cor_module_processing()}} before applying this function.}

\item{signed}{Boolean. Do you want to use the signed or unsigned version.
Defaults to \code{TRUE}.}

\item{version}{String. One of \code{c("v2", "v1")}. Defaults to \code{"v2"}.}

\item{.verbose}{Boolean. Controls verbosity of the function.}
}
\value{
The class with added data to the properties for subsequent usage.
}
\description{
This function will update the correlation matrix to a topological overlap
matrix. It defaults to \code{"v2"} and the signed version, please see
\code{\link[=calculate_tom]{calculate_tom()}} for details.
}
