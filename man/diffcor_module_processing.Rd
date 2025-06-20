% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_coexp_cor.R
\name{diffcor_module_processing}
\alias{diffcor_module_processing}
\title{Prepare differential correlation-based module detection}
\usage{
diffcor_module_processing(
  object,
  background_mat,
  cor_method = c("pearson", "spearman"),
  .verbose = TRUE
)
}
\arguments{
\item{object}{The class, see \code{\link[=bulk_coexp]{bulk_coexp()}}. Ideally, you
should run \code{\link[=preprocess_bulk_coexp]{preprocess_bulk_coexp()}} before applying this function.}

\item{background_mat}{Numerical matrix. The background data set.}

\item{cor_method}{String. Option of \code{c("pearson", "spearman")}.}

\item{.verbose}{Boolean. Controls verbosity of the function.}
}
\value{
The class with added data to the properties for subsequent usage.
}
\description{
This function will calculate the differential correlation between the stored
data set in the class and another background data set. To do so, it uses a
Fisher transformation of the correlation coefficients and calculates a Z
score based on the delta. The function will automatically subset into shared
features between the two data sets.
}
