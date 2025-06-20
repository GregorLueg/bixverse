% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_preprocessing.R
\name{preprocess_bulk_coexp}
\alias{preprocess_bulk_coexp}
\title{Process the raw data}
\usage{
preprocess_bulk_coexp(
  object,
  hvg = NULL,
  mad_threshold = NULL,
  scaling = FALSE,
  scaling_type = c("normal", "robust"),
  .verbose = TRUE
)
}
\arguments{
\item{object}{The underlying class, see \code{\link[=bulk_coexp]{bulk_coexp()}}.}

\item{hvg}{Integer or float. If an integer, the top \code{hvg} genes will be
included; if float, the float has to be between 0 and 1, representing the
percentage of genes to include.}

\item{mad_threshold}{Float. Instead of of selecting number or proportion of
genes, you can also provide a mad_threshold.}

\item{scaling}{Boolean. Shall the data be scaled.}

\item{scaling_type}{String. You have the option to use normal scaling or
robust scaling.}

\item{.verbose}{Boolean. Controls the verbosity of the function.}
}
\value{
Returns the class with the \code{processed_data} data slot populated and
applied parameters added to the \code{params} slot.
}
\description{
Function to do general pre-processing on top of the \code{\link[=bulk_coexp]{bulk_coexp()}}.
Options to do scaling, HVG selection, etc.
}
