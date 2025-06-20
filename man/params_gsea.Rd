% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_wrappers.R
\name{params_gsea}
\alias{params_gsea}
\title{Wrapper function to generate GSEA parameters}
\usage{
params_gsea(
  min_size = 5L,
  max_size = 500L,
  gsea_param = 1,
  sample_size = 121L,
  eps = 1e-50
)
}
\arguments{
\item{min_size}{Integer. Minimum number of genes per gene set.}

\item{max_size}{Integer. Maximum number of genes per gene set.}

\item{gsea_param}{Float. GSEA parameter. Defaults to \code{1.0}.}

\item{sample_size}{Integer. Number of samples to iterate through for the
multi-level implementation of fgsea.}

\item{eps}{Float. Boundary for calculating the p-value. Used for the multi-
level implementation of fgsea.}
}
\value{
List with parameters for usage in subsequent function.
}
\description{
Wrapper function to generate GSEA parameters
}
