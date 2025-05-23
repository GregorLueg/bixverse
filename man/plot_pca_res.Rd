% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_bulk.R
\name{plot_pca_res}
\alias{plot_pca_res}
\title{Plot the PCA data}
\usage{
plot_pca_res(
  object,
  cols_to_plot = c("contrast_info", "sample_source"),
  pcs_to_plot = c("PC_1", "PC_2"),
  ...
)
}
\arguments{
\item{object}{The underlying class, see \code{\link[=bulk_dge]{bulk_dge()}}.}

\item{cols_to_plot}{String vector. The columns within the meta-data to plot.
Defaults to \code{c('contrast_info', 'sample_source')}}

\item{pcs_to_plot}{String vector of length 2.
Will default to \code{c("PC_1", "PC_2")}.}

\item{...}{additional parameters}
}
\value{
A plot if the PCA information was found. \code{NULL} if no PCA was found.
}
\description{
Plot the PCA data
}
