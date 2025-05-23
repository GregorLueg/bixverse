% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_set_enrichment.R
\name{calc_gsea_traditional}
\alias{calc_gsea_traditional}
\title{Bixverse implementation of the simple fgsea algorithm}
\usage{
calc_gsea_traditional(
  stats,
  pathways,
  nperm = 2000L,
  gsea_params = params_gsea(),
  seed = 123L
)
}
\arguments{
\item{stats}{Named numeric vector. The gene level statistic.}

\item{pathways}{List. A named list with each element containing the genes for
this pathway.}

\item{nperm}{Integer. Number of permutation tests. Defaults to \code{2000L}.}

\item{gsea_params}{List. The GSEA parameters, see \code{\link[=params_gsea]{params_gsea()}}
wrapper function. This function generates a list containing:
\itemize{
\item min_size - Integer. Minimum size for the gene sets.
\item max_size - Integer. Maximum size for the gene sets.
\item gsea_param - Float. The GSEA parameter. Defaults to \code{1.0}.
}}

\item{seed}{Random seed for reproducibility.}
}
\value{
To be written.
}
\description{
Rust-based version of the traditional permutation-based GSEA algorithm.
}
