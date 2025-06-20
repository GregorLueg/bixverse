% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_go_gse.R
\name{fgsea_simple_go_elim}
\alias{fgsea_simple_go_elim}
\title{Run GO enrichment with elimination with fgsea simple}
\usage{
fgsea_simple_go_elim(
  object,
  stats,
  elim_threshold = 0.05,
  nperm = 2000L,
  gsea_params = params_gsea(max_size = 2000L),
  seed = 42L,
  .debug = FALSE
)
}
\arguments{
\item{object}{The underlying class, see \code{\link[=gene_ontology_data]{gene_ontology_data()}}.}

\item{stats}{Named numeric vector. The gene level statistic.}

\item{elim_threshold}{Float. Threshold from which p-value onwards the
elimination on the ancestors shall be conducted.}

\item{nperm}{Integer. Number of permutation tests. Defaults to \code{2000L}}

\item{gsea_params}{List. The GSEA parameters, see \code{\link[=params_gsea]{params_gsea()}}
wrapper function. This function generates a list containing:
\itemize{
\item min_size - Integer. Minimum size for the gene sets.
\item max_size - Integer. Maximum size for the gene sets.
\item gsea_param - Float. The GSEA parameter. Defaults to \code{1.0}.
\item sample_size - Integer. Number of samples to iterate through for the
multi-level implementation of fgsea.
\item eps - Float. Boundary for calculating the p-value. Used for the multi-
level implementation of fgsea.
}}

\item{seed}{Random seed for reproducibility.}

\item{.debug}{Boolean. Shall information from the Rust function be displayed.
For debugging purposes. Warning: should you run this command over a large
list, you will have a large print output!}
}
\value{
data.table with enrichment results.
}
\description{
This method takes the gene_ontology_data and a vector of gene level
statistics to perform fgsea (simple) leveraging ontological information.
It starts at the lowest levels of the ontology and tests if there is
significant enrichment for any GO terms. If the threshold of the p-value is
below the elimination threshold, the genes from this term will be removed
from all its ancestors. The function then proceeds to the next level of the
ontology and repeats the process.
}
