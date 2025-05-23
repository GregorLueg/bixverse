% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extendr-wrappers.R
\name{rs_calc_gsea_stats}
\alias{rs_calc_gsea_stats}
\title{Rust implementation of the fgsea::calcGseaStat() function}
\usage{
rs_calc_gsea_stats(stats, gs_idx, gsea_param, return_leading_edge)
}
\arguments{
\item{stats}{Numeric vector. The gene level statistic. Needs to
sorted in descending nature.}

\item{gs_idx}{Integer vector. The indices of the gene set genes.}

\item{gsea_param}{Float. The GSEA parameter. Usually defaults to 1.0.}

\item{return_leading_edge}{Boolean. Return the leading edge indices.}
}
\value{
List with the following elements
\itemize{
\item gene_stat Enrichment score for that gene set
\item leading_edge Indicies of the leading edge genes.
}
}
\description{
Rust implementation of the fgsea::calcGseaStat() function
}
