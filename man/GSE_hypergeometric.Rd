% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_gse.R
\name{gse_hypergeometric}
\alias{gse_hypergeometric}
\title{Gene set enrichment (GSE) based on a hypergeometric test.}
\usage{
gse_hypergeometric(
  target_genes,
  gene_set_list,
  gene_universe = NULL,
  threshold = 0.05,
  minimum_overlap = 3L,
  .verbose = FALSE
)
}
\arguments{
\item{target_genes}{Character vector. GeneID(s) of the target genes.}

\item{gene_set_list}{Named list of character vectors. Names should represent
the gene sets, pathways, and the elements the genes within the respective
gene set.}

\item{gene_universe}{Optional character vector. If you would like to specify
specifically the gene universe. If set to NULL, the function will default to
all represented genes in the \code{gene_set_list}.}

\item{threshold}{Float between 0 and 1 to filter on the fdr. Default: 0.05.
If 1 everything is returned.}

\item{minimum_overlap}{Number of minimum overlap between the target genes
and the respective gene set.}

\item{.verbose}{Boolean. Controls verbosity of the function.}
}
\value{
data.table with enrichment results.
}
\description{
Takes a set of target genes, a list of gene sets and calculates a p-value
(hypergeometric test) and odds ratio (OR) against all the gene sets. Also
applies a multiple hypothesis correction (BH) to the p-values.
}
