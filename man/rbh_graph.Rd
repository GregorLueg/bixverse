% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes_graphs.R
\name{rbh_graph}
\alias{rbh_graph}
\title{Reciprocal best hit graph}
\usage{
rbh_graph(module_results, dataset_col, module_col, value_col)
}
\arguments{
\item{module_results}{data.table with the all of the gene modules for which
you wish to generate the RBH graph.}

\item{dataset_col}{The column (name) which indicates from which data set/
method the gene module was derived.}

\item{module_col}{The column (name) which stores the names of the modules.}

\item{value_col}{The column (name) which stores the genes that are part of
the modules.}
}
\value{
Returns the \code{rbh_graph} class for further operations.
}
\description{
This class can be used to generate reciprocal best hit graphs between
gene modules from different origins.

Internal properties:
\itemize{
\item module_data: Nested list of the modules generated by different methods.
\item rbh_graph: The reciprocal best hit igraph.
\item rbh_edge_df: Edge data.table with the information from the reciprocal
best hits.
\item final_results: This is where the final communities will be stored.
\item params: A list with the params. This will be populated during subsequent
function calls.
}
}
