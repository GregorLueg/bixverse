% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes_ontology.R
\name{ontology}
\alias{ontology}
\title{Ontology class}
\usage{
ontology(parent_child_dt, .verbose = TRUE)
}
\arguments{
\item{parent_child_dt}{A data.table that contains the ontological information
in terms of parent child relationships. Need to contain the
\code{c("parent", "child")} columns.}

\item{.verbose}{Boolean. Controls the verbosity of the class}
}
\value{
Returns the class for subsequent usage.
}
\description{
This class is used to store any ontology and apply different methods to it.
Currently implement are semantic similarity calculations based on the
ontological information.
}
\section{Properties}{

\describe{
\item{edge_dt}{data.table. Contains the parent-child relationships.}
\item{information_content_list}{List. Contains the information content
of each individual term.}
\item{ancestor_list}{List. Contains the ancestors for each ontology term.}
\item{descendants_list}{List. Contains the descendants for each ontology
term.}
\item{semantic_similarities}{data.table. Contains the semantic similarities
if calculated.}
\item{params}{A (nested) list that will store all the parameters of the
applied function.}
\item{final_results}{Final results stored in the class.}
}
}

