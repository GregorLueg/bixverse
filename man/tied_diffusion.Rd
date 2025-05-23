% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_graphs.R
\name{tied_diffusion}
\alias{tied_diffusion}
\title{Diffuse seed genes in a tied manner over a network}
\usage{
tied_diffusion(
  object,
  diffusion_vector_1,
  diffusion_vector_2,
  summarisation = c("max", "mean", "harmonic_sum"),
  score_aggregation = c("min", "max", "mean"),
  .verbose = FALSE
)
}
\arguments{
\item{object}{The underlying class \code{\link[=network_diffusions]{network_diffusions()}}.}

\item{diffusion_vector_1}{The first named vector with values to use for the
reset parameter in the personalised page-rank diffusion. Names should
represent node names of the graph.}

\item{diffusion_vector_2}{The second named vector with values to use for the
reset parameter in the personalised page-rank diffusion. Names should
represent node names of the graph.}

\item{summarisation}{If there are duplicated names in the \code{diffusion_vector}
how to summarise
these.}

\item{score_aggregation}{How to summarise the tied scores.}

\item{.verbose}{Controls verbosity of the function.}
}
\value{
The class with added diffusion score based on a two sets of seed
genes. Additionally, the seed genes are stored in the class.
}
\description{
This function takes two sets of diffusion vector and leverages tied diffusion
to identify an intersection of influential nodes.
}
