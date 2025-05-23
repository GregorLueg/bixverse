% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_wrappers.R
\name{params_community_detection}
\alias{params_community_detection}
\title{Wrapper function to generate community detection parameters}
\usage{
params_community_detection(
  max_nodes = 300L,
  min_nodes = 10L,
  min_seed_nodes = 2L,
  initial_res = 0.5
)
}
\arguments{
\item{max_nodes}{Integer. Maximum number of nodes in a given community.}

\item{min_nodes}{Integer. Minimum number of nodes in a given community.}

\item{min_seed_nodes}{Integer. Minimum number of seed nodes within a
community.}

\item{initial_res}{Float. Initial resolution parameter to start with.}
}
\value{
List with parameters for usage in subsequent function.
}
\description{
Wrapper function to generate community detection parameters
}
