# Rust implementation of Leiden and Louvain clustering

**\[experimental\]** Modularity-based Leiden (or Louvain, with
`refine = FALSE`) on an undirected weighted graph. Mirrors
`igraph::cluster_leiden(objective_function = "modularity")`. The result
depends on the seed only, not on the number of threads.

## Usage

``` r
rs_leiden_clustering(from, to, weights, n_nodes, leiden_params)
```

## Arguments

- from:

  Integer vector. 1-based source node of each edge.

- to:

  Integer vector. 1-based target node of each edge. Same length as
  `from`.

- weights:

  Numeric vector. Edge weights, same length as `from`.

- n_nodes:

  Integer. Number of nodes in the graph.

- leiden_params:

  List. Named list with the optional fields (defaults in parentheses):

  - `resolution` numeric, modularity resolution (1.0)

  - `beta` numeric, randomness of the refinement merges (0.01)

  - `n_iterations` integer, full Leiden iterations (2)

  - `refine` logical, `TRUE` for Leiden, `FALSE` for Louvain (TRUE)

  - `seed` integer, random seed (42)

## Value

Integer vector with the 1-based community of each node; `1` is the
largest community.
