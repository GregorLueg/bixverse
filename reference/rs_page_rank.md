# Rust version of calcaluting the personalised page rank

Rust version of calcaluting the personalised page rank

## Usage

``` r
rs_page_rank(node_names, from, to, weights, personalised, undirected)
```

## Arguments

- node_names:

  String vector. Name of the graph nodes.

- from:

  String vector. The names of the `from` edges from the edge list.

- to:

  String vector. The names of the `to` edges from the edge list.

- weights:

  Optional weight vector. If NULL, defaults to 1.0 as weight for all
  edges.

- personalised:

  Numerical vector. The reset values. They must sum to 1 and be of same
  length of `node_names`!

- undirected:

  Boolean. Is this an undirected graph.

## Value

The personalised page rank values.
