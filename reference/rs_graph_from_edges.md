# Build a Rust neighbour graph from an edge list

**\[experimental\]** Builds an undirected graph. Each edge should be
listed once; duplicates are kept as separate entries, as in igraph.
Weights are stored as `f32`.

## Usage

``` r
rs_graph_from_edges(edges, weights, n_nodes)
```

## Arguments

- edges:

  Integer vector. Flattened edge pairs (1-indexed!),
  `c(from_1, to_1, from_2, to_2, ...)`.

- weights:

  Numeric vector. One weight per pair.

- n_nodes:

  Integer. Number of nodes.

## Value

External pointer to the Rust graph.
