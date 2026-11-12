# Flatten a Rust neighbour graph into an edge list

**\[experimental\]** Each undirected edge comes out once, ready for
[`igraph::add_edges()`](https://r.igraph.org/reference/add_edges.html).

## Usage

``` r
rs_graph_to_edges(ptr)
```

## Arguments

- ptr:

  External pointer to the Rust graph.

## Value

A list with

- edges - Integer vector with the flattened edge pairs (1-indexed!),
  `c(from_1, to_1, from_2, to_2, ...)`.

- weights - Numeric vector. One weight per pair.
