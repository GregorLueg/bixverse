# Summary of a Rust neighbour graph

**\[experimental\]** Reads the shape of the graph behind an external
pointer.

## Usage

``` r
rs_graph_info(ptr)
```

## Arguments

- ptr:

  External pointer to the Rust graph.

## Value

A list with

- n_nodes - Integer. Number of nodes.

- n_edges - Numeric. Number of (undirected) edges. Numeric, since a very
  large graph can hold more than `.Machine$integer.max` edges.

- directed - Logical. Whether the graph is directed.
