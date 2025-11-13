# Calculate a constrained page-rank score

This function can be used to get constrainted personalised page-rank
scores akin to Ruiz, et al. You can provide optionally `sink_nodes`
(node types that will force a reset) and/or `sink_edges` (edge types
that will force a reset).

## Usage

``` r
rs_constrained_page_rank(
  node_names,
  node_types,
  from,
  to,
  weights,
  edge_type,
  personalised,
  sink_nodes,
  sink_edges
)
```

## Arguments

- node_names:

  String vector. Name of the graph nodes.

- node_types:

  String vector. The node types.

- from:

  String vector. The names of the `from` edges from the edge list.

- to:

  String vector. The names of the `to` edges from the edge list.

- weights:

  Numerical vector. The edge weights from the edge list.

- edge_type:

  String vector. The edge types.

- personalised:

  Numerical vector. The reset values. They must sum to 1 and be of same
  length of `node_names`!

- sink_nodes:

  Optional string vector. Should these node types be seen as sinks,
  i.e., the reset occurs when this node is reached.

- sink_edges:

  Optional string vector. Shall an automatic reset occur when this edge
  type is traversed.

## Value

The personalised constrained page rank values.

## References

Ruiz, et al., Nat Commun, 2021
