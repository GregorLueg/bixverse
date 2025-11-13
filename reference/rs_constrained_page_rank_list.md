# Calculate a constrained page-rank score over a list.

This function can be used to get constrainted personalised page-rank
scores akin to Ruiz, et al. You can provide optionally `sink_nodes`
(node types that will force a reset) and/or `sink_edges` (edge types
that will force a reset). This version can take in a list of
personalisation vectors and returns a list as result.

## Usage

``` r
rs_constrained_page_rank_list(
  personalisation_list,
  node_names,
  node_types,
  from,
  to,
  weights,
  edge_type,
  sink_nodes,
  sink_edges
)
```

## Arguments

- personalisation_list:

  List. The list with the personalisation vectors. The sum must equal to
  1, otherwise the function panics!

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

- sink_nodes:

  Optional string vector. Should these node types be seen as sinks,
  i.e., the reset occurs when this node is reached.

- sink_edges:

  Optional string vector. Shall an automatic reset occur when this edge
  type is traversed.

## Value

A list of the personalised (constrained) page rank values.

## References

Ruiz, et al., Nat Commun, 2021
