# Constrained personalised page rank

This function implements a constrained personalised page-rank. You can
define `sink_nodes` or `sink_edges`. In the case of the former, the vist
of a sink_node (via a type attribute in the igraph) will automatically
cause the surfer to reset. In the case of the latter, the traversal of a
sink_edge is allowed, however, subsequently, the surfer will be reseted.
This function can be useful for deriving PPR profiles under constraints
in heterogenous graphs and has been inspired by Ruiz, et al.

## Usage

``` r
constrained_page_rank(
  graph,
  personalisation_vector,
  sink_nodes = NULL,
  sink_edges = NULL
)
```

## Arguments

- graph:

  igraph. This one needs to be directed and weighted and have the node
  attribute `type` defining the node type and the edge attribute `type`
  defining the edge type.

- personalisation_vector:

  Numerical vector. The personalisation vector for the reset of the
  personalised page rank. The function will automatically normalise it.

- sink_nodes:

  Optional String vector. The node types that should force a reset.

- sink_edges:

  Optional String vector. The edge types after which there should be a
  forced reset.

## Value

The constrained personalised page rank

## References

Ruiz, et al., Nat Commun, 2021
