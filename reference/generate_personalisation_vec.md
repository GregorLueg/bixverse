# Helper function to create personalisation vectors

Helper function to create personalisation vectors

## Usage

``` r
generate_personalisation_vec(graph, node_weights)
```

## Arguments

- graph:

  igraph. The graph for which to produce the personalisation vector.

- node_weights:

  Named numeric. The names represent the nodes and the values the
  strength of the reset.

## Value

The personalisation vector for subsequent usage in page-rank
