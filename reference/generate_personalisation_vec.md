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

## Examples

``` r
# reset weight split across two seed nodes, normalised to sum to one
g <- igraph::graph_from_data_frame(
  data.frame(from = c("a", "b", "c"), to = c("b", "c", "d")),
  directed = TRUE
)
generate_personalisation_vec(g, node_weights = c(a = 3, c = 1))
#>    a    b    c    d 
#> 0.75 0.00 0.25 0.00 
```
