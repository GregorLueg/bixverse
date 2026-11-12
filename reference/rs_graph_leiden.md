# Leiden and Louvain clustering on a Rust neighbour graph

**\[experimental\]** Same engine and parameters as
[`rs_leiden_clustering()`](https://gregorlueg.github.io/bixverse/reference/rs_leiden_clustering.md),
but on the graph behind the pointer, so nothing is flattened into R
vectors. Louvain is `refine = FALSE`.

## Usage

``` r
rs_graph_leiden(ptr, leiden_params)
```

## Arguments

- ptr:

  External pointer to the Rust graph.

- leiden_params:

  List. See
  [`rs_leiden_clustering()`](https://gregorlueg.github.io/bixverse/reference/rs_leiden_clustering.md).

## Value

Integer vector with the 1-based community of each node; `1` is the
largest community.
