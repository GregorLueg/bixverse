# Get differential correlation-based graph

Helper function to get a differential correlation-based igraph from the
class

## Usage

``` r
get_diffcor_graph(object, min_cor = 0.2, fdr_threshold = 0.05, .verbose = TRUE)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

- min_cor:

  Float. The minimum absolute correlation that needs to be present in
  either data set.

- fdr_threshold:

  Float. The maximum FDR that is tolerated for the generation of the
  graph.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A list with the following elements:

- graph - The igraph

- params - A list that contains the parameters of the graph generation
  and general graph information (node, edge numbers).
