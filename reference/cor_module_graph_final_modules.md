# Identify correlation-based gene modules via graphs

This function leverages graph-based clustering to identify gene
co-expression modules. The class has the option to sub-cluster large
communities within their respective sub graphs, akin to the approach
taken by Barrio-Hernandez, et al.

## Usage

``` r
cor_module_graph_final_modules(
  object,
  resolution = NULL,
  min_size = 10L,
  max_size = 500L,
  subclustering = TRUE,
  random_seed = 123L,
  .graph_params = params_cor_graph(),
  .max_iters = 100L,
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md).

- resolution:

  The Leiden resolution parameter you wish to use. If NULL, it will use
  the optimal one identified by
  [`cor_module_graph_check_res()`](cor_module_graph_check_res.md). If
  nothing can be found, will default to 1.

- min_size:

  Integer. Minimum size of the communities.

- max_size:

  Integer. Maximum size of the communities.

- subclustering:

  Boolean. Shall after a first clustering communities that are too large
  be further sub clustered. Defaults to `TRUE`.

- random_seed:

  Integer. Random seed.

- .graph_params:

  List. Parameters for the generation of the (differential) correlation
  graph, see [`params_cor_graph()`](params_cor_graph.md). Contains:

  - Epsilon - Defines the epsilon parameter for the radial basis
    function. Defaults to 2, but should be ideally optimised.

  - min_cor - Float. Minimum absolute correlation that needs to be
    observed in either data set. Only relevant for differential
    correlation-based graphs.

  - fdr_threshold - Float. Maximum FDR for the differential correlation
    p-value.

  - verbose - Boolean. Controls verbosity of the graph generation.

  This parameter is only relevant if you did *not* run
  [`cor_module_graph_check_res()`](cor_module_graph_check_res.md).

- .max_iters:

  Integer. If sub clustering is set to `TRUE`, what shall be the maximum
  number of iterations. Defaults to 100L.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

The class with added data to the properties.

## References

Barrio-Hernandez, et al., Nat Genet, 2023.
