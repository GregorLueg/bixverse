# Get correlation-based graph

Helper function to get a correlation-based igraph from the class

## Usage

``` r
get_cor_graph(object, epsilon, .verbose)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- epsilon:

  Float. The epsilon parameter for the RBF function, in this case the
  bump function.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the following elements:

- graph - The igraph

- params - A list that contains the parameters of the graph generation
  and general graph information (node, edge numbers).

## Examples

``` r
# igraph from the stored correlations at a fixed epsilon
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- cor_module_processing(obj, cor_method = "spearman", .verbose = FALSE)
graph_res <- get_cor_graph(obj, epsilon = 2, .verbose = FALSE)
graph_res$params
#> $epsilon
#> [1] 2
#> 
#> $no_nodes
#> [1] 217
#> 
#> $no_edges
#> [1] 4833
#> 
```
