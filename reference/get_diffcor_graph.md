# Get differential correlation-based graph

Helper function to get a differential correlation-based igraph from the
class

## Usage

``` r
get_diffcor_graph(object, min_cor = 0.2, fdr_threshold = 0.05, .verbose = TRUE)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

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

## Examples

``` r
# igraph from the differential correlations
sig <- synthetic_signal_matrix()
mat <- t(sig$mat)
target <- mat[sig$group %in% c("group1", "group2"), ]
background <- mat[sig$group == "group3", ]
meta <- data.table::data.table(sample_id = rownames(target))
obj <- BulkCoExp(target, meta)
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- diffcor_module_processing(
  obj, background, cor_method = "pearson", .verbose = FALSE
)
graph_res <- get_diffcor_graph(obj, .verbose = FALSE)
graph_res$params
#> $min_cor
#> [1] 0.2
#> 
#> $fdr_threshold
#> [1] 0.05
#> 
#> $no_nodes
#> [1] 279
#> 
#> $no_edges
#> [1] 6477
#> 
```
