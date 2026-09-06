# Get the final results from the class

Get the final results from `BixverseBaseClass` class (or child classes).

## Usage

``` r
get_results(object)

get_results.DialogueResult(object, ...)
```

## Arguments

- object:

  The underlying
  [`BixverseBaseClass()`](https://gregorlueg.github.io/bixverse/reference/BixverseBaseClass.md)
  class. The class functionality is usually inherited by other S7
  classes in `bixverse`.

- ...:

  Unused, present so the S3 methods sharing this page match the generic.

## Value

Returns the final results if any have been stored in the class.

## Examples

``` r
# communities found after a network diffusion
set.seed(42)
g <- igraph::sample_pa(15, directed = FALSE)
edges <- data.table::setDT(igraph::as_data_frame(g))[, `:=`(
  from = sprintf("node_%i", from),
  to = sprintf("node_%i", to)
)]
object <- NetworkDiffusions(edges, weighted = FALSE, directed = FALSE)
object <- diffuse_seed_nodes(object, c(node_1 = 1, node_3 = 1), "max")
object <- permute_seed_nodes(object, perm_iters = 100L, .verbose = FALSE)
object <- community_detection(
  object,
  community_params = params_community_detection(
    min_seed_nodes = 0L,
    min_nodes = 2L
  )
)
head(get_results(object))
#>    cluster_id node_id    ks_pval cluster_size seed_nodes_no diffusion_score
#>        <char>  <char>      <num>        <int>         <int>           <num>
#> 1:  cluster_1  node_1 0.00965701            5             1      0.12792857
#> 2:  cluster_1  node_2 0.00965701            5             1      0.24907561
#> 3:  cluster_1  node_5 0.00965701            5             1      0.07379672
#> 4:  cluster_1 node_13 0.00965701            5             1      0.04910154
#> 5:  cluster_1 node_10 0.00965701            5             1      0.05292857
#> 6:  cluster_2  node_3 0.16923077            3             1      0.15054927
#>    seed_node
#>       <lgcl>
#> 1:      TRUE
#> 2:     FALSE
#> 3:     FALSE
#> 4:     FALSE
#> 5:     FALSE
#> 6:      TRUE
```
