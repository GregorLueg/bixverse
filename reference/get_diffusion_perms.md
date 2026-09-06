# Get the diffusion permutations

Returns the diffusion Z-scores if you ran
[`permute_seed_nodes()`](https://gregorlueg.github.io/bixverse/reference/permute_seed_nodes.md).

## Usage

``` r
get_diffusion_perms(object)
```

## Arguments

- object:

  The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

## Value

The diffusion Z scores if found. Otherwise `NULL`.

## Examples

``` r
# Z-scores from node-degree adjusted permutations
set.seed(42)
g <- igraph::sample_pa(15, directed = FALSE)
edges <- data.table::setDT(igraph::as_data_frame(g))[, `:=`(
  from = sprintf("node_%i", from),
  to = sprintf("node_%i", to)
)]
object <- NetworkDiffusions(edges, weighted = FALSE, directed = FALSE)
object <- diffuse_seed_nodes(object, c(node_1 = 1, node_3 = 1), "max")
object <- permute_seed_nodes(object, perm_iters = 100L, .verbose = FALSE)
head(get_diffusion_perms(object))
#>     node_1     node_2     node_3     node_4     node_6     node_7 
#>  2.7710933  1.9146426  1.9256316 -0.3677003 -0.7718167 -0.6526446 
```
