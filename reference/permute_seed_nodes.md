# Generate permuation scores for the diffusion

This function generate node-degree adjusted permutations of a given
diffusion score and adds Z-scores to the object. The function will
automatically determine if the original diffusion was a single or tied
diffusion and construct permutations accordingly.

## Usage

``` r
permute_seed_nodes(
  object,
  perm_iters = 1000L,
  random_seed = 10101L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `NetworkDiffusions` object. The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

- perm_iters:

  Integer. Number of permutations to test for. Defaults to `1000L`.

- random_seed:

  Integer. Random seed for determinism.

- .verbose:

  Boolean. Controls verbosity.

## Value

The class with added diffusion score based on a single set of seed
genes. Additionally, the seed genes are stored in the class.

## Examples

``` r
# 100 node-degree adjusted permutations of a single diffusion
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
