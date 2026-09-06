# Identify privileged communities based on a given diffusion vector

Detects privileged communities after a diffusion based on seed nodes.

## Usage

``` r
community_detection(
  object,
  community_params = params_community_detection(),
  seed = 42L,
  .verbose = FALSE,
  .max_iters = 100L
)
```

## Arguments

- object:

  `NetworkDiffusions` object. The underlying class
  [`NetworkDiffusions()`](https://gregorlueg.github.io/bixverse/reference/NetworkDiffusions.md).

- community_params:

  List. Parameters for the community detection within the reduced
  network, see
  [`params_community_detection()`](https://gregorlueg.github.io/bixverse/reference/params_community_detection.md).
  A list with the following items:

  - max_nodes - Integer. Number of maximum nodes per community. Larger
    communities will be recursively subclustered.

  - min_nodes - Integer. Minimum number of nodes per community.

  - min_seed_nodes - Integer. Minimum number of seed genes that have to
    be found in a given community.

  - initial_res - Float. Initial resolution parameter for the Leiden
    clustering.

  - threshold_type - String. One of `c("prop_based", "pval_based")`. You
    can chose to include a certain proportion of the network (like in
    the original paper) with the highest diffusion scores, or use
    p-values based on permutations. Defaults to `"prop_based"`.

  - network_threshold - Float. The proportion of the network to include.
    Used if `threshold_type = "prop_based"`.

  - pval_threshold - Float. The maximum p-value for nodes to be
    included. Used if `threshold_type = "pval_based"`.

- seed:

  Random seed.

- .verbose:

  Controls the verbosity of the function.

- .max_iters:

  Controls how many iterations shall be tried for the sub-clustering. To
  note, in each iteration of the sub-clustering, the resolution
  parameter is increased by 0.05, to identify more granular communities
  within the sub communities.

## Value

The class with added diffusion community detection results (if any could
be identified with the provided parameters).

## Examples

``` r
# privileged communities in the heated part of the network
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
