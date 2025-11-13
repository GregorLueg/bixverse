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

  `network_diffusions` object. The underlying class
  [`network_diffusions()`](network_diffusions.md).

- community_params:

  List. Parameters for the community detection within the reduced
  network, see
  [`params_community_detection()`](params_community_detection.md). A
  list with the following items:

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
