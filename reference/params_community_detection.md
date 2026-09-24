# Wrapper function to generate community detection parameters

Wrapper function to generate community detection parameters

## Usage

``` r
params_community_detection(
  max_nodes = 300L,
  min_nodes = 10L,
  min_seed_nodes = 2L,
  initial_res = 0.5,
  threshold_type = c("prop_based", "pval_based"),
  network_threshold = 0.5,
  pval_threshold = 0.1
)
```

## Arguments

- max_nodes:

  Integer. Maximum number of nodes in a given community. Defaults to
  `300L`.

- min_nodes:

  Integer. Minimum number of nodes in a given community. Defaults to
  `10L`.

- min_seed_nodes:

  Integer. Minimum number of seed nodes within a community. Defaults to
  `2L`.

- initial_res:

  Numeric. Initial resolution parameter to start with. Defaults to
  `0.5`.

- threshold_type:

  String. You can chose to include a certain proportion of the network
  with the highest diffusion scores, or use p-values based on
  permutations. One of `c("prop_based", "pval_based")`. Defaults to
  `"prop_based"`.

- network_threshold:

  Numeric. The proportion of the network to include. Used if
  `threshold_type = "prop_based"`. Defaults to `0.5`.

- pval_threshold:

  Numeric. The maximum p-value for nodes to be included. Used if
  `threshold_type = "pval_based"`. Defaults to `0.1`.

## Value

A named list with the following elements:

- max_nodes - Integer. Maximum number of nodes in a given community.
  Defaults to `300L`.

- min_nodes - Integer. Minimum number of nodes in a given community.
  Defaults to `10L`.

- min_seed_nodes - Integer. Minimum number of seed nodes within a
  community. Defaults to `2L`.

- initial_res - Numeric. Initial resolution parameter to start with.
  Defaults to `0.5`.

- threshold_type - String. You can chose to include a certain proportion
  of the network with the highest diffusion scores, or use p-values
  based on permutations. One of `c("prop_based", "pval_based")`.
  Defaults to `"prop_based"`.

- network_threshold - Numeric. The proportion of the network to include.
  Used if `threshold_type = "prop_based"`. Defaults to `0.5`.

- pval_threshold - Numeric. The maximum p-value for nodes to be
  included. Used if `threshold_type = "pval_based"`. Defaults to `0.1`.
