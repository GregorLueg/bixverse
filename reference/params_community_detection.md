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

  Integer. Maximum number of nodes in a given community.

- min_nodes:

  Integer. Minimum number of nodes in a given community.

- min_seed_nodes:

  Integer. Minimum number of seed nodes within a community.

- initial_res:

  Float. Initial resolution parameter to start with.

- threshold_type:

  String. One of `c("prop_based", "pval_based")`. You can chose to
  include a certain proportion of the network with the highest diffusion
  scores, or use p-values based on permutations.

- network_threshold:

  Float. The proportion of the network to include. Used if
  `threshold_type = "prop_based"`.

- pval_threshold:

  Float. The maximum p-value for nodes to be included. Used if
  `threshold_type = "pval_based"`.

## Value

List with parameters for usage in subsequent function.
