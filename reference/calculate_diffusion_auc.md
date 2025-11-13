# Calculate the AUROC for a diffusion score

This functions can take a given `network_diffusions` object and
calculates an AUC and generates a Z-score based on random permutation of
`random_aucs` for test for statistical significance if desired.

## Usage

``` r
calculate_diffusion_auc(
  object,
  hit_nodes,
  auc_iters = 10000L,
  random_aucs = 1000L,
  permutation_test = FALSE,
  seed = 42L
)
```

## Arguments

- object:

  `network_diffusions` object. The underlying class
  [`network_diffusions()`](network_diffusions.md).

- hit_nodes:

  String vector. Which nodes in the graph are considered a 'hit'.

- auc_iters:

  Integer. How many iterations to run to approximate the AUROC.

- random_aucs:

  Integer. How many random AUROCs to calculate to estimate the Z-score.
  Only of relevance if permutation test is set to `TRUE`.

- permutation_test:

  Boolean. Shall a permutation based Z-score be calculated.

- seed:

  Integer. Random seed.

## Value

List with AUC and Z-score as the two named elements if permutations test
set to TRUE; otherwise just the AUC.
