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

  `network_diffusions` object. The underlying class
  [`network_diffusions()`](network_diffusions.md).

- perm_iters:

  Integer. Number of permutations to test for. Defaults to `1000L`.

- random_seed:

  Integer. Random seed for determinism.

- .verbose:

  Boolean. Controls verbosity.

## Value

The class with added diffusion score based on a single set of seed
genes. Additionally, the seed genes are stored in the class.
