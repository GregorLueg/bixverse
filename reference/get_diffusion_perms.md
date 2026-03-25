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
