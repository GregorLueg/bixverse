# Get the diffusion permutations

Returns the diffusion Z-scores if you ran
[`permute_seed_nodes()`](permute_seed_nodes.md).

## Usage

``` r
get_diffusion_perms(object)
```

## Arguments

- object:

  The underlying class [`network_diffusions()`](network_diffusions.md).

## Value

The diffusion Z scores if found. Otherwise `NULL`.
