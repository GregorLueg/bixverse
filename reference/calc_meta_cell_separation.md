# Meta cell separation in a diffusion map

Scores a meta cell assignment for separation, i.e. the Euclidean
distance from each meta cell centroid in diffusion space to the nearest
other meta cell centroid. Higher is better. Unlike
[`calc_manifold_metrics()`](https://gregorlueg.github.io/bixverse/reference/calc_manifold_metrics.md),
this works on any assignment, so meta cells generated outside bixverse
(or on a different object) can be scored against an existing diffusion
map.

## Usage

``` r
calc_meta_cell_separation(dcs, assignment)
```

## Arguments

- dcs:

  Numeric matrix. The diffusion map coordinates, cells x components, for
  example `S7::prop(object, "other_data")$dcs` after
  [`calc_diffusion_coordinates()`](https://gregorlueg.github.io/bixverse/reference/calc_diffusion_coordinates.md).
  Its rows follow the `used_cells` of the kNN object the diffusion map
  was built from.

- assignment:

  Atomic vector of length `nrow(dcs)` with the meta cell of each row of
  `dcs`. `NA` marks cells that belong to no meta cell. Overlapping meta
  cells (hdWGCNA-style) cannot be expressed this way; use
  [`calc_manifold_metrics()`](https://gregorlueg.github.io/bixverse/reference/calc_manifold_metrics.md)
  for those.

## Value

Named numeric vector with one separation value per meta cell, named and
ordered like [`split()`](https://rdrr.io/r/base/split.html) orders
`assignment`. A single meta cell yields `Inf`.

## References

Persad, et al., Nat Biotechnol, 2023

## Examples

``` r
set.seed(1L)
dcs <- matrix(rnorm(200L), ncol = 2L)
assignment <- rep(sprintf("mc_%i", 1:10), each = 10L)
calc_meta_cell_separation(dcs, assignment)
#>       mc_1      mc_10       mc_2       mc_3       mc_4       mc_5       mc_6 
#> 0.34729806 0.03480481 0.13514426 0.28378615 0.12786378 0.12786378 0.07421945 
#>       mc_7       mc_8       mc_9 
#> 0.34729806 0.55519325 0.03480481 
```
