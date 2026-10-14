# Meta cell compactness in a diffusion map

Scores a meta cell assignment for compactness, i.e. the mean variance
across the diffusion components over the cells of each meta cell. Lower
is better. Unlike
[`calc_manifold_metrics()`](https://gregorlueg.github.io/bixverse/reference/calc_manifold_metrics.md),
this works on any assignment, so meta cells generated outside bixverse
(or on a different object) can be scored against an existing diffusion
map.

## Usage

``` r
calc_meta_cell_compactness(dcs, assignment)
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

Named numeric vector with one compactness value per meta cell, named and
ordered like [`split()`](https://rdrr.io/r/base/split.html) orders
`assignment`.

## References

Persad, et al., Nat Biotechnol, 2023

## Examples

``` r
set.seed(1L)
dcs <- matrix(rnorm(200L), ncol = 2L)
assignment <- rep(sprintf("mc_%i", 1:10), each = 10L)
calc_meta_cell_compactness(dcs, assignment)
#>      mc_1     mc_10      mc_2      mc_3      mc_4      mc_5      mc_6      mc_7 
#> 0.6694486 0.8578757 0.6948285 0.5725448 0.5120067 0.9328036 0.9981734 1.2029001 
#>      mc_8      mc_9 
#> 0.5507736 0.7109698 
```
