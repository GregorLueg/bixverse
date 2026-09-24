# Calculates diffusion maps for density calculations for meta cells

**\[experimental\]** Builds multiscale diffusion components from the kNN
graph and uses the distance to the `k_density`-th neighbour in that
space as a density proxy. The lower quartile of these distances is
tagged high density, the upper quartile low density, the rest mid.

## Usage

``` r
rs_metacell_density(knn_data, n_dcs, k_density, knn_params, verbose, seed)
```

## Arguments

- knn_data:

  Named list. The kNN data with `indices` (0-indexed!), `dist`, `k` and
  `dist_metric`.

- n_dcs:

  Integer. The number of diffusion coordinates to return. Typically
  `10`.

- k_density:

  Integer. The k-nearest neighbour to use for the density estimation.
  Typically `150`.

- knn_params:

  List. The kNN parameters defined by
  [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md),
  used for the search in diffusion space.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. For reproducibility.

## Value

A list with the following items

- dcs - Numerical matrix of cells x `n_dcs` with the multiscale
  diffusion components.

- density_distances - Numerical vector. Distance to the `k_density`-th
  neighbour in diffusion space per cell.

- regions - Character vector. `"high"`, `"mid"` or `"low"` density per
  cell.

## References

Persad, et al., Nat. Biotechnol., 2023.
