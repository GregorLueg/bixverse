# Calculates diffusion maps for density calculations for meta cells

**\[experimental\]** Generates diffusion maps and identifies in which
density region a given cell sits (defined as distance to k-nearest
neighbours quite).

## Usage

``` r
rs_metacell_density(knn_data, n_dcs, k_density, knn_params, verbose, seed)
```

## Arguments

- knn_data:

  Named list. Needs to have the relevant data from the kNN graph.

- n_dcs:

  Integer. The number of diffusion coordinates to return. Typically
  `10`.

- k_density:

  Integer. The k-nearest neighbour to use for the density estimation.
  Typically `150`.

- knn_params:

  List. The kNN parameters defined by
  [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

- seed:

  Integer. For reproducibility.

## Value

A list with the following items

- dcs - Density coordinates

- density_distances - Density distances at `k_density` neighbours.

- regions - Region of the manifold where this given cell is.

## References

Persad, et al., Nat. Biotechnol., 2023.
