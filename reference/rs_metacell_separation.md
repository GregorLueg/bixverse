# Calculates the separation of the centroids of the MetaCells based on diffusion map coordinates.

**\[experimental\]** Calculates the separation, i.e. the Euclidean
distance from each meta cell centroid in diffusion space to the nearest
other meta cell centroid. Higher is better.

## Usage

``` r
rs_metacell_separation(dc, meta_cells)
```

## Arguments

- dc:

  Numerical matrix. The diffusion map coordinates, cells x components.

- meta_cells:

  List. Per meta cell, an integer vector with the row indices
  (1-indexed!) into `dc`.

## Value

Numerical vector with one separation value per meta cell. Empty meta
cells yield `NaN`; a lone non-empty meta cell yields `Inf`.
