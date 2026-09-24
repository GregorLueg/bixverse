# Calculates the compactness of the MetaCells based on diffusion map coordinates

**\[experimental\]** Calculates the meta cell compactness, i.e. the
average variance across the diffusion components over the cells of each
meta cell. Lower is better.

## Usage

``` r
rs_metacell_compactness(dc, meta_cells)
```

## Arguments

- dc:

  Numerical matrix. The diffusion map coordinates, cells x components.

- meta_cells:

  List. Per meta cell, an integer vector with the row indices
  (1-indexed!) into `dc`.

## Value

Numerical vector with one compactness value per meta cell. Empty meta
cells yield `NaN`.
