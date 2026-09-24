# Wrapper function to generate QC metric params for single cell

Wrapper function to generate QC metric params for single cell

## Usage

``` r
params_sc_min_quality(
  min_unique_genes = 100L,
  min_lib_size = 250L,
  min_cells = 10L,
  target_size = 10000
)
```

## Arguments

- min_unique_genes:

  Integer. Minimum number of unique genes per cell/spot to be included.
  Defaults to `100L`.

- min_lib_size:

  Integer. Minimum library size per cell/spot to be included. Defaults
  to `250L`.

- min_cells:

  Integer. Minimum number of cells a gene has to be expressed to be
  included. Defaults to `10L`.

- target_size:

  Numeric. The target size for the normalisation. Defaults to `10000.0`.

## Value

A named list with the following elements:

- min_unique_genes - Integer. Minimum number of unique genes per
  cell/spot to be included. Defaults to `100L`.

- min_lib_size - Integer. Minimum library size per cell/spot to be
  included. Defaults to `250L`.

- min_cells - Integer. Minimum number of cells a gene has to be
  expressed to be included. Defaults to `10L`.

- target_size - Numeric. The target size for the normalisation. Defaults
  to `10000.0`.
