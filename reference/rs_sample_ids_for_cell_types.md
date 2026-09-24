# Helper function to generate sample identifiers based on cells

**\[experimental\]** Extract out of
[`rs_synthetic_sc_data_with_cell_types()`](https://gregorlueg.github.io/bixverse/reference/rs_synthetic_sc_data_with_cell_types.md)
to quickly iterate over different sample to cell type patterns

## Usage

``` r
rs_sample_ids_for_cell_types(cell_type_indices, n_samples, sample_bias, seed)
```

## Arguments

- cell_type_indices:

  Integer vector. Each integer represents a cell type (0-based, as
  returned by
  [`rs_synthetic_sc_data_with_cell_types()`](https://gregorlueg.github.io/bixverse/reference/rs_synthetic_sc_data_with_cell_types.md)).

- n_samples:

  Integer. Number of different sample ids to generate.

- sample_bias:

  String. One of `c("even", "slightly_uneven", "very_uneven")`.
  Determines the cell type to sample id associations. Other values raise
  an error.

- seed:

  Integer. Random seed for reproducibility.

## Value

An integer vector with the 0-based sample per cell.
