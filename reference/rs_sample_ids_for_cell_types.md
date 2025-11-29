# Helper function to generate sample identifiers based on cells

Extract out of
[`rs_synthetic_sc_data_with_cell_types()`](rs_synthetic_sc_data_with_cell_types.md)
to quickly iterate over different sample to cell type patterns

## Usage

``` r
rs_sample_ids_for_cell_types(cell_type_indices, n_samples, sample_bias, seed)
```

## Arguments

- cell_type_indices:

  Integer vector. Each integer represents a cell type.

- n_samples:

  Integer. Number of different sample ids to generate.

- sample_bias:

  String. One of `c("even", "slightly_uneven", "very_uneven")`.
  Determins the cell type to sample id associations.

- seed:

  Integer. Random seed for reproducibility.

## Value

An integer vector representing the samples.
