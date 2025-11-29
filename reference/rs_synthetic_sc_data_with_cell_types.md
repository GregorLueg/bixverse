# Generates synthetic data for single cell

Helper function to generate synthetic single cell data with optional

## Usage

``` r
rs_synthetic_sc_data_with_cell_types(
  n_cells,
  n_genes,
  n_batches,
  n_samples,
  cell_configs,
  batch_effect_strength,
  sample_bias,
  seed
)
```

## Arguments

- n_cells:

  Integer. Number of cells to generate.

- n_genes:

  Integer. Number of genes to generate.

- n_batches:

  Integer. Number of the batches to generated.

- n_samples:

  Optional integer. Shall the cells be distributed over `n_samples`
  samples.

- cell_configs:

  A nested list that indicates which gene indices are markers for which
  cell.

- batch_effect_strength:

  String. One of `c("strong", "medium", "low")`. Defines the strength of
  the added batch effect.

- sample_bias:

  Optional string. One of `c("even", "slightly_uneven", "very_uneven")`

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with the following items.

- data - The synthetic raw counts.

- indptr - The index pointers of the cells.

- indices - The indices of the genes for the given cells.

- nrow - Number of rows.

- ncol - Number of columns

- cell_type_indices - Vector indicating which cell type this is.

- batch_indices - Vector indicating the batch.

- sample_indices - Optional sample indices if asked for.
