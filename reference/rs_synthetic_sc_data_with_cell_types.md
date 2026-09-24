# Generates synthetic data for single cell

**\[experimental\]** Helper function to generate synthetic single cell
data with optional batch effects and sample bias.

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

  Integer. Number of batches to generate.

- n_samples:

  Optional integer. Shall the cells be distributed over `n_samples`
  samples. Only used together with `sample_bias`.

- cell_configs:

  List. One element per cell type, each a list with a `marker_genes`
  integer vector of 0-based marker gene indices.

- batch_effect_strength:

  String. One of `c("strong", "medium", "weak")`. Defines the strength
  of the added batch effect. Unknown values fall back to `"strong"`.

- sample_bias:

  Optional string. One of `c("even", "slightly_uneven", "very_uneven")`.
  Other values raise an error.

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with the following items.

- data - The synthetic raw counts, CSR over cells.

- indptr - The index pointers of the cells.

- indices - The 0-based gene indices for the given cells.

- nrow - Number of cells.

- ncol - Number of genes.

- cell_type_indices - 0-based cell type per cell.

- batch_indices - 0-based batch per cell.

- sample_indices - 0-based sample per cell. `NULL` unless both
  `n_samples` and `sample_bias` are provided.
