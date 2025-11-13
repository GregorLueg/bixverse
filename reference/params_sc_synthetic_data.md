# Default parameters for generation of synthetic data

Default parameters for generation of synthetic data

## Usage

``` r
params_sc_synthetic_data(
  n_cells = 1000L,
  n_genes = 100L,
  marker_genes = list(cell_type_1 = list(marker_genes = 0:9L), cell_type_2 =
    list(marker_genes = 10:19L), cell_type_3 = list(marker_genes = 20:29L)),
  n_batches = 1L,
  batch_effect_strength = c("strong", "medium", "weak")
)
```

## Arguments

- n_cells:

  Integer. Number of cells.

- n_genes:

  Integer. Number of genes.

- marker_genes:

  List. A nested list that indicates which gene indices are markers for
  which cell.

- n_batches:

  Integer. Number of batches.

- batch_effect_strength:

  String. One of `c("strong", "medium", "weak")`. The strength of the
  batch effect to add.

## Value

A list with the parameters.
