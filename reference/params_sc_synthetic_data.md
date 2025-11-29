# Default parameters for generation of synthetic data

Default parameters for generation of synthetic data

## Usage

``` r
params_sc_synthetic_data(
  n_cells = 1000L,
  n_genes = 100L,
  n_batches = 1L,
  marker_genes = list(cell_type_1 = list(marker_genes = 0:9L), cell_type_2 =
    list(marker_genes = 10:19L), cell_type_3 = list(marker_genes = 20:29L)),
  batch_effect_strength = c("strong", "medium", "weak"),
  n_samples = NULL,
  sample_bias = NULL
)
```

## Arguments

- n_cells:

  Integer. Number of cells.

- n_genes:

  Integer. Number of genes.

- n_batches:

  Integer. Number of batches.

- marker_genes:

  List. A nested list that indicates which gene indices are markers for
  which cell.

- batch_effect_strength:

  String. One of `c("strong", "medium", "weak")`. The strength of the
  batch effect to add.

- n_samples:

  Optional integer. Shall sample membership be added to the synthetic
  data. If you want sample information you need to provide `n_samples`
  and `sample_bias`.

- sample_bias:

  Optional string. One of `c("even", "slightly_uneven", "very_uneven")`

## Value

A list with the parameters.
