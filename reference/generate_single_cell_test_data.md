# Single cell test data

This function generates synthetic data for single cell test purposes.
These data can be used for testing functionality of various single cell
functions.

## Usage

``` r
generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(),
  seed = 42L
)
```

## Arguments

- syn_data_params:

  List. Contains the parameters for the generation of synthetic data,
  see: [`params_sc_synthetic_data()`](params_sc_synthetic_data.md). Has
  the following elements:

  - n_cells - Integer. Number of cells.

  - n_genes - Integer. Number of genes.

  - marker_genes - List. A nested list that indicates which gene indices
    are markers for which cell.

  - n_batches - Integer. Number of batches.

  - batch_effect_strength - String. Indicates the strength of the batch
    effect to add.

- seed:

  Integer. The seed for the generation of the seed data.

## Value

List with the following items

- counts - dgRMatrix with cells x genes.

- obs - data.table that contains the cell information.

- var - data.table that contains the var information.
