# Single cell test data (ADT)

This function generates synthetic ADT counts for single cell test
purposes. These data can be used for testing multi-modal functionality
of various single cell functions. Pairs cell-for-cell with
[`generate_single_cell_test_data()`](https://gregorlueg.github.io/bixverse/reference/generate_single_cell_test_data.md)
when generated with matching `n_cells`, number of cell types and
`n_batches`.

## Usage

``` r
generate_single_cell_test_data_adt(
  syn_data_params = params_sc_synthetic_data_adt(),
  seed = 42L
)
```

## Arguments

- syn_data_params:

  List. Contains the parameters for the generation of synthetic ADT
  data, see:
  [`params_sc_synthetic_data_adt()`](https://gregorlueg.github.io/bixverse/reference/params_sc_synthetic_data_adt.md).
  Has the following elements:

  - n_cells - Integer. Number of cells.

  - n_proteins - Integer. Number of proteins.

  - marker_genes - List. A nested list that indicates which protein
    indices are markers for which cell type.

  - n_batches - Integer. Number of batches.

  - isotype_controls - Integer vector. Column indices (0-based) of the
    isotype controls.

  - batch_effect_strength - String. Indicates the strength of the batch
    effect to add.

- seed:

  Integer. The seed for the generation of the synthetic data.

## Value

List with the following items

- counts - Numeric matrix with cells x proteins.

- obs - data.table that contains the cell information.

- var - data.table that contains the protein information.
