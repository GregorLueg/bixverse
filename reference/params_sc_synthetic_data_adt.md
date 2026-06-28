# Default parameters for generation of synthetic single cell data (ADT)

For the generation of synthetic single cell data mostly for testing or
showcasing purposes. In this case, ADT counts to test multi-modal
integration. The default configurations generates 1000 cells x 15
proteins with probes 1:3 being cell markers for cell type 1, genes 4:6
for cell type 2 and genes 7:9 for cell type. Columns 13:15 represents
isotype controls.

## Usage

``` r
params_sc_synthetic_data_adt(
  n_cells = 1000L,
  n_proteins = 15L,
  n_batches = 1L,
  marker_genes = list(cell_type_1 = list(marker_genes = 0:2L), cell_type_2 =
    list(marker_genes = 3:5L), cell_type_3 = list(marker_genes = 6:8L)),
  isotype_controls = 12L:14L,
  batch_effect_strength = c("strong", "medium", "weak")
)
```

## Arguments

- n_cells:

  Integer. Number of cells.

- n_proteins:

  Integer. Number of proteins

- n_batches:

  Integer. Number of batches.

- marker_genes:

  List. A nested list that indicates which gene indices are markers for
  which cell.

- isotype_controls:

  Integer vector. The columns that defines the isotype controls.
  (0-indexed!)

- batch_effect_strength:

  String. One of `c("strong", "medium", "weak")`. The strength of the
  batch effect to add.

## Value

A list with the parameters.
