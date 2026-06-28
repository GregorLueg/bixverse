# Generates synthetic ADT counts with defined cell types

**\[experimental\]** Generates a dense cells x proteins matrix of
synthetic raw ADT counts for testing. Proteins are assigned roles by
column index: marker proteins are elevated in their owning cell type and
sit at background elsewhere, isotype controls only ever carry
background, and any column named as neither is a generic background-only
protein. Counts follow a negative-binomial draw with an additive
background plus per-cell-type signal, a per-cell capture efficiency
factor, and an optional per-batch staining multiplier. Cell type and
batch assignment match `rs_synthetic_sc_with_cell_types()` cell-for-cell
for matched inputs, so RNA and ADT can be paired for multi-modal tests.

## Usage

``` r
rs_synthetic_sc_adt_with_cell_types(
  n_cells,
  n_proteins,
  n_batches,
  isotype_controls,
  cell_configs,
  batch_effect_strength,
  seed
)
```

## Arguments

- n_cells:

  Integer. Number of cells (matrix rows).

- n_proteins:

  Integer. Total number of proteins (matrix columns). Must be large
  enough to cover every marker and isotype index supplied.

- n_batches:

  Integer. Number of batches. Batch 0 is unperturbed; further batches
  receive a per-protein staining multiplier.

- isotype_controls:

  Integer vector. The 0-based column indices that are isotype controls.
  These are forced to background only, even if they also appear in a
  cell type's markers.

- cell_configs:

  List. One element per cell type, each a list with a `marker_genes`
  integer vector of 0-based marker column indices for that cell type.

- batch_effect_strength:

  String. One of `c("weak", "medium", "strong")`. Controls the spread of
  the per-batch staining multiplier. Unrecognised values fall back to
  `"strong"`.

- seed:

  Integer. For reproducibility.

## Value

A list with the following items.

- data - Integer vector. The counts in row-major order, length
  `n_cells * n_proteins` (cell-major: all proteins of cell 0, then cell
  1).

- cell_type_indices - Integer vector of length `n_cells`. The 0-based
  cell type assigned to each cell.

- batch_indices - Integer vector of length `n_cells`. The 0-based batch
  assigned to each cell.
