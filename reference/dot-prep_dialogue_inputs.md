# Resolve the per-cell-type inputs DIALOGUE needs

Turns the obs table plus the user's feature matrices into the flat,
0-indexed vectors the Rust side takes. The cell types analysed are
`names(features)`, in that order: anything else in `cell_type_col` is
ignored, and that order is what the pair columns of `emp_p` follow.

## Usage

``` r
.prep_dialogue_inputs(
  object,
  obs,
  cell_type_col,
  sample_col,
  features,
  quality_col,
  cell_id_col,
  default_quality
)
```

## Arguments

- object:

  `SingleCells`, `SingleCellsSubset` or `MetaCells` class.

- obs:

  data.table. The observation table to read the labels off.

- cell_type_col:

  String. Column holding the cell type labels.

- sample_col:

  String. Column holding the sample labels.

- features:

  Named list of numeric matrices, one per cell type.

- quality_col:

  Optional string. Column to use as the quality covariate.

- cell_id_col:

  String. Column holding the cell identifiers.

- default_quality:

  Numeric vector aligned to `obs` rows. Used when `quality_col` is
  `NULL`, after a `log1p` and a z-score. Typically the library size.

## Value

A list with the following items

- cell_type_indices - List of integer vectors. 0-indexed global cell
  positions per cell type.

- features - List of numeric matrices, rows reordered to match.

- sample_ids - Integer vector. 0-indexed sample code per global cell.

- cell_quality - Numeric vector. Quality covariate per global cell.

- cell_types - Character vector. The analysed cell types, in order.

- cell_ids - List of character vectors. Cell ids per cell type.

- sample_levels - Character vector. Sample labels, position `i` being
  sample code `i - 1`.

## Details

`sample_ids` and `cell_quality` are indexed by *global* cell position,
not by position within a cell type, so both are scattered into vectors
long enough to cover the largest index in use. Cells belonging to no
analysed cell type leave holes, which Rust never reads.
