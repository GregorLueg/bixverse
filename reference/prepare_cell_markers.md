# Helper function to prepare cell markers

This function is a helper to generate a list of cell type to cell marker
annotations that can be fed into for example `[REF]`.

## Usage

``` r
prepare_cell_markers(obj, marker_df)
```

## Arguments

- obj:

  A single cell class, i.e. one of `SingleCells` or
  `SingleCellsMultiModal`.

- marker_df:

  A data.table with the columns `"cell_type"` and `"gene_id"`. You need
  to ensure that the gene id type matches the obj gene id.

## Value

A list of cell type to marker associations ready for subsequent usage.
Genes not found in the object will be automatically removed.
