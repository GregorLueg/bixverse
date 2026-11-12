# Write a cells x features count matrix as a CellRanger style h5 file

**\[experimental\]** The cells x features CSR buffers are the features x
cells CSC buffers 10x stores, so they are written as they are. Counts,
indices, pointers and shape go out as `int32`, strings as
variable-length UTF-8. An existing file is truncated.

## Usage

``` r
rs_write_tenx_h5(
  f_path,
  version,
  indptr,
  indices,
  data,
  n_cells,
  n_features,
  barcodes,
  feature_ids,
  feature_names,
  feature_types
)
```

## Arguments

- f_path:

  String. Path of the file to create.

- version:

  String. One of `c("v3", "v2")`.

- indptr:

  Integer vector. Row pointers of the cells x features CSR matrix
  (0-indexed).

- indices:

  Integer vector. Feature indices (0-indexed).

- data:

  Integer vector. The counts. Must be non-negative.

- n_cells:

  Integer. Number of cells.

- n_features:

  Integer. Number of features.

- barcodes:

  Character vector. One barcode per cell.

- feature_ids:

  Character vector. One id per feature.

- feature_names:

  Character vector. One name per feature.

- feature_types:

  Optional character vector. One type per feature, v3 only. `NULL`
  writes `"Gene Expression"` for every feature.

## Value

Invisible `NULL`.
